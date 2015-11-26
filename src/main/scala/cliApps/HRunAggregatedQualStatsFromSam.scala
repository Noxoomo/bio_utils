package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import fastFunctions.StringTools.StringToolsFastHelper
import gnu.trove.map.hash.TLongLongHashMap
import gnu.trove.procedure.TLongLongProcedure
import org.apache.commons.math3.util.FastMath
import seqUtils.Cigar
import seqUtils.Cigar._
import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqHelper, HRunSeqMeta, _}

import scala.annotation.tailrec
import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object HRunAggregatedQualStatsFromSam {
  val key = "offset\tquality\thrun_size\toffset_in_hrun\thrun_error_type\thrun_error_size"
  //10 bit + 10 bit + 2 bit + 10 bit + 10 bit  + 10bit
  val value = "count"
  val header = f"$key\t$value\n"

  val errorTypeMap = Array('E', 'M', 'I', 'D', 'C')

  @inline
  def getErrorIdx(error: Char): Int = {
    error match {
      case 'E' => 0
      case 'M' => 1
      case 'I' => 2
      case 'D' => 3
      case 'C' => 4
    }
  }

  @inline
  def compressKey(offset: Int, hrunSize: Int, offsetInHRun: Int, quality: Char, errorType: Char, errorSize: Int): Long = {
    0L | (offset.toLong | (hrunSize.toLong << 10) | (offsetInHRun.toLong << 20) | (quality.toLong << 30) | (errorType.toLong << 40) | (errorSize.toLong << 50))
  }

  @inline
  def keyToString(key: Long): String = {
    val mask = 1023
    val offset = key & mask
    val hrunSize = (key >> 10) & mask
    val offsetInHRun = (key >> 20) & mask
    val quality = qualityToProb(((key >> 30) & mask).toChar)
    val errorType = ((key >> 40) & mask).toChar
    val errorSize = key >> 50
    f"$offset\t$quality\t$hrunSize\t$offsetInHRun\t$errorType\t$errorSize"
  }

  class StatBuilder {
    val stats = new TLongLongHashMap()

    @inline
    def add(offset: Int, hrunSize: Int, offsetInHRun: Int, quality: Char, error: Char, errorSize: Int): Unit = {
      val key = compressKey(offset, hrunSize, offsetInHRun, quality, error, errorSize)
      val count = 1 + stats.get(key)
      stats.put(key, count)
    }

    override def toString(): String = {
      val builder = new StringBuilder()
      stats.forEachEntry(new TLongLongProcedure {
        override def execute(key: Long, value: Long): Boolean = {
          builder.append(keyToString(key)).append('\t').append(value).append("\n")
          true
        }
      }
      )
      builder.toString()
    }
  }

  val statBuilder = new StatBuilder()

  @inline
  final def calcStats(query: String, flag: Int, cigar: Cigar, reference: HRunSeq, read: HRunSeq, readMeta: HRunSeqMeta) = {
    @inline
    def addHRuns(ref: HRun, read: HRun, metaIdx: Int): Unit = {
      val errorType = if (base(ref) == base(read)) {
        if (size(ref) == size(read)) {
          'E'
        } else if (size(ref) > size(read)) 'D' else 'I'
      } else {
        'M'
      }

      val hrunSize = size(read)
      val qualities = readMeta.qualities(metaIdx)
      val offset = readMeta.offsets(metaIdx)
      val errorSize = FastMath.abs(size(ref) - size(read))
      for (offsetInHRun <- qualities.indices) {
        statBuilder.add(offset, hrunSize, offsetInHRun, qualities(offsetInHRun), errorType, errorSize)
      }
    }

    @tailrec
    def helper(cigarOffset: Int, refOffset: Int, readOffset: Int): Unit = {
      if (cigarOffset < cigar.operations.length) {
        val opr = cigar.operations(cigarOffset)
        val oprType = Cigar.operationType(opr)
        val oprSize = Cigar.operationSize(opr)
        if (oprType == 'M') {
          var i = 0
          while (i < oprSize) {
            addHRuns(reference(refOffset + i), read(readOffset + i), readOffset + i)
            i += 1
          }
          helper(cigarOffset + 1, refOffset + i, readOffset + i)
        } else if (oprType == 'I') {
          helper(cigarOffset + 1, refOffset + oprSize, readOffset)
        } else {
          var i = 0
          while (i < oprSize) {
            val fakeRef = hrun(base(read(readOffset + i)), 0)
            addHRuns(fakeRef, read(readOffset + i), readOffset + i)
            i += 1
          }
          helper(cigarOffset + 1, refOffset, readOffset + oprSize)
        }
      }
    }
    helper(0, 0, 0)
  }


  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))
    writer.write(header)

    def proceedRead(query: String, flag: Int, read: String, quality: String, ref: String, cigar: Cigar): Unit = {
      val inHRuns = applyAndConvertToHRuns(cigar.inverse(), ref, read)
      val readHRun = inHRuns._2
      calcStats(query, flag, inHRuns._3, inHRuns._1, readHRun, readHRun.meta(0, quality))
    }

    var i = 0

    def proceedLine(line: String) = {
      val entries = line.split("\t")
      val query = entries(0)
      val offset = Integer.parseInt(entries(3)) - 1
      val flag = Integer.parseInt(entries(1))
      if ((flag & 4) == 0) {
        //4 is unmapped
        val compl = (flag & 16) != 0
        val cigar = Cigar.parseCigar(entries(5), compl)
        val refSize = cigar.matchCount + cigar.deletionCount

        val read = if (!compl) entries(9) else entries(9).reverse
        val quality = if (!compl) entries(10) else entries(10).reverse

        val refEnd = offset + refSize
        val refStart = offset
        val ref = if (!compl) chromosome.substring(refStart, refEnd) else chromosome.rsubstring(refStart, refEnd)

        proceedRead(query, flag, read, quality, ref, cigar)
      }
      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
    }
    val startTime = System.currentTimeMillis()
    Source.fromFile(args(1)).getLines().drop(4).foreach(proceedLine)
    val endTime = System.currentTimeMillis()
    writer.write(statBuilder.toString())
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }
}