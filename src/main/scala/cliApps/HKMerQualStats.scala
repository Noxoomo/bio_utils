package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import fastFunctions.StringTools.StringToolsFastHelper
import fastFunctions.VecTools.VecToolsFastHelper
import gnu.trove.map.hash.TLongLongHashMap
import gnu.trove.procedure.TLongLongProcedure
import org.apache.commons.math3.util.FastMath
import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqHelper, HRunSeqMeta, _}
import seqUtils._

import scala.io.Source
import scala.reflect.ClassTag

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object HKMerQualStats {
  val key = "start_offset\tquality"
  //10 bit + 10 bit + 2 bit + 10 bit + 10 bit
  val value = "errors\ttotal"
  val header = f"$key\t$value\n"

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

  val qualityAccuracy = 1e3

  @inline
  def compressKey(offset: Int, quality: Double): Long = {
    0L | (offset.toLong | ((quality * qualityAccuracy).toLong << 10))
  }

  @inline
  def keyToString(key: Long): String = {
    val mask = 1023
    val offset = key & mask
    val quality = (key >> 10) * 1.0 / qualityAccuracy
    f"$offset\t$quality"
  }

  class StatBuilder {
    val stats = new TLongLongHashMap()

    @inline
    def addHRMer(offset: Int, quality: Double, error: Int): Unit = {
      val key = compressKey(offset, quality)
      val count = (1L << 32 | error) + stats.get(key)
      stats.put(key, count)
    }

    override def toString(): String = {
      val builder = new StringBuilder()
      val mask = (1L << 32) - 1
      stats.forEachEntry(new TLongLongProcedure {
        override def execute(key: Long, value: Long): Boolean = {
          builder.append(keyToString(key)).append('\t').
            append(value & mask).append("\t").
            append(value >> 32).append("\n")
          true
        }
      }
      )
      builder.toString()
    }
  }

  val statBuilder = new StatBuilder()


  class WindowStat[@specialized T <: AnyVal](windowsStat: Int)(implicit m: ClassTag[T]) {
    val data = Array.ofDim[T](windowsStat)
    private var cursor = 0

    @inline
    def add(value: T): T = {
      val removed = data(cursor % windowsStat)
      data(cursor % windowsStat) = value
      cursor += 1
      removed
    }

    @inline
    def tail(): T = data((cursor - windowsStat) % windowsStat)

    @inline
    def head(): T = data(cursor % windowsStat)

    @inline
    def updateHead(newValue: T): Unit = data(cursor % windowsStat) = newValue

    @inline
    def filled() = cursor > windowsStat
  }

  class HKMerStat(val length: Int = 16) {
    private val quality = new WindowStat[Double](length)
    private val offsets = new WindowStat[Int](length)
    private val errorsCount = new WindowStat[Int](length)
    private var totalQuality = 0.0
    private var errors = 0

    @inline
    def addErrorType(good: Boolean, updateHead: Boolean = false): Unit = {
      errors += (if (good) 0 else 1)
      if (updateHead) {
        val headSz = errorsCount.head() + (if (good) 0 else 1)
        errorsCount.updateHead(headSz)
      } else {
        errors -= errorsCount.add(if (good) 0 else 1)
      }
    }

    @inline
    def addHRun(offset: Int, qualities: Array[Double]): Unit = {
      val currentQuality = qualities.fastMap(q => FastMath.log(1 - q)).fastSum()
      offsets.add(offset)
      totalQuality += currentQuality
      totalQuality -= quality.add(currentQuality)
      if (quality.filled()) {
        statBuilder.addHRMer(offsets.head(), 1 - FastMath.exp(totalQuality), if (errors > 0) 1 else 0)
      }
    }
  }

  val hkmerBuilder = new HKMerStat(16)


  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))

    writer.write(header)
    val aligner = new SeqAligner[HRun](new HRunPenalty())


    def proceedRead(query: String, flag: Int, read: HRunSeq, readMeta: HRunSeqMeta, ref: HRunSeq): Unit = {
      val cigar = aligner(ref, read)

      @inline
      def addHRuns(ref: HRun, read: HRun, readHRunIdx: Int): Unit = {
        val good = base(ref) == base(read) && size(ref) == size(read)
        hkmerBuilder.addErrorType(good, size(read) == 0)
        if (size(read) != 0) {
          val qualities = readMeta.qualities(readHRunIdx)
          val offset = readMeta.offsets(readHRunIdx)
          hkmerBuilder.addHRun(offset, qualities.fastMap(qualityToProb))
        }
      }

      val proceeder = new CigarProceeder(addHRuns)
      proceeder.proceed(cigar, ref, read)
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
        val read = entries(9)
        val refSize = cigar.matchCount + cigar.deletionCount
        val clippingSize = cigar.clippingSize

        val readEnd = if (compl) read.length else read.length - clippingSize
        val readStart = if (compl) clippingSize else 0

        val readHRunSeq = HRunHelper.stringToHRunSeq(read, readStart, readEnd, compl)
        val readMeta = readHRunSeq.meta(0, entries(10), compl)

        val refEnd = offset + refSize
        val refStart = offset
        val refHRunSeq = HRunHelper.stringToHRunSeq(chromosome, refStart, refEnd, compl)

        proceedRead(query, flag, readHRunSeq, readMeta, refHRunSeq)
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