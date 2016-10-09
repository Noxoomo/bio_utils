package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import org.apache.commons.math3.util.FastMath
import seqUtils.Cigar
import seqUtils.Cigar._
import seqUtils.HRunHelper._
import fastFunctions.StringTools._

import scala.annotation.tailrec
import scala.io.Source

/**
  * User: Noxoomo
  * Date: 24.05.16
  * Time: 23:06
  */
object Comparator {

  class ReadStat(val id: String) {
    var insertions = 0
    var deletions = 0
    var mismatch = 0
    var distance = 0
    var hammingDistance = 0
  }

  @inline
  final def calcStats(query: String, cigar: Cigar, reference: HRunSeq, read: HRunSeq): ReadStat = {
    val readStat = new ReadStat(query)
    @inline
    def addHRuns(ref: HRun, read: HRun): Unit = {
      val errorSize = FastMath.abs(size(ref) - size(read))

      if (ref != read) {
        readStat.hammingDistance += 1
        readStat.distance += errorSize
        if (errorSize == 0) {
          readStat.distance += 1
        }

        if (size(read) > size(ref)) {
          readStat.deletions += 1
        } else if (size(read) < size(ref)) {
          readStat.insertions += 1
        } else {
          readStat.mismatch += 1
        }
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
            addHRuns(reference(refOffset + i), read(readOffset + i))
            i += 1
          }
          helper(cigarOffset + 1, refOffset + i, readOffset + i)
        } else if (oprType == 'I') {
          helper(cigarOffset + 1, refOffset + oprSize, readOffset)
        } else {
          var i = 0
          while (i < oprSize) {
            val fakeRef = hrun(base(read(readOffset + i)), 0)
            addHRuns(fakeRef, read(readOffset + i))
            i += 1
          }
          helper(cigarOffset + 1, refOffset, readOffset + oprSize)
        }
      }
    }
    helper(0, 0, 0)
    readStat
  }


  class Writer(filename: String) {
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(filename), 64 * 1024))
    var i = 0
    writer.write("id\thammingDistance\tlevDistance\tinsertions\tdeletions\tmismatch\n")

    def writeStat(stat: ReadStat) {
      val line = f"${stat.id}\t${stat.hammingDistance}\t${stat.distance}\t${stat.insertions}\t${stat.deletions}\t${stat.mismatch}\n"
      synchronized {
        writer.write(line)
        if (i % 100000 == 0)
          println(s"Write $i lines")
        i += 1
      }
    }
  }


  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val writer = new Writer(args(2))

    def proceedRead(query: String, read: String, quality: String, ref: String, cigar: Cigar): Unit = {
      val inHRuns = applyAndConvertToHRuns(cigar.inverse(), ref, read)
      val readHRun = inHRuns._2
      val stats = calcStats(query, inHRuns._3, inHRuns._1, readHRun)
      writer.writeStat(stats)
    }

    var numThreads = 4
    val blockSize = 16384

    if (args.length > 3) {
      numThreads = Integer.valueOf(args(3))
      print("Max threads count: " + numThreads)
    }

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

        proceedRead(query, read, quality, ref, cigar)
      }

    }

    val startTime = System.currentTimeMillis()
    Source.fromFile(args(1)).getLines().filterNot(_.charAt(0) == '@').grouped(numThreads * blockSize).foreach(_.grouped(blockSize).toSeq.par.foreach(_.foreach(proceedLine)))
    val endTime = System.currentTimeMillis()
    writer.writer.flush()
    writer.writer.close()
    println(f"Working time: ${endTime - startTime}")
  }
}
