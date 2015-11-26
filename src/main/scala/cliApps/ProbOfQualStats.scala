package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqHelper, HRunSeqMeta, _}
import seqUtils.{Cigar, HRunHelper, HRunPenalty, SeqAligner}

import scala.annotation.tailrec
import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object ProbOfQualStats {
  val header = "quality\terror\ttotal\tprob\n"

  val errorStats = Array.ofDim[Int](256)
  val totalStats = Array.ofDim[Int](256)


  @inline
  final def calcStats(query: String, flag: Int, cigar: Cigar, reference: HRunSeq, read: HRunSeq, readMeta: HRunSeqMeta) = {
    @inline
    def addHRuns(ref: HRun, read: HRun, metaIdx: Int): Unit = {
      for (quality <- readMeta.qualities(metaIdx)) {
        totalStats(quality.toInt) += 1
      }
      if (base(ref) == base(read)) {
        val qualities = readMeta.qualities(metaIdx)
        if (size(ref) > size(read)) {
          errorStats(qualities(size(read) - 1).toInt) += 1
        } else {
          for (i <- Range(size(ref), size(read))) {
            errorStats(qualities(i).toInt) += 1
          }
        }
      } else {
        for (quality <- readMeta.qualities(metaIdx)) {
          errorStats(quality.toInt) += 1
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
    val aligner = new SeqAligner[HRun](new HRunPenalty())

    def proceedRead(query: String, flag: Int, read: HRunSeq, readMeta: HRunSeqMeta, ref: HRunSeq): Unit = {
      val cigar = aligner(ref, read)
      calcStats(query, flag, cigar, ref, read, readMeta)
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
        val clippingSize = cigar.clippingSize
        val refSize = cigar.matchCount + cigar.deletionCount

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
    println(f"Working time: ${endTime - startTime}")
    for (i <- totalStats.indices) {
      if (totalStats(i) > 0) {
        val total = totalStats(i)
        val error = errorStats(i)
        val prob = qualityToProb(i.toChar)
        writer.write(f"$prob\t$error\t$total\t${error * 1.0 / total}\n")
      }
    }
    writer.flush()
    writer.close()
  }
}