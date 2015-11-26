package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import fastFunctions.StringTools.StringToolsFastHelper
import fastFunctions.VecTools.VecToolsFastHelper
import org.apache.commons.math3.util.FastMath
import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqHelper, HRunSeqMeta, _}
import seqUtils.{Cigar, HRunHelper, HRunPenalty, SeqAligner}

import scala.annotation.tailrec
import scala.io.Source
import scala.util.Random

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object LearnTestQualitiesDatasets {

  def main(args: Array[String]) {

    val hrunSize = Integer.parseInt(args(2))
    val learnWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(f"${args(3)}_${hrunSize}_learn.gz"), 64 * 1024))
    val testWriter = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(f"${args(3)}_${hrunSize}_test.gz"), 64 * 1024))
    val learnFraction = 0.1
    val random = new Random()

    @inline
    def write(lines: String) = {
      val writer = if (random.nextDouble() < learnFraction) {
        learnWriter
      } else {
        testWriter
      }
      writer.write(lines)
      writer.flush()
    }

    @inline
    def errorProb(qualities: Array[Double]): Double = {
      var noErrorProb = 0.0
      for (quality <- qualities) {
        noErrorProb += FastMath.log(1 - quality)
      }
      1 - FastMath.exp(noErrorProb)
    }

    @inline
    def calcStats(query: String, flag: Int, cigar: Cigar, reference: HRunSeq, read: HRunSeq, readMeta: HRunSeqMeta): String = {
      val builder = new StringBuilder()
      @inline
      def addHRuns(ref: HRun, read: HRun, metaIdx: Int): Unit = {
        if (size(read) == hrunSize) {
          val errorType = if (base(ref) == base(read)) {
            if (size(ref) == size(read)) {
              "E"
            } else if (size(ref) > size(read)) "D" else "I"
          } else {
            "M"
          }

          val qualities = readMeta.qualities(metaIdx).fastMap(qualityToProb)
          builder.
            append(errorType).append("\t").
            append(readMeta.offsets(metaIdx)).append("\t").
            append(base(read)).append("\t").
            append(errorProb(qualities)).append("\t").
            append(qualities.mkStringFast('\t', '\n'))
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
      builder.toString()
    }



    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val aligner = new SeqAligner[HRun](new HRunPenalty())

    def proceedRead(query: String, flag: Int, read: HRunSeq, readMeta: HRunSeqMeta, ref: HRunSeq): Unit = {
      val cigar = aligner(ref, read)
      write(calcStats(query, flag, cigar, ref, read, readMeta))
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
    learnWriter.close()
    testWriter.close()
  }
}