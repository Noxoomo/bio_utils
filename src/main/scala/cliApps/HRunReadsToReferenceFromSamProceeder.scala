package cliApps

import fastFunctions.StringTools.CharSeqHelper
import seqUtils.Cigar._
import seqUtils.HRunHelper._
import seqUtils._

import scala.io.Source

/**
 * User: Noxoomo
 * Date: 21.11.15
 * Time: 16:25
 */

class HRunReadsToReferenceFromSamProceeder(val referenceFilename: String, val readFilename: String) {
  val chromosome = Source.fromFile(referenceFilename).getLines().drop(1).mkString("")
  val chromosomeReverse = chromosome.reverse

  def proceed(func: (String, Int, Cigar, HRunSeq, HRunSeqMeta, HRunSeq) => Unit) {
    var i = 0
    def proceedLine(line: String) = {
      val entries = line.split("\t")
      val query = entries(0)
      val offset = Integer.parseInt(entries(3)) - 1
      val flag = Integer.parseInt(entries(1))
      if ((flag & 4) == 0) {
        val compl = (flag & 16) != 0
        val cigarStr = Cigar.parseCigar(entries(5), compl).withoutClipping()
        val refSize = cigarStr.matchCount + cigarStr.deletionCount

        val read = if (!compl) entries(9) else entries(9).reverse
        val quality = if (!compl) entries(10) else entries(10).reverse

        val refEnd = offset + refSize
        val refStart = offset
        val ref = if (!compl) chromosome.viewSlice(refStart, refEnd) else chromosome.viewSlice(refStart, refEnd).reverse

        val (refHRun, readHRun, cigar) = applyAndConvertToHRuns(cigarStr.inverse(), ref, read)
        val readMeta = readHRun.meta(0, quality)

        func(query, flag, cigar, readHRun, readMeta, refHRun)
      }
      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
    }
    Source.fromFile(readFilename).getLines().drop(4).foreach(proceedLine)
  }
}
