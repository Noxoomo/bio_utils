package cliApps

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

  def proceed(func: (String, Int, HRunSeq, HRunSeqMeta, HRunSeq) => Unit) {
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
        func(query, flag, readHRunSeq, readMeta, refHRunSeq)
      }
      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
    }
    Source.fromFile(readFilename).getLines().drop(4).foreach(proceedLine)
  }
}
