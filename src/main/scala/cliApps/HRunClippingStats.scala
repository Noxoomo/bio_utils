package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqHelper, HRunSeqMeta, _}
import seqUtils.{Cigar, HRunHelper}

import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object HRunClippingStats {

  def main(args: Array[String]) {

    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(1)), 64 * 1024))
    writer.write("qid\toffset\tbase\tsize\tqualities\n")

    def proceedRead(query: String, read: HRunSeq, readMeta: HRunSeqMeta): Unit = {
      val statsBuilder = new StringBuilder()
      @inline
      def addHRuns(read: HRun, quality: String, offset: Int): Unit = {
        statsBuilder.append(query).append("\t").append(offset).append("\t").
          append(base(read)).append("\t").
          append(size(read)).append("\t").
          append(quality).append("\n")
      }
      var i = 0
      while (i < read.length) {
        addHRuns(read(i), readMeta.qualities(i), readMeta.offsets(i))
        i += 1
      }
      writer.write(statsBuilder.toString())
      writer.flush()
    }

    var i = 0

    def proceedLine(line: String) = {
      val entries = line.split("\t")
      val flag = Integer.parseInt(entries(1))
      if ((flag & 4) == 0) {
        //4 is unmapped
        val compl = (flag & 16) != 0
        val cigar = Cigar.parseCigar(entries(5), compl)
        val read = entries(9)
        val clippingSize = cigar.clippingSize
        if (clippingSize > 0) {
          val readEnd = if (compl) clippingSize else read.length
          val readStart = if (compl) 0 else read.length - clippingSize
          val readHRunSeq = HRunHelper.stringToHRunSeq(read, readStart, readEnd, compl)
          val readMeta = readHRunSeq.meta(read.length - clippingSize, entries(10), compl)
          proceedRead(entries(0), readHRunSeq, readMeta)
        }
      }
      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
    }

    val startTime = System.currentTimeMillis()
    Source.fromFile(args(0)).getLines().drop(4).foreach(proceedLine)
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }

}