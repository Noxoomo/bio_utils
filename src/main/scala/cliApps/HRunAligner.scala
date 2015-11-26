import java.io.{BufferedWriter, FileWriter}

import seqUtils.HRunHelper.{HRun, HRunSeq}
import seqUtils.{Cigar, HRunHelper, HRunPenalty, SeqAligner}

import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object HRunAligner {


  def main(args: Array[String]) {

    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val writer = new BufferedWriter(new FileWriter(args(2))) // new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))


    val aligner = new SeqAligner[HRun](new HRunPenalty())

    def proceedRead(query: String, read: HRunSeq, ref: HRunSeq, flag: Int): Unit = {
      val cigar = aligner(ref, read)
      writer.write(query)
      writer.write("\t")
      writer.write(flag.toString)
      writer.write("\t")
      writer.write(cigar.toString())
      writer.write("\n")
      val result = HRunHelper.mkString(cigar, ref, read)
      writer.write(result._1)
      writer.write("\n")
      writer.write(result._2)
      writer.write("\n")
      writer.flush()
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
        val refEnd = offset + refSize
        val refStart = offset
        proceedRead(query, HRunHelper.stringToHRunSeq(read, readStart, readEnd, compl), HRunHelper.stringToHRunSeq(chromosome, refStart, refEnd, compl), flag)
      }

      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
    }
    val startTime = System.currentTimeMillis()
    Source.fromFile(args(1)).getLines().drop(4).foreach(proceedLine)
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }
}