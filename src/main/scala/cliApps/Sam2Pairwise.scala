package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import seqUtils.Cigar
import seqUtils.Cigar._

import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object Sam2Pairwise {


  def main(args: Array[String]) {

    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))

    def proceedRead(query: String, read: String, ref: String, cigar: Cigar, flag: Int): Unit = {
      writer.write(query)
      writer.write("\t")
      writer.write(flag.toString)
      writer.write("\t")
      writer.write(cigar.toString())
      writer.write("\n")
      val result = applyCigar(cigar, read, ref)
      writer.write(result._1)
      writer.write("\n")
      writer.write(result._3)
      writer.write("\n")
      writer.write(result._2)
      writer.write("\n")
      writer.flush()
    }

    var i = 0

    def proceedLine(line: String) = {
      val entries = line.split("\t")
      val query = entries(0)
      val refStart = Integer.parseInt(entries(3)) - 1
      val flag = Integer.parseInt(entries(1))
      if ((flag & 4) == 0) {
        //4 is unmapped
        val cigar = Cigar.parseCigar(entries(5))
        val read = entries(9)
        val refSize = cigar.matchCount + cigar.deletionCount
        val refEnd = refStart + refSize
        proceedRead(query, read, chromosome.substring(refStart, refEnd), cigar, flag)
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