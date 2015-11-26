package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import fastFunctions.StringTools.StringToolsFastHelper
import fastFunctions.VecTools.VecToolsFastHelper
import seqUtils.HRunHelper._

import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object QualitiesExtractor {

  def main(args: Array[String]) {
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(1)), 256 * 1024))

    var i = 0
    def proceedLine(line: String) = {
      val entries = line.split("\t")
      val quality = entries(10).fastMap(qualityToProb).mkStringFast('\n', '\n')
      writer.write(quality)
      if (i % 1000 == 0)
        println(s"Proceed $i lines")
      i += 1
      writer.flush()
    }

    val startTime = System.currentTimeMillis()
    Source.fromFile(args(0)).getLines().drop(4).foreach(proceedLine)
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }
}