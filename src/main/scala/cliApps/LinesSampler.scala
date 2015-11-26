package cliApps

import java.io._
import java.util.zip.{GZIPInputStream, GZIPOutputStream}

import scala.io.Source
import scala.util.Random

/**
 * User: Noxoomo
 * Date: 26.10.15
 * Time: 1:20
 */
object LinesSampler {

  def main(args: Array[String]) {
    val takenProb = args(0).toDouble
    val random = new Random()
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))
    @inline
    def proceedLine(line: String) = {
      if (random.nextDouble() < takenProb) {
        writer.write(line)
        writer.write("\n")
      }
    }

    val startTime = System.currentTimeMillis()
    Source.fromInputStream(new GZIPInputStream(new FileInputStream(args(1)))).getLines().drop(4).foreach(proceedLine)
    val endTime = System.currentTimeMillis()

    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }

}
