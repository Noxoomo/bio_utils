package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import gnu.trove.map.hash.TDoubleLongHashMap
import gnu.trove.procedure.TDoubleLongProcedure
import seqUtils._
import fastFunctions.StringTools._
import scala.io.Source

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object GenomicClusterCenters {
  val header = f"quality\tgenomic\ttotal\n"
  val stats = new TDoubleLongHashMap()

  @inline
  def add(quality: Double, genomic: Boolean): Unit = {
    stats.put(quality, stats.get(quality) + (1L << 32 | (if (genomic) 1 else 0)))
  }

  @inline
  def toCsvLine(quality: Double, stat: Long): String = {
    f"$quality\t${stat & ((1L << 32) - 1)}\t${stat >> 32}\n"
  }

  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val kmersSet = new SimpleTrie(chromosome, 100)
    val kmersComplSet = new SimpleTrie(chromosome.reverseComplementary(), 100)
    val startTime = System.currentTimeMillis()
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))
    writer.write(header)
    val discretizationLevel = args(3).toInt
    Source.fromFile(args(1)).getLines().grouped(2).foreach(entry => {
      val quality = (entry.head.split("_")(2).toDouble * discretizationLevel).toInt * 1.0 / discretizationLevel
      add(quality, kmersSet.contains(entry(1)) || kmersComplSet.contains(entry(1)))
    })
    stats.forEachEntry(new TDoubleLongProcedure {
      override def execute(quality: Double, stats: Long): Boolean = {
        writer.write(toCsvLine(quality, stats))
        true
      }
    })
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }
}