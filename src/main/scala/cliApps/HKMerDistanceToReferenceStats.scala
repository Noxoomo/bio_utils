package cliApps

import java.io.{FileOutputStream, OutputStreamWriter}
import java.util.zip.GZIPOutputStream

import seqUtils.HRunHelper.{HRun, HRunSeq, HRunSeqMeta, _}
import seqUtils._

/**
 * User: Noxoomo
 * Date: 18.10.15
 * Time: 20:19
 */


object HKMerDistanceToReferenceStats {

  def wrap(dist: HKMerDistance, genomHKMers: Set[HRunSeq]): HKMerDistance = {
    new HRunDistanceWithGenomicFlag(dist, genomHKMers)
  }

  def createMetrics(genomHKMers: Set[HRunSeq], size: Int): HKMerDistanceStat = {
    new HKMerDistanceStat(List(wrap(new HRunHamming, genomHKMers), wrap(new HRunLevenstein, genomHKMers)), hkmerSize = size)
  }

  val header = f"hkmerSize\tgoodProb\tmetric\tgenomic\tvalue\tcount\n"
  val sizes = Array(8, 12, 16, 20)


  def main(args: Array[String]) {
    val aligner = new SeqAligner[HRun](new HRunPenalty())
    val proceeder = new HRunReadsToReferenceFromSamProceeder(args(0), args(1))
    val kmersSet = new RadixTreeHRunSet(proceeder.chromosome, sizes.max * 2)
    val hkmerStatBuilders = sizes.map(createMetrics(kmersSet, _))

    @inline
    def proceedRead(query: String, flag: Int, cigar: Cigar, read: HRunSeq, readMeta: HRunSeqMeta, ref: HRunSeq): Unit = {
      val (alignedRef, alignedRead) = align(cigar, ref, read)
      val alignedQualities = HRunSeqQuality.alignQuality(alignedRead, readMeta.qualities)
      hkmerStatBuilders.par.foreach(_.proceedAlignedRead(alignedRef, alignedRead, alignedQualities))
    }
    val startTime = System.currentTimeMillis()
    proceeder.proceed(proceedRead)
    val endTime = System.currentTimeMillis()
    val writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args(2)), 64 * 1024))
    writer.write(header)
    hkmerStatBuilders.zip(sizes).foreach({ case (builder, size) => writer.write(builder.toCsv(prefix = f"$size")) })
    println(f"Working time: ${endTime - startTime}")
    writer.flush()
    writer.close()
  }
}