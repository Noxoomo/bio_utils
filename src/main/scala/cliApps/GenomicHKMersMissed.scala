package cliApps

import seqUtils._

/**
  * User: Noxoomo
  * Date: 18.10.15
  * Time: 20:19
  */


object GenomicHKMersMissed {

  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")
    val centersSet = new SimpleTrie("", chromosome.length)
    val startTime = System.currentTimeMillis()

    Source.fromFile(args(1)).getLines().grouped(2).foreach(entry => {
      val quality = (entry.head.split("_")(2).toDouble * 100000).toInt * 0.00001
      if (quality > 0.0) {
        centersSet.+(entry(1))
        centersSet.+(entry(1).reverseComplementary())
      }
    })



    val hrunChromosome = HRunHelper.stringToHRunSeq(chromosome).map(HRunHelper.mkString)

    @inline
    def makeSliceStr(offset: Int): String = {
      val builder = new StringBuilder
      var i = offset
      while (i < (offset + 16)) {
        builder.append(hrunChromosome(i))
        i += 1
      }
      builder.toString()
    }

    var missed = 0
    var total = 0
    for (i <- 0 until hrunChromosome.length - 16) {
      val hkmer = makeSliceStr(i)
      if (!centersSet.contains(hkmer)) {
        missed += 1
      }
      total += 1
    }

    println(f"missed hkmers: $missed; total hkmers: $total")
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
  }
}