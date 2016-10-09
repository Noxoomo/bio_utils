package cliApps

import seqUtils._

/**
  * User: Noxoomo
  * Date: 18.10.15
  * Time: 20:19
  */


object ClusterCentersQuality {

  def main(args: Array[String]) {
    val chromosome = Source.fromFile(args(0)).getLines().drop(1).mkString("")

    //    val centersSet = new SimpleTrie("", chromosome.length)
    val centerKmers = new mutable.HashSet[String]()
    val refKmers = new mutable.HashSet[String]()

    val startTime = System.currentTimeMillis()

    val border = args(2).toDouble

    Source.fromFile(args(1)).getLines().grouped(2).foreach(entry => {
      val quality = 1 - entry.head.split("_")(2).toDouble
      if (quality <= border) {
        centerKmers.add(entry(1))
        centerKmers.add(entry(1).reverseComplementary())
      }
    })

    val hrunChromosome = HRunHelper.stringToHRunSeq(chromosome).map(HRunHelper.mkString)

    @inline
    def nextHrunSize(offset: Int): Int = {
      var i = offset + 1
      while (i < hrunChromosome.length && hrunChromosome(i) == hrunChromosome(i - 1)) {
        i += 1
      }
      i - offset
    }

    @inline
    def nextHKMer(offset: Int): String = {
      val builder = new StringBuilder
      var i = offset
      var k = 0
      while (i < hrunChromosome.length && k < 16) {
        var sz = nextHrunSize(i)
        val nextStart = i + sz
        while (sz > 0) {
          builder.append(hrunChromosome(i))
          sz = sz - 1
        }
        k = k + 1
        i = nextStart
      }
      if (k == 16) {
        builder.toString()
      } else {
        ""
      }
    }


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

    for (i <- hrunChromosome.indices) {
      val hkmer = nextHKMer(i)
      if (hkmer != "") {
        refKmers.add(hkmer)
        refKmers.add(hkmer.reverseComplementary())
      }
    }

    for (kmer <- refKmers) {
      if (!centerKmers.contains(kmer)) {
        missed += 1
      }
    }

    var wrongCenters = 0
    for (kmer <- centerKmers) {
      if (!refKmers.contains(kmer)) {
        wrongCenters += 1
      }
    }
    val totalGenom = refKmers.size
    val totalCenter = centerKmers.size
    val genomicPart = 1.0 * (totalCenter - wrongCenters) / totalGenom
    val missedPart = missed * 1.0 / totalGenom



    println(f"total genom hkmers: $totalGenom")
    println(f"total centers hkmers: $totalCenter")
    println(f"genomic hkmers: $genomicPart")
    println(f"wrong centers hkmers: $wrongCenters")
    println(f"missed hkmers: $missed")
    println(f"missed hkmers part: $missedPart")
    val endTime = System.currentTimeMillis()
    println(f"Working time: ${endTime - startTime}")
  }
}