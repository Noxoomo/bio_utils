package seqUtils

import fastFunctions.VecTools
import gnu.trove.map.hash.TIntLongHashMap
import gnu.trove.procedure.TIntLongProcedure
import org.apache.commons.math3.util.FastMath
import seqUtils.HRunHelper.{HRun, HRunSeq, _}

/**
 * User: Noxoomo
 * Date: 21.11.15
 * Time: 14:16
 */


trait HKMerDistance extends ((HRunSeq, HRunSeq) => Int) {
  def name(): String

  def toCsvLine(key: Int, value: Long, sep: Char = '\t'): String = f"${name()}$sep$key$sep$value"
}

class HRunHamming extends HKMerDistance {
  override def name(): String = "Hamming"

  override def apply(read: HRunSeq, ref: HRunSeq): Int = {
    VecTools.zipAggregate(read, ref, (errors: Int, readHRun: HRun, refHRun: HRun) => {
      if (readHRun != refHRun) {
        errors + 1
      } else errors
    })
  }
}

class HRunLevenstein extends HKMerDistance {

  override def name(): String = "Levenshtein"

  override def apply(read: HRunSeq, ref: HRunSeq): Int = {
    assert(read.length == ref.length)
    var result = 0
    @inline
    def errors(count: Int, readHRun: HRun, refHRun: HRun): Int = {
      count + (if (base(readHRun) == base(refHRun)) {
        FastMath.abs(size(readHRun) - size(refHRun))
      } else {
        FastMath.max(size(readHRun), size(refHRun))
      })
    }
    VecTools.zipAggregate(read, ref, errors)
  }
}


class HRunDistanceWithGenomicFlag(base: HKMerDistance, genom: Set[HRunSeq]) extends HKMerDistance {
  override def name(): String = f"genomic_${base.name()}"

  override def apply(read: HRunSeq, ref: HRunSeq): Int = {
    if (genom.contains(read)) {
      1
    } else {
      0
    } << 31 | base.apply(read, ref)
  }

  override def toCsvLine(key: Int, value: Long, sep: Char = '\t'): String = {
    val distance = ((1 << 31) - 1) & key
    val genomic = key & 1 << 31
    f"${name()}$sep$genomic$sep$distance$sep$value"
  }
}

class LongAdditiveStatistic {
  private val stats = new TIntLongHashMap()

  @inline
  def addStat(stat: Int, inc: Long = 1): Unit = {
    val count = inc + stats.get(stat)
    stats.put(stat, count)
  }

  def foreach(func: (Int, Long) => Unit): Unit = {
    stats.forEachEntry(new TIntLongProcedure {
      override def execute(key: Int, value: Long): Boolean = {
        func(key, value)
        true
      }
    })
  }
}

class HKMerDistanceStat(val metrics: Seq[HKMerDistance], val hkmerSize: Int = 16) {
  private val refCache = Array.ofDim[HRun](hkmerSize)
  private val readCache = Array.ofDim[HRun](hkmerSize)
  private val stats = {
    val data = Array.ofDim[LongAdditiveStatistic](metrics.size)
    for (i <- data.indices) {
      data(i) = new LongAdditiveStatistic
    }
    data
  }

  @inline
  def fill(buffer: HRunSeq, from: HRunSeq, offset: Int, size: Int): Unit = {
    var i = 0
    while (i < size) {
      buffer(i) = from(offset + i)
      i += 1
    }
  }

  @inline
  private def calcMetrics(ref: HRunSeq, read: HRunSeq): Unit = {
    metrics.zip(stats.indices).foreach({ case (metric, id) => stats(id).addStat(metric(read, ref)) })
  }

  def proceedAlignedRead(ref: HRunSeq, read: HRunSeq): Unit = {
    for (offset <- 0 until (read.length - hkmerSize)) {
      fill(refCache, ref, offset, hkmerSize)
      fill(readCache, read, offset, hkmerSize)
      calcMetrics(refCache, readCache)
    }
  }


  def toCsv(sep: Char = '\t', prefix: String = ""): String = {
    val builder = new StringBuilder
    @inline
    def proceedMetric(metric: HKMerDistance, stats: LongAdditiveStatistic): Unit = {
      stats.foreach({ case (key, count) => {
        builder.append(f"$prefix$sep${metric.toCsvLine(key, count, sep)}\n")
      }
      })
    }
    metrics.zip(stats).foreach({ case (metric, stat) => proceedMetric(metric, stat) })
    builder.mkString
  }
}
