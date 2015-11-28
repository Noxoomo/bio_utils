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

  def apply(read: HRunSeq, ref: HRunSeq): Int

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
    val result = (if (genom.contains(read)) {
      1
    } else {
      -1
    }) * base.apply(read, ref)
    result
  }

  override def toCsvLine(key: Int, value: Long, sep: Char = '\t'): String = {
    val distance = FastMath.abs(key)
    val genomic = key >= 0
    f"${name()}$sep$genomic$sep$distance$sep$value"
  }
}

class CountStatistics {
  private val stats = new TIntLongHashMap()

  @inline
  def addStat(stat: Int, readQuality: Array[Double], inc: Long = 1): Unit = {
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


trait StatsHolder extends (Array[Double] => CountStatistics) {
  def foreach(func: (CountStatistics, String) => Unit)
}

class SimpleStatHolder extends StatsHolder {
  val stat = new CountStatistics

  override def apply(quality: Array[Double]): CountStatistics = stat

  override def foreach(func: (CountStatistics, String) => Unit): Unit = func(stat, "")
}

object SimpleStatHolder {
  def apply(u: Unit): StatsHolder = new SimpleStatHolder
}

class SimpleStatHolderFactory extends ((Unit) => StatsHolder) {
  override def apply(v1: Unit): StatsHolder = SimpleStatHolder()
}

class SimpleQualityHolder extends StatsHolder {
  val stats = {
    val res = Array.ofDim[CountStatistics](101)
    for (i <- res.indices) {
      res(i) = new CountStatistics
    }
    res
  }

  def quality(qual: Array[Double]): Int = {
    (FastMath.exp(VecTools.sum(qual)) * 100).toInt
  }

  override def foreach(func: (CountStatistics, String) => Unit): Unit = {
    for (i <- 0 to 100) {
      val id = i * 0.01
      func(stats(i), f"$id.2")
    }
  }

  override def apply(qual: Array[Double]): CountStatistics = stats(quality(qual))
}

class SimpleQualityStatHolderFactory extends ((Unit) => StatsHolder) {
  override def apply(v1: Unit): StatsHolder = new SimpleQualityHolder
}

class HKMerDistanceStat(val metrics: Seq[HKMerDistance], val hkmerSize: Int = 16, val statHolderFactory: ((Unit) => StatsHolder) = new SimpleQualityStatHolderFactory) {
  private val refCache = Array.ofDim[HRun](hkmerSize)
  private val readCache = Array.ofDim[HRun](hkmerSize)
  private val qualityCache = Array.ofDim[Double](hkmerSize)
  private val metricStats = {
    val data = Array.ofDim[StatsHolder](metrics.size)
    for (i <- data.indices) {
      data(i) = statHolderFactory()
    }
    data
  }

  @inline
  def fill[@specialized(Int) T](buffer: Array[T], from: Array[T], offset: Int, size: Int): Unit = {
    var i = 0
    while (i < size) {
      buffer(i) = from(offset + i)
      i += 1
    }
  }


  @inline
  private def calcMetrics(ref: HRunSeq, read: HRunSeq, readQuality: Array[Double]): Unit = {
    metrics.zip(metricStats).foreach({ case (metric, statHolder) => statHolder(readQuality).addStat(metric(read, ref), readQuality) })
  }


  def proceedAlignedRead(ref: HRunSeq, read: HRunSeq, readQuality: Array[Double] = null): Unit = {
    for (offset <- 0 until (read.length - hkmerSize)) {
      fill(refCache, ref, offset, hkmerSize)
      fill(readCache, read, offset, hkmerSize)
      if (readQuality != null) {
        fill(qualityCache, readQuality, offset, hkmerSize)
      }
      calcMetrics(refCache, readCache, if (readQuality != null) qualityCache else null)
    }
  }


  def toCsv(sep: Char = '\t', prefix: String = ""): String = {
    val builder = new StringBuilder
    @inline
    def proceedMetric(prefix: String, metric: HKMerDistance, stat: CountStatistics): Unit = {
      stat.foreach({ case (key, count) => {
        builder.append(f"$prefix$sep${metric.toCsvLine(key, count, sep)}\n")
      }
      })
    }
    metrics.zip(metricStats).foreach({ case (metric, stats) => stats.foreach({ case (stat, id) => proceedMetric(f"$prefix$sep${id.toString}", metric, stat) }) })
    builder.mkString
  }
}
