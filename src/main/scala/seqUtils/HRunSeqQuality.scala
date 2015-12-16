package seqUtils

import gnu.trove.list.array.TDoubleArrayList
import org.apache.commons.math3.util.FastMath
import seqUtils.HRunHelper.{HRunSeq, _}

import scala.annotation.tailrec


/**
 * User: Noxoomo
 * Date: 28.11.15
 * Time: 15:12
 */


object HRunSeqQuality {
  @inline
  def qualityToLogCorrectProb(quality: Char): Double = FastMath.log(1.0 - qualityToProb(quality)) // i don't like 10-base scale, it's unnatural :)


  @inline
  def alignQuality(read: HRunSeq, qualities: Array[String]): Array[Double] = {
    val quality = new TDoubleArrayList()
    @tailrec
    def helper(offset: Int, qualityOffset: Int): Unit = {
      if (offset < read.length) {
        val len = size(read(offset))
        var i = 0
        var qual = 0.0
        if (len > 0) {
          val row = qualities(qualityOffset)
          while (i < row.length) {
            qual += qualityToLogCorrectProb(row.charAt(i))
            i += 1
          }
        }
        quality.add(qual)
        helper(offset + 1, qualityOffset + (if (len > 0) 1 else 0))
      }
    }
    helper(0, 0)
    quality.toArray
  }

}