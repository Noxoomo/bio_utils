package seqUtils

import gnu.trove.list.array.TCharArrayList
import org.apache.commons.math3.util.FastMath
import seqUtils.HRunHelper._

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 19:54
 */

abstract class Penalty[T] extends ((T, T) => Int) {
  val gapPenalty: Int = 1
}

class NaivePenalty[@specialized T] extends Penalty[T] {
  final override def apply(v1: T, v2: T): Int = if (v1 == v2) -1 else 1
}

class HRunPenalty extends Penalty[Int] {
  override def apply(left: Int, right: Int): Int = {
    val leftBase = base(left)
    val rightBase = base(right)
    val diff: Int = FastMath.abs(size(left) - size(right))
    if (leftBase == rightBase || diff > gapPenalty) {
      diff
    } else {
      gapPenalty
    }
  }

  final override val gapPenalty: Int = 5
}

//our aligner minimizes score/
class SeqAligner[@specialized T](val score: Penalty[T]) extends ((Array[T], Array[T]) => Cigar) {

  @inline
  private def createScoreMatrix(rows: Int, cols: Int): Matrix = {
    val scores = new Matrix(rows, cols)
    var i = 0
    while (i < cols) {
      scores(0, i) = score.gapPenalty * i
      i += 1
    }

    i = 0
    while (i < rows) {
      scores(i, 0) = score.gapPenalty * i
      i += 1
    }
    scores
  }

  @inline
  private def min(a: Int, b: Int, c: Int): Int = {
    if (a < b) {
      if (a < c) a else c
    } else {
      if (b < c) b else c
    }
  }


  final override def apply(reference: Array[T], read: Array[T]): Cigar = {
    val scores = createScoreMatrix(reference.length + 1, read.length + 1)
    @inline
    def matchScore(i: Int, j: Int): Int = scores(i - 1, j - 1) + score(reference(i - 1), read(j - 1))
    @inline
    def delScore(i: Int, j: Int): Int = scores(i - 1, j) + score.gapPenalty //delete from reference
    @inline
    def insScore(i: Int, j: Int): Int = scores(i, j - 1) + score.gapPenalty //insert to reference

    {
      var i = 1
      while (i <= reference.length) {
        var j = 1
        while (j <= read.length) {
          scores(i, j) = min(matchScore(i, j), delScore(i, j), insScore(i, j))
          j += 1
        }
        i += 1
      }
    }

    val operations: TCharArrayList = new TCharArrayList(1 + FastMath.max(reference.length, read.length))

    {
      var i = reference.length
      var j = read.length
      //note: Cigar — how to edit read string
      while (i > 0 | j > 0) {
        if (i > 0 && j > 0 && (scores(i, j) == matchScore(i, j))) {
          operations.add('M')
          i -= 1
          j -= 1
        }
        else if (i > 0 && (scores(i, j) == delScore(i, j))) {
          operations.add('I')
          i -= 1
        } else {
          operations.add('D')
          j -= 1
        }
      }
    }

    {
      val aligmentBuilder = new CigarBuilder()
      var i = operations.size() - 1
      while (i >= 0) {
        aligmentBuilder.add(operations.get(i))
        i -= 1
      }
      aligmentBuilder.build()
    }
  }

}
