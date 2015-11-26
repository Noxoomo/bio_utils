package seqUtils

import com.googlecode.concurrenttrees.radix.ConcurrentRadixTree
import com.googlecode.concurrenttrees.radix.node.concrete.SmartArrayBasedNodeFactory
import fastFunctions.StringTools._
import seqUtils.HRunHelper._

/**
 * User: Noxoomo
 * Date: 26.11.15
 * Time: 21:39
 */

class RadixTreeHRunSet(val genom: String, val maxKMerSize: Int = 128) extends Set[HRunSeq] {

  case class Value()

  val emptyValue = Value()
  val trie = {
    val tree = new ConcurrentRadixTree[Value](new SmartArrayBasedNodeFactory)
    val movingSlice = genom.movingSlice(0, maxKMerSize)
    for (i <- 0 until genom.length - maxKMerSize) {
      tree.putIfAbsent(movingSlice, emptyValue)
      movingSlice.moveNext()
    }
    tree
  }

  override def contains(hkmer: HRunSeq): Boolean = trie.getKeysStartingWith(HRunHelper.mkString(hkmer)).iterator().hasNext

  override def +(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def -(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def iterator: Iterator[HRunSeq] = throw new UnsupportedOperationException
}