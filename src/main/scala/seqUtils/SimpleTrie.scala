package seqUtils

import com.googlecode.concurrenttrees.radix.ConcurrentRadixTree
import com.googlecode.concurrenttrees.radix.node.concrete.SmartArrayBasedNodeFactory
import fastFunctions.StringTools.CharSeqHelper

/**
 * User: Noxoomo
 * Date: 26.11.15
 * Time: 21:39
 */


class SimpleTrie(val genom: String, val maxKMerSize: Int = 32) extends Set[CharSequence] {

  case class Value()

  val emptyValue = Value()

  val trie = {
    val tree = new ConcurrentRadixTree[Value](new SmartArrayBasedNodeFactory)

    def add(str: CharSequence) {
      for (i <- 0 until str.length) {
        val end = if (i + maxKMerSize < str.length()) i + maxKMerSize else str.length
        tree.putIfAbsent(str.viewSlice(i, end), emptyValue)
      }
    }
    add(genom)
    tree
  }

  def contains(seq: CharSequence): Boolean = {
    if (seq.length == 0) {
      true
    } else if (seq.length() <= maxKMerSize) trie.getKeysStartingWith(seq).iterator().hasNext
    else {
      print(f"Long hkmer: ${seq.length}\n")
      genom.contains(seq)
    }
  }

  override def +(elem: CharSequence): Set[CharSequence] = throw new UnsupportedOperationException

  override def -(elem: CharSequence): Set[CharSequence] = throw new UnsupportedOperationException

  override def iterator: Iterator[CharSequence] = throw new UnsupportedOperationException
}