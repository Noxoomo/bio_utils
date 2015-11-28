package seqUtils

import com.googlecode.concurrenttrees.radix.ConcurrentRadixTree
import com.googlecode.concurrenttrees.radix.node.concrete.SmartArrayBasedNodeFactory
import fastFunctions.StringTools.CharSeqHelper
import seqUtils.HRunHelper.{stringToHRunSeq, _}
/**
 * User: Noxoomo
 * Date: 26.11.15
 * Time: 21:39
 */


class RadixTreeHRunSet(val genom: String, val maxKMerSize: Int = 16) extends Set[HRunSeq] {

  val genomSeq = mapToCharSeq(stringToHRunSeq(genom))
  val genomSeqReverse = mapToCharSeq(stringToHRunSeq(genom, reverse = true))

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
    add(genomSeq)
    add(genomSeqReverse)
    tree
  }

  def mapToCharSeq(hkmerSeq: HRunSeq): CharSequence = {
    new ArrayCharSequence(hkmerSeq.map(hrun => {
      val b = base(hrun)
      val sz = HRunHelper.size(hrun)
      assert(sz < 32)
      val id: Int = if (b == 'A') {
        0
      } else if (b == 'T') {
        1
      } else if (b == 'C') {
        2
      } else if (b == 'G') {
        3
      } else {
        throw new RuntimeException("wrong argument")
      }
      ((sz << 2) | id).toChar
    }))
  }

  override def contains(hkmer: HRunSeq): Boolean = {
    val seq = mapToCharSeq(hkmer)
    if (seq.length <= maxKMerSize) trie.getKeysStartingWith(seq).iterator().hasNext else throw new RuntimeException("wrong hkmer size")
  }

  override def +(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def -(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def iterator: Iterator[HRunSeq] = throw new UnsupportedOperationException
}