package seqUtils

import com.googlecode.concurrenttrees.radix.ConcurrentRadixTree
import com.googlecode.concurrenttrees.radix.node.concrete.SmartArrayBasedNodeFactory
import fastFunctions.StringTools.CharSeqHelper
import gnu.trove.list.array.TCharArrayList
import seqUtils.HRunHelper._
/**
 * User: Noxoomo
 * Date: 26.11.15
 * Time: 21:39
 */


class RadixTreeHRunSet(val genom: String, val maxKMerSize: Int = 32) extends Set[HRunSeq] {

  val genomSeq = genom
  val genomSeqReverse = genom.reverse

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
    val result = new TCharArrayList()
    for (hrun <- hkmerSeq) {
      val sz = HRunHelper.size(hrun)
      val b = base(hrun)
      for (i <- 0 until sz) {
        result.add(b)
      }
    }
    new ArrayCharSequence(result.toArray)
  }

  def contains(seq: CharSequence): Boolean = {
    if (seq.length == 0) {
      true
    } else if (seq.length() <= maxKMerSize) trie.getKeysStartingWith(seq).iterator().hasNext
    else {
      print(f"Long hkmer: ${seq.length}")
      genomSeq.contains(seq) || genomSeqReverse.contains(seq)
    }
  }

  override def contains(hkmer: HRunSeq): Boolean = contains(mapToCharSeq(hkmer))

  override def +(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def -(elem: HRunSeq): Set[HRunSeq] = throw new UnsupportedOperationException

  override def iterator: Iterator[HRunSeq] = throw new UnsupportedOperationException
}