package seqUtils

import org.scalatest.{FlatSpec, Matchers}
import seqUtils.HRunHelper._

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 22:10
 */
class SeqAlignerTest extends FlatSpec with Matchers {

  "SeqAligment" should "should correctly align strings" in {
    val seq = "AAAGCCCCTTATATTGGGGGG"
    val res = Array[Int](hrun('A', 3), hrun('G',1), hrun('C',4), hrun('T',2), hrun('A',1), hrun('T',1),hrun('A',1),hrun('T',2),hrun('G',6))
    assertResult(stringToHRunSeq(seq, 0, seq.length, reverse = false))(res)
    assertResult(stringToHRunSeq(seq.reverse, 0, seq.length, reverse = true))(res)

    assertResult(stringToHRunSeq(seq, 3, seq.length, reverse = false))(res.slice(1, res.length))
    assertResult(stringToHRunSeq(seq, 19, seq.length, false))(Array[Int](hrun('G',2)))
    assertResult(stringToHRunSeq(seq, 19, seq.length, reverse = true))(Array[Int](hrun('G',2)))
    assertResult(stringToHRunSeq(seq, 0, 2, reverse = true))(Array[Int](hrun('A',2)))
    assertResult(stringToHRunSeq(seq, 0, 4, reverse = true))(Array[Int](hrun('G',1), hrun('A',3)))

  }

}
