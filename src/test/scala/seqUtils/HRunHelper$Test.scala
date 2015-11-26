package seqUtils

import org.scalatest.{FlatSpec, Matchers}
import seqUtils.HRunHelper._

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 17:36
 */
class HRunHelper$Test extends FlatSpec with Matchers {


  "HRunHelper" should "should correctly convert to HRunSeq" in {
    val seq = "AAAGCCCCTTATATTGGGGGG"
    val qual = "987210100145678987654"


    val res = Array[Int](hrun('A', 3), hrun('G',1), hrun('C',4), hrun('T',2), hrun('A',1), hrun('T',1),hrun('A',1),hrun('T',2),hrun('G',6))
    assertResult(stringToHRunSeq(seq, 0, seq.length, false))(res)
    assertResult(stringToHRunSeq(seq.reverse, 0, seq.length, reverse = true))(res)

    assertResult(stringToHRunSeq(seq, 3, seq.length, false))(res.slice(1, res.length))
    assertResult(stringToHRunSeq(seq, 19, seq.length, false))(Array[Int](hrun('G',2)))
    assertResult(stringToHRunSeq(seq, 19, seq.length, reverse = true))(Array[Int](hrun('G',2)))
    assertResult(stringToHRunSeq(seq, 0, 2, reverse = true))(Array[Int](hrun('A',2)))
    assertResult(stringToHRunSeq(seq, 0, 4, reverse = true))(Array[Int](hrun('G',1), hrun('A',3)))


    val meta0 = stringToHRunSeq(seq, 0, seq.length, false).meta(0, qual)
    val qual0 = Array("987", "2", "1010", "01", "4", "5", "6", "78", "987654")
    assertResult(qual0)(meta0.qualities)
    val meta1 = stringToHRunSeq(seq, 2, seq.length, false).meta(2, qual)
    val qual1 = Array("7", "2", "1010", "01", "4", "5", "6", "78", "987654")
    assertResult(qual1)(meta1.qualities)


    val meta2 = stringToHRunSeq(seq, 0, seq.length, true).meta(0, qual, true)
    val qual2 = Array("456789", "87", "6", "5", "4", "10", "0101", "2", "789")
    assertResult(qual2)(meta2.qualities)
    val meta3 = stringToHRunSeq(seq, 2, seq.length, true).meta(0, qual, true)
    val qual3 = Array("456789", "87", "6", "5", "4", "10", "0101", "2", "7")
    assertResult(qual3)(meta3.qualities)


  }


}
