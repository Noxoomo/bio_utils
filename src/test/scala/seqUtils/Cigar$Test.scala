package seqUtils

import org.scalatest._
import seqUtils.Cigar._

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 22:23
 */
class Cigar$Test extends FlatSpec with Matchers {
  "Cigar" should "parse cigar-string correctly" in {
    val cigar = "35I69M1I15M1I12M1I29M1I3M3D1M"
    val trueCigar = Array[Int](Cigar.createOperation('I',35), Cigar.createOperation('M',69), Cigar.createOperation('I',1),
    createOperation('M',15), createOperation('I',1), createOperation('M',12), createOperation('I',1), createOperation('M',29),
    createOperation('I',1), createOperation('M',3), createOperation('D',3), createOperation('M',1))

    val parsedCigar = Cigar.parseCigar(cigar)
    assertResult(trueCigar)(parsedCigar.operations)
    assertResult(cigar)(parsedCigar.toString)
  }
  "Cigar" should "create pairwise aligment correctly" in {
    val cigar = "35I69M1I15M1I12M1I29M1I3M3D1M"
    //really reference is second, but in our interfaces Cigar should be operations, which should be applied to read, not to reference
    val referencePaired = "CCACCAAATAAAAAACGCCTTAGTAAGTATTTTTCAGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCCTTCTGAACTGGTTACCCTGCCGTGAGTAAAATTAAAATTTTATTGACTTAGGTCACTAAAA---T"
    val readPaired = "-----------------------------------AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAG-CTTCTGAACTGGTTA-CCTGCCGTGAGT-AAATTAAAATTTTATTGACTTAGGTCACT-AAAGGGT"
    val reference = referencePaired.replace("-", "")
    val read = readPaired.replace("-", "")
    val parsedCigar = Cigar.parseCigar(cigar)

    assertResult((referencePaired, readPaired))(applyCigar(parsedCigar, reference, read))
  }

  "Penalty" should "correctly align strings" in {
    val first = "GCATGCU"
    val second = "GATTACA"
    val aligner = new SeqAligner[Char](new Penalty[Char] {
      final override def apply(v1: Char, v2: Char): Int = if (v1 == v2) -1 else 1
    })
    val aligment: Cigar = aligner(first.toCharArray, second.toCharArray)
    //note: there are more, then one correct aligment
    assertResult(("GCA-TGCU", "G-ATTACA"))(applyCigar(aligment, first, second))
  }
}
