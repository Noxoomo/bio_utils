package seqUtils

import seqUtils.HRunHelper._

import scala.annotation.tailrec

/**
 * User: Noxoomo
 * Date: 12.11.15
 * Time: 19:51
 */


class CigarProceeder(hrunProceeder: (HRun, HRun, Int) => Unit) {

  @inline
  final def proceed(cigar: Cigar, reference: HRunSeq, read: HRunSeq) = {
    @tailrec
    def helper(cigarOffset: Int, refOffset: Int, readOffset: Int): Unit = {
      if (cigarOffset < cigar.operations.length) {
        val opr = cigar.operations(cigarOffset)
        val oprType = Cigar.operationType(opr)
        val oprSize = Cigar.operationSize(opr)
        if (oprType == 'M') {
          var i = 0
          while (i < oprSize) {
            hrunProceeder(reference(refOffset + i), read(readOffset + i), readOffset + i)
            i += 1
          }
          helper(cigarOffset + 1, refOffset + i, readOffset + i)
        } else if (oprType == 'I') {
          var i = 0
          while (i < oprSize) {
            val fakeRead = hrun(base(reference(refOffset + i)), 0)
            hrunProceeder(reference(refOffset + i), fakeRead, readOffset + i)
            i += 1
          }
          helper(cigarOffset + 1, refOffset + oprSize, readOffset)
        } else {
          var i = 0
          while (i < oprSize) {
            val fakeRef = hrun(base(read(readOffset + i)), 0)
            hrunProceeder(fakeRef, read(readOffset + i), readOffset + i)
            i += 1
          }
          helper(cigarOffset + 1, refOffset, readOffset + oprSize)
        }
      }
    }
    helper(0, 0, 0)
  }

}
