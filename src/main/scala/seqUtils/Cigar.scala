package seqUtils

import fastFunctions.VecTools.VecToolsFastHelper
import gnu.trove.list.array.TIntArrayList

import scala.annotation.tailrec

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 20:12
 */

object Cigar {

  @inline
  final def operationType(i: Int) = (i & 1023).toChar

  @inline
  final def operationSize(i: Int) = i >> 10

  @inline
  final def createOperation(opr: Char, size: Int) = (size << 10) | opr

  private val dict = Set('M', 'N', 'D', 'I', 'S', 'H')
  private val zeroCode = '0'.toInt

  @inline
  final def parseCigar(cigar: String, complementary: Boolean = false) = {
    val entries = new TIntArrayList()
    cigar.view.foldLeft(0) { case (count, c) => {
      if (dict.contains(c)) {
        entries.add(createOperation(c, count))
        0
      } else {
        count * 10 + c.toInt - zeroCode
      }
    }
    }
    if (complementary)
      entries.reverse()
    new Cigar(entries.toArray)
  }

  @inline
  final def applyCigar(cigar: Cigar, reference: String, read: String): (String, String, String) = applyCigar(cigar.operations, reference, read)

  final def applyCigar(cigar: Array[Int], reference: String, read: String) = {
    val refBuilder = new StringBuilder()
    val readBuilder = new StringBuilder()
    val bars = new StringBuilder()
    @tailrec
    def helper(cigarOffset: Int, refOffset: Int, readOffset: Int): Unit = {
      if (cigarOffset < cigar.length) {
        val opr = cigar(cigarOffset)
        operationType(opr) match {
          case 'M' => {
            var i = 0
            val size = operationSize(opr)
            while (i < size) {
              refBuilder.append(reference.charAt(refOffset + i))
              readBuilder.append(read.charAt(readOffset + i))
              bars.append('|')
              i += 1
            }
            helper(cigarOffset + 1, refOffset + i, readOffset + i)
          }
          case ('I' | 'S') => {
            val c = if (operationType(opr) == 'I') '-' else 'N'
            var i = 0
            val size = operationSize(opr)
            while (i < size) {
              refBuilder.append(reference.charAt(refOffset + i))
              readBuilder.append(c)
              bars.append(' ')
              i += 1
            }
            helper(cigarOffset + 1, refOffset + i, readOffset)
          }
          case 'D' => {
            var i = 0
            val size = operationSize(opr)
            while (i < size) {
              readBuilder.append(read.charAt(readOffset + i))
              refBuilder.append('-')
              bars.append(' ')
              i += 1
            }
            helper(cigarOffset + 1, refOffset, readOffset + i)
          }
        }
      }
    }
    helper(0, 0, 0)
    (refBuilder.toString(), readBuilder.toString(), bars.toString())
  }

  import HRunHelper._

  @inline
  def equalNucl(left: Char, right: Char): Boolean = {
    left == '-' || right == '-' || left == right
  }

  class THRunSeqBuilder {
    val buffer = new TIntArrayList()
    var hRun = 0


    @inline
    def endHRun(): Int = {
      val sz = size(hRun)
      if (sz > 0) {
        buffer.add(hRun)
      }
      hRun = 0
      sz
    }

    @inline
    def reset(): Unit = {
      hRun = 0
    }

    @inline
    def addNucl(nucl: HRun, endNucl: Boolean = false): Int = {
      if (base(hRun) != base(nucl) || endNucl) {
        val size = endHRun()
        hRun = nucl
        size
      } else {
        val sz = size(nucl) + size(hRun)
        hRun = hrun(base(nucl), sz)
        0
      }
    }

    @inline
    def lastBase(): Char = if (buffer.size() > 0) base(buffer.get(buffer.size() - 1)) else '-'

    @inline
    def lastSize(): Int = if (buffer.size() > 0) size(buffer.get(buffer.size() - 1)) else 0

    @inline
    def build(): HRunSeq = buffer.toArray
  }

  final def applyAndConvertToHRuns(cgr: Cigar, reference: String, read: String): (HRunSeq, HRunSeq, Cigar) = {
    val cigar = cgr.operations
    val refBuilder = new THRunSeqBuilder()
    val readBuilder = new THRunSeqBuilder()
    val cigarBuilder = new CigarBuilder()

    @inline
    def addCigar(refSize: Int, readSize: Int): Unit = {
      if (refSize != 0 || readSize != 0) {
        if (refSize == 0) {
          cigarBuilder.add('D')
        } else if (readSize == 0) {
          cigarBuilder.add('I')
        } else {
          cigarBuilder.add('M')
        }
      }
    }
    @inline
    def add(ref: HRun, read: HRun): Unit = {
      val equal = base(ref) == base(read)
      val refSize = refBuilder.addNucl(ref, !equal)
      val readSize = readBuilder.addNucl(read, !equal)
      addCigar(refSize, readSize)
    }

    @inline
    def nextBase(nucl: String, offset: Int) = if (offset < nucl.length) nucl(offset) else '-'

    @tailrec
    def helper(cigarOffset: Int, refOffset: Int, readOffset: Int): Unit = {
      if (cigarOffset < cigar.length) {
        val opr = cigar(cigarOffset)
        var i = 0
        val size = operationSize(opr)
        operationType(opr) match {
          case 'M' => {
            while (i < size) {
              add(hrun(reference.charAt(refOffset + i), 1), hrun(read.charAt(readOffset + i), 1))
              i += 1
            }
            helper(cigarOffset + 1, refOffset + i, readOffset + i)
          }
          case 'I' => {
            while (i < size) {
              add(hrun(reference.charAt(refOffset + i), 1), hrun(readBuilder.lastBase(), 0))
              i += 1
            }
            helper(cigarOffset + 1, refOffset + size, readOffset)
          }
          case 'D' => {
            while (i < size) {
              add(hrun(refBuilder.lastBase(), 0), hrun(read.charAt(readOffset + i), 1))
              i += 1
            }
            helper(cigarOffset + 1, refOffset, readOffset + size)
          }
        }
      }
    }


    helper(0, 0, 0)

    val refSize = refBuilder.endHRun()
    val readSize = readBuilder.endHRun()
    addCigar(refSize, readSize)

    val newCigar = cigarBuilder.build()
    (refBuilder.build(), readBuilder.build(), newCigar)
  }
}

class Cigar(val operations: Array[Int]) extends Iterable[Int] {

  final override def iterator: Iterator[Int] = new SpecializedIterator[Int] {
    var offset: Int = 0
    var oprOffset: Int = 0

    final override def next(): Int = {
      val current = operations(offset)
      val opr = Cigar.operationType(current)
      val size = Cigar.operationSize(current)
      oprOffset += 1
      if (oprOffset == size) {
        oprOffset = 0
        offset += 1
      }
      opr
    }

    final override def hasNext: Boolean = offset < operations.length &&
      (offset == operations.length && oprOffset < Cigar.operationSize(operations(offset)))
  }

  lazy val unclippedLengths: Int = {
    var refSize = 0
    var readSize = 0
    @tailrec
    def helper(cigarOffset: Int): Unit = {
      if (cigarOffset < operations.length) {
        val opr = Cigar.operationType(operations(cigarOffset))
        val size = Cigar.operationSize(operations(cigarOffset))
        if (opr == 'M') {
          refSize += size
          readSize += size
          helper(cigarOffset + 1)
        } else if (opr == 'I') {
          readSize += size
          helper(cigarOffset + 1)
        } else if (opr == 'D') {
          helper(cigarOffset + 1)
        }
      }
    }
    helper(0)
    (refSize << 16) | readSize
  }


  @tailrec
  private def counter(cigarOffset: Int, size: Int, c: Char): Int = {
    if (cigarOffset < operations.length) {
      if (Cigar.operationType(operations(cigarOffset)) == c) {
        counter(cigarOffset + 1, size + Cigar.operationSize(operations(cigarOffset)), c)
      } else {
        counter(cigarOffset + 1, size, c)
      }
    } else {
      size
    }
  }

  lazy val matchCount: Int = counter(0, 0, 'M')

  lazy val insertionCount: Int = counter(0, 0, 'I')

  lazy val deletionCount: Int = counter(0, 0, 'D')

  lazy val clippingSize = {
    val opr = operations(operations.length - 1)
    if (Cigar.operationType(opr) == 'S') {
      Cigar.operationSize(opr)
    } else {
      0
    }
  }

  val insertionSize = {
    val opr = operations(0)
    if (Cigar.operationType(opr) == 'I') {
      Cigar.operationSize(opr)
    } else {
      0
    }
  }

  final override def toString() = {
    val builder = new StringBuilder()
    operations.foreach(opr => {
      builder.append(Cigar.operationSize(opr))
      builder.append(Cigar.operationType(opr))
    })
    builder.toString()
  }

  final def inverse(): Cigar = {
    val inversed = operations.fastMap(opr => {
      val oprType = Cigar.operationType(opr)
      val size = Cigar.operationSize(opr)
      if (oprType == 'I') {
        Cigar.createOperation('D', size)
      } else if (oprType == 'D') {
        Cigar.createOperation('I', size)
      } else {
        Cigar.createOperation(oprType, size)
      }
    })
    new Cigar(inversed)
  }
}

class CigarBuilder {
  private val cigar = new TIntArrayList(100)
  private var current: Char = '?'
  private var currentCount: Int = 0

  @inline
  private def endOpr(): Unit = {
    if (currentCount > 0) {
      cigar.add(Cigar.createOperation(current, currentCount))
    }
    currentCount = 1
  }

  @inline
  final def add(opr: Char): Unit = {
    if (opr == current) {
      currentCount += 1
    } else {
      endOpr()
      current = opr
    }
  }

  def size() = cigar.size()

  final def build(): Cigar = {
    endOpr()
    new Cigar(cigar.toArray)
  }
}



