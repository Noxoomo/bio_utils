package seqUtils


import fastFunctions.VecTools


object HRunHelper {
  type HRun = Int
  type HRunSeq = Array[HRun]

  @inline
  def qualityToProb(quality: Char): Double = FastMath.pow(10, (quality.toInt - 33.0) * (-0.1))

  @inline
  def probToQuality(quality: Double): Int = (FastMath.log10(quality) * (-10) + 33).toInt

  @inline
  def qualityToIdx(quality: Char): Int = quality.toInt - 33

  class HRunSeqMeta(val qualities: Array[String], val offsets: Array[Int])

  @inline
  final def hrun(nucl: Char, size: Int): HRun = (size << 10) | nucl.toInt


  @inline
  final def inchrun(run: HRun): HRun = run + (1 << 10)

  @inline
  final def base(hrun: HRun): Char = (hrun & 1023).toChar

  @inline
  final def size(hrun: HRun): Int = hrun >> 10

  final class HRunStringIterator(seq: String, left: Int, right: Int) extends SpecializedIterator[HRun] {
    private var cursor = left

    @tailrec
    private def forwardHelper(count: Int): HRun = {
      val cur = seq.charAt(cursor)
      cursor += 1
      if (cursor < right && cur == seq.charAt(cursor)) {
        forwardHelper(count + 1)
      } else {
        hrun(cur, count + 1)
      }
    }

    override def hasNext: Boolean = cursor < right

    override def next(): Int = forwardHelper(0)
  }

  final class HRunStringReverseIterator(seq: String, left: Int, right: Int) extends SpecializedIterator[HRun] {
    var cursor = right - 1

    @tailrec
    private def backwardHelper(count: Int): HRun = {
      val cur = seq.charAt(cursor)
      cursor -= 1
      if (cursor >= left && cur == seq.charAt(cursor)) {
        backwardHelper(count + 1)
      } else {
        hrun(cur, count + 1)
      }
    }

    override def hasNext: Boolean = cursor >= left

    override def next(): HRun = backwardHelper(0)
  }

  @inline
  final def stringToHRunSeq(seq: String, reverse: Boolean = false): HRunSeq = stringToHRunSeq(seq, 0, seq.length, reverse)

  final def stringToHRunSeq(seq: String, left: Int, right: Int, reverse: Boolean): HRunSeq = {
    val buffer: TIntArrayList = new TIntArrayList(right - left)
    val it = if (reverse) new HRunStringReverseIterator(seq, left, right) else new HRunStringIterator(seq, left, right)
    it.foreach(buffer.add)
    buffer.toArray
  }

  @inline
  final def hrunIterator(seq: String, left: Int, right: Int) = new HRunStringIterator(seq, left, right)

  @inline
  final def hrunReverseIterator(seq: String, left: Int, right: Int) = new HRunStringReverseIterator(seq, left, right)

  implicit class HRunSeqHelper(seq: HRunSeq) {
    @inline
    final def meta(skipFirst: Int, quality: String, reverse: Boolean = false) = {
      val qualities = Array.ofDim[String](seq.length)
      val offsets = Array.ofDim[Int](seq.length)

      @tailrec
      def fillMeta(idx: Int, offset: Int): Unit = {
        if (idx < seq.length) {
          val hrunSize = size(seq(idx))
          qualities(idx) = if (reverse) {
            quality.rsubstring(quality.length - offset - hrunSize - skipFirst, quality.length - offset - skipFirst)
          } else {
            quality.substring(skipFirst + offset, skipFirst + offset + hrunSize)
          }
          offsets(idx) = skipFirst + offset
          fillMeta(idx + 1, offset + hrunSize)
        }
      }
      fillMeta(0, 0)
      new HRunSeqMeta(qualities, offsets)
    }
  }

  final def mkString(cigar: Cigar, reference: HRunSeq, read: HRunSeq) = {
    val refBuilder = new StringBuilder()
    val readBuilder = new StringBuilder()
    val (alignedRef, alignedRead) = align(cigar, reference, read)
    @inline
    def addHRuns(ref: HRun, read: HRun): Unit = {
      val min = FastMath.min(size(ref), size(read))
      val max = FastMath.max(size(ref), size(read))
      var refBase = base(ref)
      var readBase = base(read)
      var i = 0
      while (i < min) {
        refBuilder.append(refBase)
        readBuilder.append(readBase)
        i += 1
      }

      refBase = if (size(ref) == min) '-' else refBase
      readBase = if (size(read) == min) '-' else readBase
      while (i < max) {
        refBuilder.append(refBase)
        readBuilder.append(readBase)
        i += 1
      }
    }
    VecTools.zipForeach(alignedRef, alignedRead, addHRuns)
    (refBuilder.mkString, readBuilder.mkString)
  }


  final def align(cigar: Cigar, reference: HRunSeq, read: HRunSeq): (HRunSeq, HRunSeq) = {
    val refBuilder = new TIntArrayList()
    val readBuilder = new TIntArrayList()

    def addHRuns(ref: HRun, read: HRun, offset: Int): Unit = {
      refBuilder.add(ref)
      readBuilder.add(read)
    }
    val proceeder = new CigarProceeder(addHRuns)
    proceeder.proceed(cigar, reference, read)
    (refBuilder.toArray, readBuilder.toArray)
  }

  final def mkString(hrun: HRun): String = {
    val builder = new StringBuilder
    val hrunBase = base(hrun)
    val hrunSize = size(hrun)
    for (j <- 0 until hrunSize) {
      builder.append(hrunBase)
    }
    builder.toString
  }

  final def mkString(seq: HRunSeq) = {
    val builder = new StringBuilder
    for (hrun <- seq) {
      val hrunBase = base(hrun)
      val hrunSize = size(hrun)
      for (j <- 0 until hrunSize) {
        builder.append(hrunBase)
      }
    }
    builder.mkString
  }

}





