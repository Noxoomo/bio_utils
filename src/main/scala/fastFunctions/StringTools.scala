package fastFunctions

import scala.reflect.ClassTag

/**
 * User: Noxoomo
 * Date: 25.10.15
 * Time: 0:42
 */
object StringTools {

  implicit class StringToolsFastHelper(str: String) {

    @inline
    final def fastMap[@specialized T <: AnyVal](func: Char => T)(implicit m: ClassTag[T]): Array[T] = {
      val result = new Array[T](str.length)
      var i = 0
      while (i < result.length) {
        result(i) = func(str(i))
        i += 1
      }
      result
    }

    @inline
    final def rsubstring(start: Int, end: Int): String = {
      val result = new StringBuilder()
      var i = end - 1
      while (i >= start) {
        result.append(str(i))
        i -= 1
      }
      result.toString()
    }
  }

  implicit class CharSeqHelper(seq: CharSequence) {

    final def viewSlice(start: Int, end: Int): CharSequence = new CharSequenceSlice(seq, start, end - start)

  }

}

class CharSequenceSlice(val owner: CharSequence, val offset: Int, val length: Int) extends CharSequence {

  override def charAt(index: Int): Char = owner.charAt(index + offset)

  override def subSequence(start: Int, end: Int): CharSequence = new CharSequenceSlice(owner, offset + start, end - start)
}


//char seq should be immutable
trait CharSeqMovingSlice {
  def moveNext(): Unit

  def charAt(index: Int): Char

  def subSequence(start: Int, end: Int): CharSequence
}

class StringMovingSlice(val owner: CharSequence, val startOffset: Int, val length: Int) extends CharSeqMovingSlice {
  private var offset = 0

  override def charAt(index: Int): Char = owner.charAt(startOffset + offset + index)

  override def subSequence(start: Int, end: Int): CharSequence = new CharSequenceSlice(owner, startOffset + offset + start, end - start)

  override def moveNext(): Unit = {
    offset += 1
  }

  override def toString() = owner.subSequence(startOffset + offset, startOffset + offset + length).toString
}
