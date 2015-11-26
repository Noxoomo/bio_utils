package seqUtils

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 17:20
 */

abstract class SpecializedIterator[@specialized T] extends Iterator[T] {

  @inline
  final override def foreach[@specialized U](f: (T) => U): Unit = {
    while (hasNext) f(next())
  }
}
