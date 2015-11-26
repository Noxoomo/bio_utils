package fastFunctions

import scala.reflect.ClassTag

/**
 * User: Noxoomo
 * Date: 25.10.15
 * Time: 0:42
 */
object VecTools {

  @inline
  final def sum(arr: Array[Int]): Int = {
    var res = arr(0)
    var i: Int = 1
    while (i < arr.length) {
      res += arr(i)
      i += 1
    }
    res
  }

  @inline
  final def sum(arr: Array[Double]): Double = {
    var res = arr(0)
    var i: Int = 1
    while (i < arr.length) {
      res += arr(i)
      i += 1
    }
    res
  }


  @inline
  final def zipForeach[@specialized T <: AnyVal](left: Array[T], right: Array[T], func: (T, T) => Unit): Unit = {
    assert(left.length == right.length)
    var i = 0
    while (i < left.length) {
      func(left(i), right(i))
      i += 1
    }
  }

  @inline
  final def zipAggregate[@specialized T <: AnyVal](left: Array[T], right: Array[T], func: (T, T, T) => T)(implicit num: Numeric[T]): T = {
    assert(left.length == right.length)
    var i = 0
    var result: T = num.zero
    while (i < left.length) {
      result = func(result, left(i), right(i))
      i += 1
    }
    result
  }

  implicit class VecToolsFastHelper[@specialized T <: AnyVal](arr: Array[T]) {
    @inline
    final def fastMap[@specialized U <: AnyVal](func: T => U)(implicit m: ClassTag[U]): Array[U] = {
      val result = new Array[U](arr.length)
      var i = 0
      while (i < result.length) {
        result(i) = func(arr(i))
        i += 1
      }
      result
    }

    @inline
    final def fastMax()(implicit ord: Ordering[T]): T = {
      var max = arr(0)
      var i = 1
      while (i < arr.length) {
        max = if (ord.gteq(arr(i), max)) arr(i) else max
        i += 1
      }
      max
    }

    @inline
    final def fastMin()(implicit ord: Ordering[T]): T = {
      var min = arr(0)
      var i = 1
      while (i < arr.length) {
        min = if (ord.lteq(arr(i), min)) arr(i) else min
        i += 1
      }
      min
    }

    @inline
    final def fastSum()(implicit num: Numeric[T]): T = {
      var sum = arr(0)
      var i = 1
      while (i < arr.length) {
        sum = num.plus(sum, arr(i))
        i += 1
      }
      sum
    }

    @inline
    final def mkStringFast(sep: Char, endWith: Char = 0): String = {
      val result = new StringBuilder
      var i = 0
      while (i < arr.length - 1) {
        result.append(arr(i)).append(sep)
        i += 1
      }
      if (i < arr.length) {
        result.append(arr(i))
      }
      if (endWith != 0) {
        result.append(endWith)
      }
      result.mkString
    }
  }

}
