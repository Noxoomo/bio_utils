package seqUtils

/**
 * User: Noxoomo
 * Date: 24.10.15
 * Time: 19:54
 */
class Matrix(rows: Int, columns: Int) {
  private val data = Array.ofDim[Int](columns*rows)

  @inline
  private def idx(i: Int, j: Int) = i * columns + j

  @inline
  final def apply(i: Int, j: Int): Int = data(idx(i, j))

  @inline
  final def update(i: Int, j: Int, value: Int): Unit = {
    data(idx(i, j)) = value
  }
}
