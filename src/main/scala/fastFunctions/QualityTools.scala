package fastFunctions

import org.apache.commons.math3.util.FastMath

/**
 * User: Noxoomo
 * Date: 02.11.15
 * Time: 1:16
 */
object QualityTools {
  @inline
  def errorProb(qualities: Array[Double]): Double = {
    var noErrorProb = 0.0
    for (quality <- qualities) {
      noErrorProb += FastMath.log(1 - quality)
    }
    1 - FastMath.exp(noErrorProb)
  }

}
