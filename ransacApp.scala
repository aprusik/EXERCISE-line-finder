import scala.collection.mutable.ArrayBuffer
import scala.io.StdIn
import scala.util.Random

object ransacApp {
  val rand = new Random()
  var height: Double = 0
  var width: Double = 0

  def explainedBy(
                   p: (Double, Double),
                   sl: Double,
                   in: Double,
                   maxX: Double = height,
                   maxY: Double = width,
                   minX: Double = 0,
                   minY: Double = 0
                 ): Boolean = {
    val test = math.abs(sl * p._1 - p._2 + in) / math.sqrt(sl * sl + 1)
    val d = 1
    val foo = d >= test &&
      ((p._1 >= minX + d && p._1 <= maxX + d) &&
        (p._2 >= minY + d && p._2 <= maxY + d))
    foo
  }

  def demReg(inliers: Array[(Double, Double)]): (Double, Double) = {
    val n = inliers.length
    val x = inliers.unzip._1
    val y = inliers.unzip._2

    val xBar = (1 / n) * x.sum
    val yBar = (1 / n) * y.sum
    val sXX = (1 / n-1) * x.reduce((x1, x2) => x1 + (x2 - xBar)*(x2 - xBar))
    val sYY = (1 / n-1) * y.reduce((y1, y2) => y1 + (y2 - xBar)*(y2 - xBar))

    var sum: Double = 0
    for (i <- inliers.indices) {
      sum += (inliers(i)._1 - xBar)*(inliers(i)._2 - yBar)
    }
    val sXY = (1 / n-1) * sum

    val slope =
      (sYY - sXX + math.sqrt((sYY - sXX)*(sYY - sXX) + 4*(sXY*sXY))) / (2*sXY)
    val inter = yBar - slope*xBar
    (slope, inter)
  }

  def ransac(startPoints: Array[(Double, Double)]):
  Array[((Double, Double), (Double, Double))] = {
    var points: ArrayBuffer[(Double, Double)] =
      new ArrayBuffer(startPoints.length)
    points ++= startPoints

    var finalSegments: Array[((Double, Double), (Double, Double))] = Array.empty
    def getPoints: ((Double, Double), (Double, Double)) = {
      val p1 = points(rand.nextInt(points.length))
      val p2 = points(rand.nextInt(points.length))

      if (p2._1 == p1._2) getPoints
      else (p1, p2)
    }

    // repeat
    for (count <- 1 to 1000) if (points.nonEmpty) {
      var lineSegs: Array[((Double, Double), (Double, Double), Int)] = Array.empty
      var bm: (Double, Double) = (0, 0)
      //repeat many times
      for (_ <- 1 to 500) {
        //randomly pick minimum data to fit model
        val rp = getPoints
        var sl = (rp._2._2 - rp._1._2) / (rp._2._1 - rp._1._1)
        var in = rp._1._2 - (sl * rp._1._1)

        // find inliers
        var inliers: Array[(Double, Double)] =
          points.filter(explainedBy(_, sl, in)).toArray

        // repeat until no change
        if (inliers.length > 0) {
          var oldSl = sl
          var oldIn = in
          while (oldSl - sl > 0 && oldIn - in > 0) {
            // fit model to inliers
            val dr = demReg(inliers)
            sl = dr._1
            in = dr._2
            oldSl = sl
            oldIn = in
            // find new inliers
            inliers = points.filter(explainedBy(_, sl, in)).toArray
          }
        }

        bm = (sl, in)

//        var lp: (Double, Double) = null
//        var seg: Array[(Double, Double)] = Array.empty
//        var segs: Array[Array[(Double, Double)]] = Array.empty
//        inliers.foreach(p => {
//          if (lp != null && dist(lp, p) > 4) {
//            segs = segs :+ seg
//            seg = Array.empty
//          }
//          seg = seg :+ p
//          lp = p
//        })

//        // if best model has enough inliers
//        segs = segs.filter(seg => seg.length > 7 )

//        if (segs.nonEmpty) {
//          val bestFit = segs.maxBy(_.length)
//          val p1 = bestFit(0)
//          val p2 = bestFit(bestFit.length - 1)
//          lineSegs = (p1, p2, lineSegs.length) +: lineSegs
//        }
      }

      var inliers = points.filter(explainedBy(_, bm._1, bm._2))

      inliers = inliers.map(p => {
        val lpx = p._1 + (bm._1 / (bm._1 * bm._1 + 1)) * (p._2 - bm._2 - bm._1 * p._1)
        val lpy = bm._1 * lpx + bm._2
        (lpx, lpy)
      })

      inliers = inliers.sortWith((p1, p2) => p1._1 < p2._1)

      var lp: (Double, Double) = null
      var seg: Array[(Double, Double)] = Array.empty
      var segs: Array[Array[(Double, Double)]] = Array.empty
      inliers.foreach(p => {
        if (lp != null && dist(lp, p) > 4) {
          segs = segs :+ seg
          seg = Array.empty
        }
        seg = seg :+ p
        lp = p
      })
      if (seg.length > 0) segs = segs :+ seg

      segs = segs.filter(s => s.length >= 8 && {
        s.length / dist(s(0), s(s.length-1)) >= 0.8
      })

      if (segs.nonEmpty) {
        val bestSeg = segs.maxBy(_.length)
        val p1 = bestSeg(0)
        val p2 = bestSeg(bestSeg.length-1)
//        val bestSlope = (p2._2 - p1._2) / (p2._1 - p1._1)
//        val bestInter = p2._2 - (bestSlope * p1._1)
        val f = 5
        val explained = points.filter(p => {
          explainedBy(
            p,
            bm._1,
            bm._2,
            Array(p1._1, p2._1).max,
            Array(p1._2, p2._2).max,
            Array(p1._1, p2._1).min,
            Array(p1._2, p2._2).min
          )
        })
        points --= explained
        finalSegments = (bestSeg(0), bestSeg(bestSeg.length-1)) +: finalSegments
      }
    }

    finalSegments
  }

  def dist(p1: (Double, Double), p2: (Double, Double)): Double = {
    val x = p1._1 - p2._1
    val y = p1._2 - p2._2

    math.sqrt(x*x + y*y)
  }

  def nextLine(): Array[String] = {
    val line = StdIn.readLine()
    if (line == null) Array.empty
    else if (line.charAt(0) == '/') nextLine()
    else line.split(" ")
  }

  def main(args: Array[String]): Unit = {
    height = nextLine()(0).toDouble
    width = nextLine()(0).toDouble

    var points: Array[(Double, Double)] = Array.empty
    var line = nextLine()
    while (line.nonEmpty) {
      points = points :+ (line(0).toDouble, line(1).toDouble)
      line = nextLine()
    }

    val lines = ransac(points)

    println("number of circles: 0")
    println("number of lines: " + lines.length)
    lines.foreach(ls => {
      val x1 = ls._1._1
      val y1 = ls._1._2

      val x2 = ls._2._1
      val y2 = ls._2._2

      println(s"$x1 $y1 $x2 $y2")
    })
  }
}
