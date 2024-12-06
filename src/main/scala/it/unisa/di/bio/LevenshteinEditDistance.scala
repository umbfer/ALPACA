package it.unisa.di.bio


import org.apache.hadoop.fs.{FileSystem, Path}
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.commons.io.FilenameUtils.{getBaseName, getFullPath}

import java.io.{BufferedReader, BufferedWriter, File, FileNotFoundException, FileWriter, IOException, OutputStreamWriter, StringReader}
import scala.io.BufferedSource
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.{LongWritable, Text}
import org.apache.hadoop.mapreduce.lib.input.NLineInputFormat

import java.net.URI
import scala.collection.mutable.ArrayBuffer



object LevenshteinEditDistance {

  var local = true

  var sc:  SparkContext = null

  var hadoopConf: org.apache.hadoop.conf.Configuration = null

  class SparseMatrix(var size: Int) {

    private val map = new Long2IntOpenHashMap( size)

    def get(x: Int, y: Int ): Int = {
      val key : Long = (x.toLong << 32) + y
      val ret = map.getOrDefault(key, -1)
      return ret
    }
    def set(x: Int, y: Int, value: Int ) = {
      val key : Long = (x.toLong << 32) + y
      map.put( key, value)
    }
  }


  def main(args: Array[String]) {

    if (args.length < 3) {
      System.err.println("Errore nei parametri sulla command line")
      System.err.println("Usage:\nit.unisa.di.bio.LevenshteinEditDistance inputDataSet [[local|yarn]")
      throw new IllegalArgumentException(s"illegal number of argouments. ${args.length} should be at least 2")
      return
    }

    val inputPath = args(0)
    val outDir = args(1)
    local = (args(2) == "local")

    val sparkConf = new SparkConf().setAppName("LevenshteinEditDistance")
      .setMaster(if (local) "local" else "yarn")
    sc = new SparkContext(sparkConf)
    hadoopConf = sc.hadoopConfiguration

    println(s"***App ${this.getClass.getCanonicalName} Started***")

    val meta = "-1000.10000."
    val AM1Prefix = "MotifRepl-U"
    val AM2Prefix = "PatTransf-U"
    val NM = List(s"${inputPath}/Uniform${meta}fasta")

    val gValues = List("010", "050", "100")

    val AM1 = gValues.map(x => s"${inputPath}/${AM1Prefix}${meta}G=0.${x}.fasta")
    val AM2 = gValues.map(x => s"${inputPath}/${AM2Prefix}${meta}G=0.${x}.fasta")

    val dsList = NM ++ AM1 ++ AM2

    hadoopConf.setInt("mapreduce.input.lineinputformat.linespermap", 4);

    val records = sc.newAPIHadoopFile(inputPath,
        classOf[NLineInputFormat], classOf[LongWritable], classOf[Text], hadoopConf)

    // val records = sc.wholeTextFiles(inputPath)

    println(s"Num Partition = ${records.partitions.size}")

    var ris = records.mapPartitions( processDataset).saveAsTextFile(outDir)

//    ris.foreach( x => {
//      println(s"records: ${x._1} - ${x._2}")
//    })
  }


  def processDataset( it: Iterator[(LongWritable, Text)]) : Iterator[(String, Double)] = {

    val buf = ArrayBuffer[(String, Double)]()
    val params = Array.ofDim[String](4)
    // println(s"*** Slave Started len = ${it.length} ***")

    var p = 0
    while( it.hasNext) {
      params(p) = it.next()._2.toString
      p += 1
    }
    if (p != 4) {
      println("Split failed")
      throw( new Exception("Split Error"))
    }
    val d = if (params(1).length > 10000) distanceSparse(params(1), params(3)) else distanceMatrix( params(1), params(3))
    val distance = d.toDouble / (params(1).length + params(3).length)

    buf += ((params(0) + params(2), distance))

    println(s"*** Slave Finished ***")
    return buf.iterator
  }



  def LevenshtainDistanceSequential( ds: String, local: Boolean) : Int = {

    var seq1: String = null
    var seq2: String = null
    val outputPath = s"${getFullPath(ds)}${getBaseName(ds)}.csv"

    try {
      val writer = if (local)
            new BufferedWriter(new FileWriter(new File(outputPath)))
          else
            new BufferedWriter(
              new OutputStreamWriter(FileSystem.get(URI.create(outputPath), hadoopConf)
                .create(new Path(outputPath))))

      val reader: BufferedSource =
        if (local)
          // legge dal filesystem locale
          scala.io.Source.fromFile(ds)
        else
          // solo questo codice legge dal HDFS Source.fromFile legge solo da file system locali
          // val hdfs = FileSystem.get(new URI("hdfs://master:8020/"), new Configuration())
          scala.io.Source.fromInputStream(FileSystem.get(hadoopConf)
            .open(new Path( ds)))

      var i = 1
      val it : Iterator[String] = reader.getLines()
      // per tutte le sequenze nel file input
      while (it.hasNext) {
        var seqHeader = it.next()
        seq1 = it.next().toString
        seqHeader = it.next();
        seq2 = it.next().toString
        i += 1

        var startTime = System.currentTimeMillis()
        var d = distanceMatrix(seq1, seq2)
        var totTime = System.currentTimeMillis() - startTime
        println (s"The edit distance (matrix) between seq1 and seq2 (len=${seq1.length}) is ${d}, delay:${totTime} msec.")
        val distance = s"${(seq1.length+seq2.length)/d.toDouble}\n"
        writer.write( distance)
        writer.flush()

        // startTime = System.currentTimeMillis()
        // d = distanceSparse(seq1, seq2)
        // totTime = System.currentTimeMillis() - startTime
        // println (s"The edit distance (sparse) between seq1 and seq2 (len=${seq2.length}) is ${d}, delay:${totTime/1000} sec.")
      }
      writer.close()
      println(s"*** Slave Finished ***")
      return i
    }
    catch {
      case x: FileNotFoundException => {
        println(s"Exception: ${x.getMessage}")
        x.printStackTrace()
        println(s"dataset ${ds} or ${outputPath} not found")
        return -1
      }
      case x: IOException   => {
        println(s"Exception: ${x.getMessage}")
        x.printStackTrace()
        println("Input/output Exception")
        return -1
      }
    }
  }



  def  distanceSparse (word1: String,word2: String) : Int = {

     val matrix = new SparseMatrix(word1.length * 10)

    for (i <- 0 to word1.length) {
      // matrix(i)(0) = i
      matrix.set(i, 0, i)
    }
    for (i <- 0 to word2.length) {
      // matrix(0)(i) = i
      matrix.set(0, i, i)
    }

    for (i <- 1 to word1.length) {
      var c1: Char = 0

      c1 = word1.charAt(i - 1)
      for (j <- 1 to word2.length) {
        var c2: Char = 0

        c2 = word2.charAt(j - 1)
        if (c1 == c2) {
          matrix.set(i, j, matrix.get(i - 1, j - 1))
        }
        else {
          val delete = matrix.get(i - 1, j) + 1
          val insert = matrix.get(i, j - 1) + 1
          val substitute = matrix.get(i - 1, j - 1) + 1
          var minimum = delete

          if (insert < minimum) {
            minimum = insert
          }
          if (substitute < minimum) {
            minimum = substitute
          }
          matrix.set(i, j, minimum)
        }
      }
    }
    return matrix.get(word1.length, word2.length)
  }


  def  distanceMatrix (word1: String,word2: String) : Int = {

    val matrix = Array.ofDim[Int](word1.length + 1, word2.length + 1)

    for (i <- 0 to word1.length) {
      matrix(i)(0) = i
    }
    for (i <- 0 to word2.length) {
      matrix(0)(i) = i
    }

    for (i <- 1 to word1.length) {
      var c1: Char = 0

      c1 = word1.charAt(i - 1)
      for (j <- 1 to word2.length) {
        var c2: Char = 0

        c2 = word2.charAt(j - 1)
        if (c1 == c2) {
          matrix(i)(j) = matrix(i-1)(j-1)
        }
        else {
          val delete = matrix(i-1)(j) + 1
          val insert = matrix(i)(j-1) + 1
          val substitute = matrix(i-1)(j-1) + 1
          var minimum = delete

          if (insert < minimum) {
            minimum = insert
          }
          if (substitute < minimum) {
            minimum = substitute
          }
          matrix(i)(j) = minimum
        }
      }
    }
    return matrix(word1.length)(word2.length)
  }
}
