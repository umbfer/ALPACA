//
// $Id: DatasetBuilder.scala 1716 2020-08-07 15:25:19Z cattaneo@dia.unisa.it $
//


package it.unisa.di.bio

import java.io.{BufferedWriter, File, FileNotFoundException, IOException, OutputStreamWriter}
import java.net.URI
import java.time.LocalDateTime
import java.util.Properties
import java.time.format.DateTimeFormatter
import org.apache.spark.{SparkConf, SparkContext, sql}
import it.unisa.di.bio.Misc._
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.fs.{FileSystem, Path}
import org.apache.spark.rdd.RDD

import scala.collection.mutable.{ArrayBuffer, ListBuffer}
import scala.io.BufferedSource
import scala.math._
import scala.util.control.Breaks._
import com.concurrentthought.cla._





object DatasetBuilder {

  var local = true
  var debug = true

  var sc:  SparkContext = null
  val fileExt: String = "fasta"
  var savePath: String = "data/powerstatistics"
  var appProperties: Properties = null

  var numberOfPairs = 0
  var seqLen = 0
  var patternLen = 5
  var geneSize = 3

  var motif = Array('A', 'C', 'C', 'C', 'G')

  var parsed: Args = null

  var gValues: Array[Double] = null
  val uniformDist: Array[Double] = new Array[Double](4)
  val GCReachDist: Array[Double] = new Array[Double](4)
  val mitocondriDist: Array[Double] = new Array[Double](4)
  val shigellaDist: Array[Double] = new Array[Double](4)

  var hadoopConf: org.apache.hadoop.conf.Configuration = null

  var multifile : Boolean = false
  var writer : BufferedWriter = null
  var reader : BufferedSource = null
  var inputIter : Iterator[String] = null
  val rng: scala.util.Random = new scala.util.Random(System.currentTimeMillis / 1000)



  def main(args: Array[String]) {

    val initialArgs: Args =
      """
        |powerstatistics-1.3-SNAPSHOT.jar -o destDir -g method -m yarn|local -f start -t last -s step -n numPairs [-b geneSize -p patternSize]
        |  [-f | --multiFile flag]   specify the format of target dataset -f implies one file for each pair
        |  [-d | --debug flag]    specify debug mode some statistical information are printed (slower)
        |  -o | --output     string  Path to destination Directory
        |  -m | --mode       string  Spark cluster mode local|yarn
        |  -g | --generator  string  Generator selection uniform|eColiShuffled|synthetic|mitocondri|shigella
        |  -f | --from-len   int=10000 from sequence length
        |  -t | --to-len     int     to sequence length
        |  -s | --step       int     step size
        |  -n | --pairs      int     number of pairs
        |  [-b | --geneSize  int=3]  size of imported block from echerichiaColi
        |  [-p | --patternSize int]  size of pattern for both alternative models
        |""".stripMargin.toArgs


    parsed = initialArgs.process(args)

    // If here, successfully parsed the args and none where "--help" or "-h".
    parsed.printAllValues()

    savePath = parsed.getOrElse( "output", "out")
    val datasetType = parsed.getOrElse("generator", "uniform")
    // gamma = args(4).toDouble
    // motif = getMotif(args(4))

    local = parsed.getOrElse("mode", "").compareTo("local") == 0

    appProperties = new Properties
    appProperties.load(this.getClass.getClassLoader.getResourceAsStream("PowerStatistics.properties"))

    val sparkConf = new SparkConf().setAppName(appProperties.getProperty("powerstatistics.datasetBuilder.appName"))
      .setMaster(if (local) "local" else "yarn")
    sc = new SparkContext(sparkConf)
    hadoopConf = sc.hadoopConfiguration

    println(s"***App ${this.getClass.getCanonicalName} Started***")

    numberOfPairs = parsed.getOrElse( "pairs", appProperties.getProperty("powerstatistics.datasetBuilder.numberOfPairs").toInt)

    val params = Array[Int]( parsed.getOrElse("from-len", appProperties.getProperty("powerstatistics.datasetBuilder.startSequenceLength").toInt),
                             parsed.getOrElse("to-len", appProperties.getProperty("powerstatistics.datasetBuilder.maxSequenceLength").toInt),
                             parsed.getOrElse("step", appProperties.getProperty("powerstatistics.datasetBuilder.lengthStep").toInt))

    patternLen = parsed.getOrElse("patternSize", appProperties.getProperty("powerstatistics.datasetBuilder.replacePatternLength").toInt)
    geneSize = parsed.getOrElse( "geneSize", appProperties.getProperty("powerstatistics.datasetBuilder.geneSize").toInt)
    multifile = parsed.getOrElse( "multiFile", false)
    debug = parsed.getOrElse("debug", false)

    // Array(0.25, 0.25, 0.25, 0.25) // P(A), P(C), P(G), P(T) probability
    var stVal = appProperties.getProperty("powerstatistics.datasetBuilder.uniformDistribution").split(",")
    for (i <- 0 to 3) {
      uniformDist(i) = stVal(i).toDouble
    }

    // Array(0.166666666666667, 0.333333333333333, 0.333333333333333, 0.166666666666667) // P(A), P(C), P(G), P(T) probability
    stVal = appProperties.getProperty("powerstatistics.datasetBuilder.GCReachDistribution").split(",")
    for (i <- 0 to 3) {
      GCReachDist(i) = stVal(i).toDouble
    }
    stVal = appProperties.getProperty("powerstatistics.datasetBuilder.mitocondriEmpiricalDistribution").split(",")
    for (i <- 0 to 3) {
      mitocondriDist(i) = stVal(i).toDouble
    }
    stVal = appProperties.getProperty("powerstatistics.datasetBuilder.shigellaEmpiricalDistribution").split(",")
    for (i <- 0 to 3) {
      shigellaDist(i) = stVal(i).toDouble
    }

    // Array(0.001, 0.005, 0.01, 0.05, 0.1)
    stVal = appProperties.getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",")
    gValues = new Array[Double](stVal.length)
    for (i <- 0 to stVal.length - 1) {
      gValues(i) = stVal(i).toDouble
    }

    motif = appProperties.getProperty("powerstatistics.datasetBuilder.motif").toCharArray

    val bs = appProperties.getProperty("powerstatistics.datasetBuilder.defaultBlockSize").toInt
    hadoopConf.setInt("dfs.blocksize ", bs)

    datasetType match {
      case x if (x.compareTo("eColiShuffled") == 0) => datasetEColi(
        appProperties.getProperty("powerstatistics.datasetBuilder.ecoliPrefix"), geneSize, params)

      case x if (x.compareTo("uniform") == 0) => datasetUniform(
        appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"), uniformDist, params)

      case x if (x.compareTo("synthetic") == 0) => dataset2(
        appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"), uniformDist, params)

      case x if (x.compareTo("syntheticMitocondri") == 0) => dataset2(
        appProperties.getProperty("powerstatistics.datasetBuilder.synthMitocondriPrefix"), mitocondriDist, params)

      case x if (x.compareTo("syntheticShigella") == 0) => dataset2(
        appProperties.getProperty("powerstatistics.datasetBuilder.synthShigellaPrefix"), shigellaDist, params)

      case y if (y.compareTo("mitocondri") == 0) => semiSynthetic(
        appProperties.getProperty("powerstatistics.datasetBuilder.syntheticNullModelPrefix.mitocondri"), params)

      case y if (y.compareTo("shigella") == 0) => semiSynthetic(
        appProperties.getProperty("powerstatistics.datasetBuilder.syntheticNullModelPrefix.shigella"), params)

      case z => println(s"${z} must be in: detailed | eColiShuffled | synthetic | syntheticMitocondri | syntheticShigella | Mitocondri | Shigella");
        sc.stop()
    }
  }


  def datasetEColi( nullModelPrefix: String, geneSize: Int, lengths: Array[Int]) : Unit = {

    val fromLen =   lengths(0)
    val maxSeqLen = lengths(1)
    val step =      lengths(2)

    seqLen = fromLen
    while (seqLen <= maxSeqLen ) {
      buildEColiDataset(seqLen, nullModelPrefix, geneSize)
      seqLen = seqLen + step;
    }
  }


  def datasetUniform(nullModelPrefix: String, distribution: Array[Double], lengths: Array[Int]) : Unit = {

    val fromLen =   lengths(0)
    val maxSeqLen = lengths(1)
    val step =      lengths(2)

    seqLen = fromLen
    while (seqLen <= maxSeqLen ) {
      buildSyntenthicDataset(seqLen, nullModelPrefix, distribution)
      seqLen = seqLen + step;
    }
  }


  def dataset2( nullModelPrefix: String, distribution: Array[Double], lengths: Array[Int]) : Unit = {

    val fromLen =   lengths(0)
    val maxSeqLen = lengths(1)
    val step =      lengths(2)

    seqLen = fromLen
    val meanValues = Array(10, 25, 50, 75)
    var bl = 10000
    while (bl <= maxSeqLen * 10) { // l'ultimo sarà il 10% del successore del maxSeqSize

      if (bl >= 1000000)
          hadoopConf.setInt("dfs.blocksize", 67108864)

      for (m <- meanValues) {
        seqLen = bl / 100 * m // prima la divisione altrimenti va in overflow l'integer

        if (seqLen <= maxSeqLen) {
          buildSyntenthicDataset(seqLen, nullModelPrefix, distribution)
        }
      }
      bl = bl * 10;
    }
  }

  //
  // Build the NullModel shuffling a real sequence (EscherichiaColi) usign subsequence of length geneSize
  //
  def buildEColiDataset( targetLen: Int, nullModelPrefix: String, geneSize: Int) : Unit = {

    val st = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").format(LocalDateTime.now)
    println(s"${st} *** Building EColi Derived Dataset (Null Model + Alternate Model) <- ${nullModelPrefix} #pairs: ${numberOfPairs} for len: ${targetLen} ***")

    seqLen = targetLen
    // crea the NullModel
    buildDatasetWithShuffle(numberOfPairs, seqLen, geneSize, nullModelPrefix)

    // crea un altro dataset uniform per il type 1 check
    buildDatasetWithShuffle(numberOfPairs, seqLen, geneSize,
      nullModelPrefix + appProperties.getProperty("powerstatistics.datasetBuilder.Type1CheckSuffix"))

    for (g <- gValues) {

      // build both the Alternative Models using the NullModel just created
      buildDatasetMotifReplace(seqLen, motif, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))

      buildDatasetPatternTransfer(seqLen, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))

    }
  }


  def buildSyntenthicDataset( targetLen: Int, nullModelPrefix: String, distribution: Array[Double]) : Unit = {

    val st = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").format(LocalDateTime.now)
    println(s"${st} *** Building Synthetic Dataset (Null Model + Alternate Model) <- ${nullModelPrefix} #pairs: ${numberOfPairs} for len: ${targetLen} ***")

    seqLen = targetLen
    // crea il dataset uniform
    buildDatasetWithDistribution(numberOfPairs, seqLen, distribution, nullModelPrefix)

    // crea un altro dataset uniform per il type 1 check
    buildDatasetWithDistribution(numberOfPairs, seqLen, distribution,
        nullModelPrefix + appProperties.getProperty("powerstatistics.datasetBuilder.Type1CheckSuffix"))

    for (g <- gValues) {

      // generiamo entrambi i dataset dipendenti dal nullModel Uniform
      buildDatasetMotifReplace(seqLen, motif, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))

      buildDatasetPatternTransfer(seqLen, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))

    }
  }


  def semiSynthetic( nullModelPrefix: String, lengths: Array[Int]) : Unit = {

    val fromLen =   lengths(0)
    val maxSeqLen = lengths(1)
    val step =      lengths(2)

    val meanValues = Array(10, 25, 50, 75)
    var bl = fromLen
    while (bl <= maxSeqLen * 10) { // l'ultimo sarà il 10% del successore del maxSeqSize

      if (bl >= 1000000)
        hadoopConf.setInt("dfs.blocksize", 67108864)

      for (m <- meanValues) {
        seqLen = bl * m / 100

        if (seqLen <= maxSeqLen) {
          buildAlternateModelsOnly(seqLen: Int, nullModelPrefix)
        }
      }
      bl = bl * 10;
    }
  }


  def buildAlternateModelsOnly( targetLen: Int, nullModelPrefix: String) : Unit ={

    val st = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").format(LocalDateTime.now)
    println(s"${st} *** Building Alternate Models <- ${nullModelPrefix} for len: ${targetLen} ***")

    seqLen = targetLen

    for (g <- gValues) {

      // generiamo entrambi i dataset dipendenti dal nullModel (rnd_mitocondri o rnd_shigella)
      buildDatasetMotifReplace(seqLen, motif, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))

      buildDatasetPatternTransfer(seqLen, g, nullModelPrefix,
        appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))

    }

  }


  def reinertDataset(): Unit = {

    for (j <- 1 to 8) {

      seqLen = scala.math.pow(2, j).toInt * 100
      // crea il dataset uniform

      println(s"*** Process starting for len: ${seqLen} ***")

      buildDatasetWithDistribution(numberOfPairs, seqLen, uniformDist,
        appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"))

      // genera GCRich Dataset
      buildDatasetWithDistribution(numberOfPairs, seqLen, GCReachDist,
        appProperties.getProperty("powerstatistics.datasetBuilder.GCReachPrefix"))

      for (g <- gValues) {

        buildDatasetMotifReplace(seqLen, motif, g,
          appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"),
          appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))

        buildDatasetPatternTransfer(seqLen, g,
          appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"),
          appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))

        buildDatasetMotifReplace(seqLen, motif, g,
          appProperties.getProperty("powerstatistics.datasetBuilder.GCReachPrefix"),
          appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))

        buildDatasetPatternTransfer(seqLen, g,
          appProperties.getProperty("powerstatistics.datasetBuilder.GCReachPrefix"),
          appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))
      }
    }
  }


  // This function builds the dataset as a sequence of closures evaluation ... to improve parallelism
  // n.b. useless is I/O Bounded
//  def dataset3() : Unit = {
//
//    val oom = Array(1000, 10000, 100000, 1000000, 10000000)
//    // numberOfPairs = 100 // legge dalle properties
//    savedTask = new Saver.SaveClosure()
//
//    for( bl <- oom) {
//      for (j <- -1 to 1) {
//
//        seqLen = bl + (bl * j / 10)
//
//        var p = Seq[Any](numberOfPairs, seqLen, uniformDist, appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"))
//
//        savedTask.save(buildDatasetWithDistributionEx, p)
//      }
//    }
//
//    for( bl <- oom) {
//      for (j <- -1 to 1) {
//
//        seqLen = bl + (bl * j / 10)
//        for (g <- gValues) {
//
//          // generiamo entrambi i dataset dipendenti dal nullModel Uniform
//          var p = Seq[Any]( seqLen, motif, g,
//            appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"),
//            appProperties.getProperty("powerstatistics.datasetBuilder.altMotifPrefix"))
//          savedTask.save(buildDatasetMotifReplaceEx, p)
//
//          p = Seq[Any]( seqLen, g,
//            appProperties.getProperty("powerstatistics.datasetBuilder.uniformPrefix"),
//            appProperties.getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix"))
//          savedTask.save(buildDatasetPatternTransferEx, p)
//        }
//      }
//    }
//  }


  def buildDatasetWithShuffle(numberOfPairs: Int, sequenceLen: Int,
                              geneSize: Int, prefix: String): Unit = {

    for( i <- 1 to numberOfPairs) {

      val seq1 = buildShuffledSequence( sequenceLen, geneSize)
      val seq2 = buildShuffledSequence( sequenceLen, geneSize)

      if (debug) {
        val st1 = getDistribution(seq1)
        print(s"SEQ1(${i},${prefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println

        val st2 = getDistribution(seq2)
        print(s"SEQ2(${i},${prefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println
      }

      val outputPath = getNullModelFilename( prefix, sequenceLen, i)
      saveSequencePair( outputPath, getSequenceHeader( prefix, i, 0.0, 1), seq1,
                                    getSequenceHeader( prefix, i, 0.0, 2), seq2)
    }
    closeOutput()
  }

  // build a null model sequence randomly shuffling a real sequence powerstatistics.datasetBuilder.escherichiaColi
  def buildShuffledSequence( sequenceLen: Int, geneSize: Int): Array[Char] = {

    val seq: Array[Char] = new Array[Char](sequenceLen)
    // initialize the random source
    var len = 0
    val inSeq = appProperties.getProperty("powerstatistics.datasetBuilder.escherichiaColi")
    val max = inSeq.length() / geneSize
    var i = 0 : Int
    var ndx = 0 : Int

    while(len < sequenceLen) {
      ndx = rng.nextInt(max) * geneSize
      breakable {
        for (i <- 0 until geneSize) {
          val dst = len + i
          if (dst >= sequenceLen)
            break
          else
            seq(dst) = inSeq.charAt(ndx + i)
        }
      }
      len = len + geneSize
    }
    return seq
  }


  def buildDatasetWithDistribution(numberOfPairs: Int, sequenceLen: Int,
                                   distribution: Array[Double], prefix: String): Unit = {

    val seq1: Array[Char] = new Array[Char](sequenceLen)
    val seq2: Array[Char] = new Array[Char](sequenceLen)
    var st1: Array[Int] = null // Array.fill(4)(0)
    var st2: Array[Int] = null // Array.fill(4)(0)
    var n1 = 0
    var n2 = 0

    val rg = new randomNumberGenerator(distribution)

    for( i <- 1 to numberOfPairs) {

      for (c <- 0 to sequenceLen - 1) {
        n1 = rg.getNextBase()
        seq1(c) = nucleotideRepr(n1)

        n2 = rg.getNextBase()
        seq2(c) = nucleotideRepr(n2)
      }

      if (debug) {
        st1 = getDistribution(seq1)
        print(s"SEQ1(${i},${prefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println

        st2 = getDistribution(seq2)
        print(s"SEQ2(${i},${prefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println
      }

      val outputPath = getNullModelFilename( prefix, sequenceLen, i)
      saveSequencePair( outputPath, getSequenceHeader( prefix, i, 0.0, 1), seq1,
                                    getSequenceHeader( prefix, i, 0.0, 2), seq2)
    }
    closeOutput()
  }



  def buildDatasetMotifReplace(sequenceLen: Int, motif: Array[Char], probSubstitution: Double,
                               inPrefix: String, outPrefix: String): Unit = {

    val rg = new randomNumberGenerator(probSubstitution)

    // motifLen variabile in funzione della lunghezza della sequenza
    // val motifLen = kValueFromSeqlen(sequenceLen)
    // in alternativa motifLen costante >= del valore massimo di k utilizzato.
    val motifLen = parsed.getOrElse("patternSize", appProperties.getProperty("powerstatistics.datasetBuilder.replacePatternLength").toInt)

    if (debug) println(f"Motif Replace, gamma=$probSubstitution%.3f\n")

    // per tutte le coppie nel file input
    for( i <- 1 to numberOfPairs) {
      val sq = readSequencePair(inPrefix, sequenceLen, i)
      val seq1Header = sq._1
      if (debug) println(s"header: ${seq1Header.mkString("")}")
      val seq1Body = sq._2
      val seq2Header = sq._3
      if (debug) println(s"header: ${seq2Header.mkString("")}")
      val seq2Body = sq._4

      val ris = replaceMotif(seq1Body, seq2Body, motifLen, motif, rg)

      val AMseq1 = ris._1
      val AMseq2 = ris._2

      if (debug) {
        val st1 = getDistribution(AMseq1)
        print(s"SEQ1(${i},${outPrefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println

        val st2 = getDistribution(AMseq2)
        print(s"SEQ2(${i},${outPrefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println
      }
      val outputPath = getAltModelFilename( outPrefix, sequenceLen, i, probSubstitution, inPrefix)
      saveSequencePair( outputPath, getSequenceHeader(outPrefix, i, probSubstitution, 1), AMseq1,
                                    getSequenceHeader(outPrefix, i, probSubstitution, 2), AMseq2)
    }
    closeInput()
    closeOutput()
  }



  def replaceMotif( inputSeq1: Array[Char], inputSeq2: Array[Char], motifLen: Int, motif: Array[Char],
                    rg: randomNumberGenerator) : (Array[Char], Array[Char]) = {

    var bd: Boolean = false
    var c: Int = 0
    var cnt: Int = 0
    var k : Char = '\0'

    // Randomly selects a motif from a set of motifs to keep the number of motif replacements
    // equal to the P(gamma) * 25.000 / 5 (Reinert's dataset has MotifLength= 5 e SeqLength= 25.000)
    // val maxNumberOfMotifs : Int =  max( 1, (inputSeq1.length / (5000 * motifLen)).toInt)  // SeqLen 10.000.000 => circa 200
    var maxNumberOfMotifs : Int =  max( 1, (inputSeq1.length / 25000 ).toInt)  // SeqLen 10.000.000 => circa 200
    maxNumberOfMotifs = min( maxNumberOfMotifs, motif.length/motifLen) // to avoid arrayOutOfBoundError on motif source string
    val motifCnt = Array.fill[Int](maxNumberOfMotifs)(0)        // to count effective use of each motif
    var motifIndex : Int = 0

    while( c < inputSeq1.length - motifLen) {
      // distribuzione di bernoulli
      bd = rg.getNextBernoulliValue()
      if (bd) {
        // sceglie un motif a caso di lunghezza motifLen dalla sequenza motif
        motifIndex = if (maxNumberOfMotifs <= 1) 0 else rg.rng.nextInt(maxNumberOfMotifs)
        val ndx = motifIndex * motifLen
        motifCnt(motifIndex) = motifCnt(motifIndex) + 1
        // replace the characters specified by motif (c is incremented no overlap)
        for (j <- ndx to ndx + motifLen - 1) {
          k = motif(j)
          inputSeq1(c) = k
          inputSeq2(c) = k
          c = c + 1
        }
        cnt = cnt + 1
      }
      else
        c = c + 1 // jump to next base
    }
    if (debug) {
      println(s"${cnt}/${inputSeq1.length} motif substitutions from ${maxNumberOfMotifs} motifs")
      for( i <- 0 until motifCnt.length)
        println(s"${i} -> ${motifCnt(i)}")
    }

    return (inputSeq1, inputSeq2)
  }



  def buildDatasetPatternTransfer(sequenceLen: Int, probSubstitution: Double,
                                  inPrefix: String, outPrefix: String): Unit = {

    val rg = new randomNumberGenerator(probSubstitution)


    // len variabile in funzione della lunghezza della sequennza
    // val patternLen = kValueFromSeqlen( sequenceLen)
    // in alternativa si puo' usare un valore costante per tutte le lunghezze >= del massimo valore di k utilizzato nei test
    val patternLen = parsed.getOrElse("patternSize", appProperties.getProperty("powerstatistics.datasetBuilder.replacePatternLength").toInt)

    if (debug) println(f"Pattern Transfer, gamma = $probSubstitution%.3f\n")
    // per tutte le coppie
    for( i <- 1 to numberOfPairs) {

      val sq = readSequencePair(inPrefix, sequenceLen, i)
      val seq1Header = sq._1
      if (debug) println(s"header: ${seq1Header.mkString("")}")
      val seq1Body = sq._2
      val seq2Header = sq._3
      if (debug) println(s"header: ${seq2Header.mkString("")}")
      val seq2Body = sq._4

      val AMseq2 = patternTransfer(seq1Body, seq2Body, patternLen, rg)

      if (debug) {
        val st1 = getDistribution(seq1Body)
        print(s"SEQ1(${i},${outPrefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println

        val st2 = getDistribution(AMseq2)
        print(s"SEQ2(${i},${outPrefix}): ")
        for (c <- 0 to 3)
          printf("%c=%.1f%%, ", nucleotideRepr(c), st1(c) * 100.toDouble / sequenceLen)
        println
      }
      val outputPath = getAltModelFilename( outPrefix, sequenceLen, i, probSubstitution, inPrefix)
      saveSequencePair( outputPath, getSequenceHeader(outPrefix, i, probSubstitution, 1), seq1Body,
        getSequenceHeader(outPrefix, i, probSubstitution, 2), AMseq2)
    }
    closeInput()
    closeOutput()
  }



  def patternTransfer( inputSeq: Array[Char], trgtSeq: Array[Char], patternLen: Int,
                       rg: randomNumberGenerator) : Array[Char] = {

    var bd: Boolean = false
    var c: Int = 0
    var cnt: Int = 0

    while( c < trgtSeq.length - patternLen) {
      // distribuzione di bernoulli
      bd = rg.getNextBernoulliValue()
      if (bd) {
        if (debug) printf("%d, ", c)
        // replace the characters extracted from the inputSeq in the targetSeq (c is incremented no overlap)
        for (k <- 1 to patternLen) {
          trgtSeq(c) = inputSeq(c)
          c = c + 1
        }
        cnt = cnt + 1
      }
      else
        c = c + 1 // jump to next base
    }
    if (debug)  println(s"${cnt}/${trgtSeq.length} substitutions")
    return trgtSeq
  }


  def getDistribution(seq: Array[Char]) : Array[Int] = {
    var dist: Array[Int] = Array.fill(5)(0)

    for( b <- seq) {
      b match {
        case 'A' => dist(0) = dist(0) +1
        case 'C' => dist(1) = dist(1) +1
        case 'G' => dist(2) = dist(2) +1
        case 'T' => dist(3) = dist(3) +1
        case 'N' => dist(4) = dist(4) +1
        case x => print(s"unknown character: ${x}")
      }
    }
    return dist
  }



  //
  // legge una coppia di sequenze
  //
  def readSequencePair(prefix: String, sequenceLen: Int, seqId : Int) :
                            (Array[Char], Array[Char], Array[Char], Array[Char]) = {

    var retVal : (Array[Char], Array[Char], Array[Char], Array[Char]) = null
    val inputPath: String = getNullModelFilename(prefix, sequenceLen, seqId)

    try {
      if (multifile || reader == null) {
        if (local) {
          // legge dal filesystem locale
          reader = scala.io.Source.fromFile(inputPath)
        }
        else {
          // Questo codice legge (o scrive) dal HDFS.
          // Source.fromFile legge solo da file system locali
          // val hdfs = FileSystem.get(new URI("hdfs://master:8020/"), new Configuration())
          val hdfs = FileSystem.get(sc.hadoopConfiguration)
          val path = new Path(inputPath)
          val stream = hdfs.open(path)
          reader = scala.io.Source.fromInputStream(stream)
        }
        inputIter = reader.getLines()
      }
    }
    catch {
      // Case statement-1
      case x: FileNotFoundException => {
        println(s"Exception: ${inputPath} File missing")
      }
      case x: IOException   => {
        println("Input/output Exception")
      }
    }
    if (inputIter.hasNext) {
      val s1 = inputIter.next
      val seq1Header = new Array[Char](s1.length)
      s1.copyToArray(seq1Header)
      val seq1Body = new Array[Char](sequenceLen)
      inputIter.next.copyToArray(seq1Body)

      val s2 = inputIter.next
      val seq2Header = new Array[Char](s2.length)
      s2.copyToArray(seq2Header)
      val seq2Body = new Array[Char](sequenceLen)
      inputIter.next.copyToArray(seq2Body)

      retVal = (seq1Header, seq1Body, seq2Header, seq2Body)
      }

    if (multifile)
      closeInput()

    return retVal
  }


  def saveSequencePair(outputPath: String, sequenceHeader1: String, seq1: Array[Char],
                                            sequenceHeader2: String, seq2: Array[Char]): Unit = {

    if (multifile || writer == null) {
      openOutput( outputPath)
    }
    writer.write( s">${sequenceHeader1}\n")
    writer.write( seq1.mkString)
    writer.newLine()

    writer.write( s">${sequenceHeader2}\n")
    writer.write( seq2.mkString)
    writer.newLine()

    if (multifile) {
      closeOutput
    }
  }


  def openOutput( outputPath: String) {
    writer = new BufferedWriter(
      new OutputStreamWriter(FileSystem.get(URI.create(outputPath),
        new Configuration()).create(new Path(outputPath))))
  }


  def closeOutput() : Unit = {
    if (writer != null) {
      writer.close
      writer = null
    }
  }


  def closeInput() : Unit = {
    if (reader != null) {
      reader.close
      reader = null
    }
  }


  def getSequenceHeader(prefix: String, ndx: Int, g: Double, c: Int) : String = {

    val pair = if (c == 1) 'A' else 'B'
    val pg = if (g > 0) ".G=" + "%.3f".format( g).substring(2) else ""

    val name = "%s.%05d%s-%c".format(prefix, ndx, pg, pair)

    return name
  }


  def getNullModelFilename(dataSetName: String, seqLength: Int, seqId: Int) : String = {

    val fld2 =  if (multifile) seqId else numberOfPairs
    val outputPath =  f"$savePath%s/$dataSetName%s-$fld2%04d.$seqLength%d.$fileExt%s"

    return outputPath
  }



  def  getAltModelFilename(alternateModel: String, seqLength: Int, seqId: Int, g: Double, nullModel: String) : String = {

    val fld2 =  if (multifile) seqId else numberOfPairs
    val pg = ".G=%.3f".format( g).replace(',','.')
    var sfx: String = null

    if (nullModel.charAt(0) == 'U')
      // from Uniform
        sfx = nullModel.substring(0,1)
      else if (nullModel.substring(0, 4).compareTo("rnd_") == 0)
        // rnd_mito e rnd_shigella ...
        sfx = nullModel.substring(4,6)
      else
        // GCReach
        sfx = nullModel.substring(0,2)

    val dataset = alternateModel + "-" + sfx

    val outputPath = f"$savePath%s/$dataset%s-$fld2%04d.$seqLength%d$pg%s.$fileExt%s"
    return outputPath
  }


  def getMotif( motif: String) : Array[Char] = {

    val res = new Array[Char]( motif.length)
    var i = 0
    motif.foreach( c => {
      res(i) = c match {
        case 'A' | 'a' => 'A'
        case 'C' | 'c' => 'C'
        case 'G' | 'g' => 'G'
        case 'T' | 't' => 'T'
        case invChar => System.err.println(s"Wrong motif specification: ${motif} ${invChar}")
                        throw new IllegalArgumentException(s"illegal base in motif specificaton. ${invChar} ")
      }
      i = i + 1
    })
    return res
  }


  def kValueFromSeqlen( seqLen: Int) : Int = {
    //  le n sequenze hanno tutte la stessa lunghezza => somma(seqLen(i)) / n = seqLen
    val k: Int = (log10(seqLen.toDouble) / log10(4.toDouble)).ceil.toInt - 1
    return k
  }
}

object DatasetTypes extends Enumeration
{
  type DatasetType = Value

  val Uniform, GCReach, AlternateMotifReplace, AlternatePatternTransfer = Value
}
