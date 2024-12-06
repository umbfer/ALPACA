package it.unisa.di.bio



import scala.collection.mutable.ListBuffer

object Saver {

  class SaveClosure {
    // val closures: ListBuffer[(f: Any* => Unit, params: Seq[Any])]
    var closures = ListBuffer[((Seq[Any]) => Unit, Seq[Any])]()

    // a method that takes a function and a string, and passes the string into
    // the function, and then executes the function
    def exec(f:(Seq[Any]) => Unit, params: Seq[Any]) : Int = {
      f( params)
      return 1
    }

    def save(f: (Seq[Any]) => Unit, params: Seq[Any]): Unit = {
      val c = (f, params)
      closures += c
    }

    def printClosure() = {
      for( t <- closures) {
        printf(s"functor: ${t._1}, parameters: ${t._2}")
      }
    }
  }
}
