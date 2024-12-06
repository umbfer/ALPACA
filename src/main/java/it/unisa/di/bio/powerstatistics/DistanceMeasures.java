package it.unisa.di.bio.powerstatistics;

import java.util.HashMap;


public class DistanceMeasures {

    private HashMap<String, Boolean> distances = new HashMap<String, Boolean>();

    private static DistanceMeasures instance;

    private DistanceMeasures() {

        for( Measures x : Measures.values()) {
            distances.put( x.toString().toLowerCase(), true);
        }
        distances.put(Measures.Canberra.toString().toLowerCase(), true);        // Distance
        distances.put(Measures.Chebyshev.toString().toLowerCase(), true);       // Distance
        distances.put(Measures.ChiSquare.toString().toLowerCase(), true);       // Distance
        distances.put(Measures.D2.toString().toLowerCase(), false);             // Similarity
        distances.put(Measures.D2s.toString().toLowerCase(), false);            // Similarity
        distances.put(Measures.D2star.toString().toLowerCase(), false);         // Similarity
        distances.put(Measures.D2z.toString().toLowerCase(), false);            // Similarity
        distances.put(Measures.Euclidean.toString().toLowerCase(), true);       // Distance
        distances.put(Measures.HarmonicMean.toString().toLowerCase(), false);   // Similarity
        distances.put(Measures.Intersection.toString().toLowerCase(), false);   // Similarity
        distances.put(Measures.Jaccard.toString().toLowerCase(), false);        // Similarity
        distances.put(Measures.Jeffrey.toString().toLowerCase(), true);         // Distance
        distances.put(Measures.JensenShannon.toString().toLowerCase(), true);   // Distance
        distances.put(Measures.Kulczynski2.toString().toLowerCase(), false);    // Similarity
        distances.put(Measures.Manhattan.toString().toLowerCase(), true);       // Distance
        distances.put(Measures.SquaredChord.toString().toLowerCase(), true);    // Distance
        distances.put(Measures.Mash.toString().toLowerCase(), true);            // Distance
    }


    public static DistanceMeasures getInstance(){
        if(instance == null){
            instance = new DistanceMeasures();
        }
        return instance;
    }

    public String getMeasureType(String measure) {
        return distances.get(measure.toLowerCase()) ? "distance" : "similarity";
    }

    public boolean isDistance( String measure) {

        return distances.get( measure.toLowerCase());
    }

    //
    public enum Measures {
        Canberra, Chebyshev, ChiSquare,
        D2, D2s, D2star, D2z,
        Euclidean, HarmonicMean, Intersection,
        Jaccard, Jeffrey, JensenShannon,
        Kulczynski2, Manhattan, SquaredChord, Mash
    }


//    public enum Measures {
//        CanberraDistance,  ChebyshevDistance,
//        ChiSquareDistance, D2Distance, D2zDistance, EuclideanDistance,
//        HammingDistance, HarmonicMeanDistance, IntersectionDistance, JaccardDistance,
//        JeffreysDistance, JensenShannonDistance, Kulczynski1Distance, Kulczynski2Distance,
//        Kulczynski2NoPseudoDistance, ManhattanDistance, SquaredChordDistance, Mash
//    }
}
