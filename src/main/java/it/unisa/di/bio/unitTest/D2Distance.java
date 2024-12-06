package it.unisa.di.bio.unitTest;

import java.util.Enumeration;

public class D2Distance extends DistanceEvaluator {

    public D2Distance(KMersStatistics h1, KMersStatistics h2) {
        super( h1, h2);
        distanceName = "D2Distance";
    }

    @Override
    public double evaluateDistance() {

        double value = 0.;

        Enumeration<String> kmerms1 = histo1.histogram.keys();

        while (kmerms1.hasMoreElements()) {
            String kmer = kmerms1.nextElement();
            int v1 = histo1.histogram.get(kmer).intValue();
            Integer v2 = histo2.histogram.get(kmer);
            if (v2 != null) {
                value += v1 * v2.intValue();
                System.out.printf("kmer=%s -> v1=%d, v2=%d\n", kmer, v1, v2.intValue());
            }
        }
        // final step
        return value;
    }
}
