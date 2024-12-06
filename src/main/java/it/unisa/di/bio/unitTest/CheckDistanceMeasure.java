package it.unisa.di.bio.unitTest;

public class CheckDistanceMeasure {

    public static void main(String[] args) {

        if (args.length != 3) {
            System.out.println("Usage:\n" + "DistanceMeasures path file1 file2" );
            System.exit(-1);
        }
        String resultsPath = args[0];
        String seq1Filename = resultsPath + "/" + args[1];
        String seq2Filename = resultsPath + "/" + args[2];

        KMersStatistics h1 = new KMersStatistics( seq1Filename);
        KMersStatistics h2 = new KMersStatistics( seq2Filename);

        DistanceEvaluator d2 = new D2Distance(h1, h2);

        System.out.printf("%s value = %.3f\n",
                d2.getDistanceName(), d2.evaluateDistance());
    }
}