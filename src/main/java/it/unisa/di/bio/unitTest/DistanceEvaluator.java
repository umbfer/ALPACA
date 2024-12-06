package it.unisa.di.bio.unitTest;

import it.unisa.di.bio.unitTest.KMersStatistics;

public abstract class DistanceEvaluator {

    protected KMersStatistics histo1;
    protected KMersStatistics histo2;
    protected String distanceName;


    public DistanceEvaluator( KMersStatistics h1, KMersStatistics h2) {
        histo1 = h1;
        histo2 = h2;
    }

    public abstract double evaluateDistance();
    public String getDistanceName() {return distanceName;}
}
