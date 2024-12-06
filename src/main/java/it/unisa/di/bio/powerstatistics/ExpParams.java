package it.unisa.di.bio.powerstatistics;

import java.util.Comparator;

public class ExpParams {

    int k = 0;
    int numPairs = 0;
    int seqLength = 0;
    double gamma = 0.;
    double minTr = 0.;
    double maxTr = 0.;
    double varTr = 0.;

    String inputFolder = "";
    String filename = "";

    public ExpParams() {}

    public ExpParams(int numPairs, int seqLength) {
        this.numPairs = numPairs;
        this.seqLength = seqLength;
    }

    static class SortBySequenceLegth implements Comparator<ExpParams>
    {
        public int compare(ExpParams a, ExpParams b)
        {
            if (a.seqLength < b.seqLength)
                return -1;
            else if (a.seqLength > b.seqLength)
                return 1;
            else
                return 0;
        }
    }
}
