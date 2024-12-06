package it.unisa.di.bio.powerstatistics;


import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;


public class ReadDistances {

    String[]    measureNames;
    private int refMeasure;
    DatasetType nullModel;
    DatasetType alternateModel;

    ArrayList<DistanceValues> nullModelValues;
    ArrayList<DistanceValues> alternateModelValues;
    ArrayList<DistanceValues> type1CheckValues;

    PowerProps  appProps;

    String      inPath;
    String      fnamePrefix = "dist-"; // "distances" nella precedente versione con i dataset di Reinert
    String      fileExt = ".csv";
    int         numberOfPairs = 0;
    int         seqLen;
    int         k;
    double      gamma;

    double      threshold = 0;

    DistanceMeasures distances = null;


    public ReadDistances(DatasetType nullModel, DatasetType alternateModel, String inPath,
                         int numberOfPairs, int k, int seqLen, double gamma){

        this.nullModel      = nullModel;
        this.alternateModel = alternateModel;
        this.inPath         = inPath;
        this.numberOfPairs  = numberOfPairs;
        this.k              = k;
        this.seqLen         = seqLen;
        this.gamma          = gamma;

        this.refMeasure = -1;

        appProps = PowerProps.getInstance();

        distances = DistanceMeasures.getInstance();

        readNullModel();
    }



    public ReadDistances(DatasetType nullModel, DatasetType alternateModel, String filename, double gamma){
        this.nullModel      = nullModel;
        this.alternateModel = alternateModel;

        this.refMeasure = -1;

        appProps = PowerProps.getInstance();
        distances = DistanceMeasures.getInstance();

        ExpParams p = parseFilename( filename);
        this.inPath         = p.inputFolder;
        this.numberOfPairs  = p.numPairs;
        this.k              = p.k;
        this.seqLen         = p.seqLength;
        this.gamma          = gamma;

        readNullModel( filename);
    }



    public int getRefMeasure() {
        return refMeasure;
    }


    public int setRefMeasure(int refMeasure) {
        if (refMeasure < measureNames.length)
            this.refMeasure = refMeasure;

        // ordina i risultati rispetto alla misura scelta
        sortValues(nullModelValues);
        sortValues(alternateModelValues);
        sortValues(type1CheckValues);

        return this.refMeasure;
    }



    public void printMeasures() {
        for(int i = 0; i < measureNames.length; i++) {
            System.out.printf("%s (%s)%s", measureNames[i], distances.getMeasureType(measureNames[i]),
                    ((i == measureNames.length - 1) ? "\n" : ", "));
        }
    }

    public String getCurrentMeasureName() {
        if (refMeasure >= 0 && refMeasure < measureNames.length)
            return measureNames[refMeasure];
        else
            return null;
    }


    public double getThreshold() {
        return threshold;
    }

    public void setThreshold(double alfa) {
        int ndx = (int) (nullModelValues.size() * alfa);

        this.threshold = nullModelValues.get( ndx).v[getRefMeasure()];
        return;
    }


    private void sortValues(ArrayList<DistanceValues> lst) {

        Comparator<DistanceValues> cmp = null;   // dal migliore al peggiore
        if (distances.isDistance(measureNames[refMeasure]))
            cmp = new SortByDistance();     // distanza => dal piu' piccolo al piu grande
        else
            cmp = new SortByDistanceRev();  // similarità => dal piu' grande al  piu' piccola

        lst.sort(cmp);
    }



    public int readNullModel() {
        String filename = getNullModelFileName( nullModel);

        return readNullModel( filename);
    }

    public int readNullModel(String filename) {

        try
        {
           // Prepare list.
            Stream<String> s = Files.lines(Paths.get( filename));
            int initialCapacity = (int) s.count();
            s.close();
            nullModelValues = new ArrayList<DistanceValues>( initialCapacity );

            BufferedReader reader = Files.newBufferedReader( Paths.get( filename));

            Iterable<CSVRecord> records = CSVFormat.RFC4180.withFirstRecordAsHeader().parse( reader);

            for (CSVRecord record : records) {

                if (record.getRecordNumber() == 1) {
                    int nm = record.getParser().getHeaderNames().size() - 2; // meno i descrittori delle sequenze
                    measureNames = new String[nm];

                    for(int i = 0; i < nm; i++) {
                        // qui vengono settate per tutta la classe gli altri controllano che siano allineate
                        measureNames[i] = record.getParser().getHeaderNames().get(i+2);
                    }
                }
                DistanceValues rec = new DistanceValues();

                Iterator<String> it = record.iterator();
                int i = 0;
                while( it.hasNext()) {
                    String val = it.next();
                    // System.out.printf("field %d : %s\n", i, val);
                    switch (i++) {
                        case 0:
                            rec.testId = val.substring(0, val.length() - 2);
                            break;
                        case 1:
                            if (rec.testId.compareTo( val.substring(0, val.length() - 2)) != 0)
                                System.out.printf("Wrong sequenceId %s vs %s\nIgnored.\n", rec.testId, val);
                            break;
                        default:
                            rec.v[i - 3] = Double.parseDouble(val);
                    }
                }
                nullModelValues.add( rec );
            }
            if (nullModelValues.size() != numberOfPairs) {
                System.out.printf("wrong number of results: expected:%d, read:%d\n", numberOfPairs, nullModelValues.size());
            }

            reader.close();
            return nullModelValues.size();
        }
        catch ( IOException e ) {
            e.printStackTrace();
        }
        return 0;
    }



    public int readAlternateModel() {

        String filename = getAlternateModelFileName(nullModel, alternateModel);

        return readAlternateModel(filename);
    }


    public int readAlternateModel(String filename) {

        try
        {
            // Prepare list.
            Stream<String> s = Files.lines(Paths.get( filename));
            int initialCapacity = (int) s.count();
            s.close();
            alternateModelValues = new ArrayList<DistanceValues>( initialCapacity );

            BufferedReader reader = Files.newBufferedReader( Paths.get( filename));

            Iterable<CSVRecord> records = CSVFormat.RFC4180.withFirstRecordAsHeader().parse( reader);

            for (CSVRecord record : records) {

                if (record.getRecordNumber() == 1) {
                    int nm = record.getParser().getHeaderNames().size() - 2; // meno i descrittori delle sequenze
                    for(int i = 0; i < nm; i++) {
                        if (measureNames[i].compareTo(record.getParser().getHeaderNames().get(i + 2)) != 0) {
                            System.out.printf("Misalligned results %s vs %s (i = %d)\n", measureNames[i],
                                    record.getParser().getHeaderNames().get(i + 2));
                        }
                    }
                }
                DistanceValues rec = new DistanceValues();

                Iterator<String> it = record.iterator();
                int i = 0;
                while( it.hasNext()) {
                    String val = it.next();
                    // System.out.printf("field %d : %s\n", i, val);
                    i++;
                    if (i <= 2) {
                        rec.testId += val;
                    }
                    else if (i > 2) {
                        rec.v[i - 3] = Double.parseDouble(val);
                        if (Double.isNaN(rec.v[i - 3]))
                            rec.v[i - 3] = Double.MAX_VALUE;
                    }
                }
                alternateModelValues.add( rec );
            }
            if (alternateModelValues.size() != numberOfPairs) {
                System.out.printf("wrong number of results: expected:%d, read:%d\n", numberOfPairs, alternateModelValues.size());
            }
            // sortAlternateModelValues(); rispetto a quale misura ??? verrà ordinato quando si chiama setMeasure()
            reader.close();
            return alternateModelValues.size();
        }
        catch ( IOException e ) {
            e.printStackTrace();
        }
        return 0;
    }




    public double getPowerValue() {

        int tot = 0;

        if (distances.isDistance(measureNames[refMeasure])) {
            for (DistanceValues dist : alternateModelValues) {
                if (dist.v[refMeasure] < threshold)     // distanza => risultato migliore della soglia
                    tot++;
            }
        }
        else {
            for (DistanceValues dist : alternateModelValues) {
                if (dist.v[refMeasure] > threshold)     // similarità => risultato migliore della soglia
                    tot++;
            }
        }
        return (double) tot / (double) numberOfPairs;
    }




    public int readType1Check() {

        String filename = getType1FileName(nullModel);

        return readType1Check(filename);
    }


    public int readType1Check(String filename) {

        try
        {
            // Prepare list.
            Stream<String> s = Files.lines(Paths.get( filename));
            int initialCapacity = (int) s.count();
            s.close();
            type1CheckValues = new ArrayList<DistanceValues>( initialCapacity );

            BufferedReader reader = Files.newBufferedReader( Paths.get( filename));

            Iterable<CSVRecord> records = CSVFormat.RFC4180.withFirstRecordAsHeader().parse( reader);

            for (CSVRecord record : records) {

                if (record.getRecordNumber() == 1) {
                    int nm = record.getParser().getHeaderNames().size() - 2; // meno i descrittori delle sequenze
                    for (int i = 0; i < nm; i++) {
                        if (measureNames[i].compareTo(record.getParser().getHeaderNames().get(i + 2)) != 0) {
                            System.out.printf("Misallined results %s vs %s (i = %d)\n", measureNames[i],
                                    record.getParser().getHeaderNames().get(i + 2));
                        }
                    }
                }
                DistanceValues rec = new DistanceValues();

                Iterator<String> it = record.iterator();
                int i = 0;
                while( it.hasNext()) {
                    String val = it.next();
                    // System.out.printf("field %d : %s\n", i, val);
                    i++;
                    if (i <= 2) {
                        rec.testId += val;
                    }
                    else if (i > 2) {
                        rec.v[i - 3] = Double.parseDouble(val);
                        if (Double.isNaN(rec.v[i - 3]))
                            rec.v[i - 3] = Double.MAX_VALUE;
                    }
                }
                type1CheckValues.add( rec );
            }
            if (type1CheckValues.size() != numberOfPairs) {
                System.out.printf("wrong number of results: expected:%d, read:%d\n", numberOfPairs, type1CheckValues.size());
            }
            // sortType1Values(); rispetto a quale misura ??? verrà ordinato quando si chiama setMeasure()
            reader.close();
            return type1CheckValues.size();
        }
        catch ( IOException e ) {
            e.printStackTrace();
        }
        return 0;
    }




    public double getType1Value() {

        int tot = 0;

        if (distances.isDistance(measureNames[refMeasure])) {
            for (DistanceValues dist : type1CheckValues) {
                if (dist.v[refMeasure] < threshold)     // distanza => risultato inferiore alla soglia
                    tot++;
            }
        }
        else {
            for (DistanceValues dist : type1CheckValues) {
                if (dist.v[refMeasure] > threshold)     // similarità => risultato superiore alla soglia
                    tot++;
            }
        }
        return (double) tot / (double) numberOfPairs;
    }


    private String getNullModelFileName( DatasetType nullModel) {

        String fname = nullModel.toString();

        if (fname.length() == 0) {
            System.err.println("Unknown null model");
            return null;
        }
        fname = String.format("%s/%sk=%d_%s-%d.%d%s", inPath, fnamePrefix, k, fname, numberOfPairs, seqLen, fileExt);
        return fname;
    }


    private String getAlternateModelFileName( DatasetType nullModel, DatasetType alternateModel) {

        String nm = "";

        switch (nullModel) {
            case Uniform: nm = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.uniformPrefix").substring(0,1);
                    break;
            case GCReach: nm = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.GCReachPrefix").substring(0,2);
                    break;

            default:
                nm = nullModel.toString().substring(0, 2).toLowerCase();
        }

        String fname = "";

        switch (alternateModel) {
            case MotifReplace: fname = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.altMotifPrefix");
                break;
            case PatternTransfer: fname = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.altPatTransfPrefix");
                break;
            default:
                System.err.println("Unknown alternate model");
        }
        String gf = String.format("G=%.3f", gamma).replace(",", ".");
        fname = String.format("%s/%sk=%d_%s-%s-%d.%d.%s%s", inPath, fnamePrefix, k, fname, nm,
                numberOfPairs, seqLen, gf, fileExt);
        return fname;
    }


    private String getType1FileName( DatasetType nullModel) {

        String src = nullModel.toString();

        String gf = String.format("G=%.3f", gamma).replace(",", ".");
        String fname = String.format("%s/%sk=%d_%s%s-%d.%d%s", inPath, fnamePrefix, k,
                src, appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.Type1CheckSuffix"),
                numberOfPairs, seqLen, fileExt);
        return fname;
    }



    static public ExpParams parseFilename( String filename) {

        ExpParams param = new ExpParams();
        Path path = Paths.get(filename);
        param.inputFolder = path.getParent().toString();
        param.filename = filename;
        // distancesk=6_Uniform-10000.900.fasta.csv
        // Pattern p = Pattern.compile(".*k=(\\d+)_\\S-(\\d+)\\.(\\d+)\\.fasta\\.csv");
        Pattern pattern = Pattern.compile(".*k=(\\d+).*-(\\d+)\\.(\\d+).*$");
        Matcher match = pattern.matcher(filename);
        if (!match.find()) {
            System.out.printf("Malformed csv file name %s\n", filename);
            return null;
        }
        param.k = Integer.parseInt(match.group(1));
        param.numPairs = Integer.parseInt( match.group(2));
        param.seqLength = Integer.parseInt(match.group(3));

        return param;
    }


    class DistanceValues {
        String testId;
        double[] v;

        public DistanceValues() {
            testId = "";
            v = new double[DistanceMeasures.Measures.values().length];
        }
    }


    class SortByDistance implements Comparator<DistanceValues>
    {
        public int compare(DistanceValues a, DistanceValues b)
        {
            if (a.v[refMeasure] > b.v[refMeasure])
                return 1;
            else if (a.v[refMeasure] < b.v[refMeasure])
                return -1;
            else
                return 0;
        }
    }

    class SortByDistanceRev implements Comparator<DistanceValues>
    {
        public int compare(DistanceValues a, DistanceValues b)
        {
            if (a.v[refMeasure] < b.v[refMeasure])
                return 1;
            else if (a.v[refMeasure] > b.v[refMeasure])
                return -1;
            else
                return 0;
        }
    }

    enum Order{ asc, desc};

}
