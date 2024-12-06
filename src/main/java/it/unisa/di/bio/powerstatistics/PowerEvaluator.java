package it.unisa.di.bio.powerstatistics;


import java.io.*;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.util.ArrayList;


import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import static it.unisa.di.bio.powerstatistics.DatasetType.*;



public class PowerEvaluator {

    static String resultsPath = "";

    static PowerProps appProps = null;


    static DatasetType nullModel = Uniform;
    // static DatasetType nullModel = GCReach;
    static DatasetType alternateModel = PatternTransfer;
    // static DatasetType alternateModel = MotifReplace;

    private static final String AREAS_CSV_FILE = "MeasureAreas";
    private static String Current_CSV_File = "";

    public static void main(String[] args) {

        resultsPath = args[0];
        // alternateModel = DatasetType.valueOf(args[1]);
        // int k = Integer.parseInt(args[1]);
        // double alfa = Double.parseDouble(args[2]);

        appProps = PowerProps.getInstance();

        System.out.printf("Reading results from dir: %s\n", resultsPath);
        // AnalyzeDataset1( k, alfa);
        // AnalyzeDataset2( alfa);

        double[] alphaValues = {0.010, 0.050, 0.100};
        DatasetType[] dtValues;
        double[] av = new double[1];
        if (args[1].compareTo("Both") == 0) {
            dtValues = new DatasetType[2];
            dtValues[0] = MotifReplace;
            dtValues[1] = PatternTransfer;
        }
        else {
            dtValues = new DatasetType[1];
            dtValues[0] = DatasetType.valueOf(args[1]);
        }

        for(DatasetType m : dtValues) {
            alternateModel = m;
            for (double v : alphaValues) {
                switch(args[2])
                {
                    case "syntheticAllK":       // Dataset4-1000
                        AnalyzeDataset4(v, Uniform);
                        break;
                    case "syntheticAllLen":     // ultimo data set distanza lineare tra le lunghezze
                        AnalyzeDataset5(v, Uniform, 0, -1); // 0, -1
                        break;
                    case "syntheticShortLen":     // ultimo data set distanza lineare tra le lunghezze
                        AnalyzeDataset6(v, Uniform, 0, -1); // 0, -1
                        break;
                    default:
                        AnalyzeDataset3(v, DatasetType.valueOf(args[2]));
                }
            }
        }
    }


    // Analisi del dataset small len in {1000, 5000, 10000, 50000, 100000}
    static void AnalyzeDataset6(double alphaValue, DatasetType nullModel, int firstMeasure, int lastMeasure) {

        String[] gValues = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma, minDst, maxDst, avg, area[] = new double[gValues.length];

        ExpParams par = null;
        String measure = "";

        int m,  kValues[] = {4, 6, 8, 10};

        ArrayList<ExpParams> lenList = getNullModelSeqLengths(nullModel, Paths.get( resultsPath, String.format("k=%d", kValues[0])));

        double[] tresh = new double[lenList.size()];
        System.out.printf("Distances for alternate model: %s\n", alternateModel.toString());

        ReadDistances ar = new ReadDistances(nullModel, alternateModel, lenList.get(0).filename, 0.); // solo per leggere la lista di misure disponibili
        lastMeasure = lastMeasure >= 0 ? lastMeasure : ar.measureNames.length;

        // per tutte le misure previste dall'esperimento
        for(int d = firstMeasure; d < lastMeasure; d++) {

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();
            m = 0;
            // per tutti i valori di k
            for(int k : kValues) {
                System.out.printf("Analyzing distance: %s for alpha = %.3f and k = %d\n", ar.measureNames[d], alphaValue, k);
                Path  kdir = Paths.get(resultsPath, String.format("k=%d", k));

                // per tutti i gamma presenti nelle properties (devono coincidere con quelli effettivamente utilizzati)
                for (int i = 0; i < gValues.length; i++) {

                    gamma = Double.parseDouble(gValues[i]);

                    // per tutte le lunghezze definite per il null model (cf. getNullModelSeqLengths())
                    for (int j = 0; j < lenList.size(); j++) {
                        par = lenList.get(j);

                        ar = new ReadDistances(nullModel, alternateModel, kdir.toString(), par.numPairs, k, par.seqLength, gamma);
                        // il costruttore legge automaticamente il nullModel
                        ar.readAlternateModel();
                        ar.readType1Check();

                        ar.setRefMeasure(d);        // ordina anche i risultati infunzione della misura scelta

                        ar.setThreshold(alphaValue);
                        tresh[j] = ar.getThreshold();
                        System.out.printf("Threshold(alpha = %.3f): %.3f\t", alphaValue, ar.getThreshold());

                        measure = ar.getCurrentMeasureName();
                        double power = ar.getPowerValue();

                        double t1 = ar.getType1Value();

                        System.out.printf("%s Power(k = %d, gamma = %.3f, n = %d) = %.4f, T1 = %.3f\n",
                                ar.getCurrentMeasureName(), k, gamma, par.seqLength, power, t1);

                        list.add(new PowerValue(k, par.seqLength, alphaValue, gamma, power, t1));
                    } // foreach len
                    double tot = 0;
                    minDst = Double.MAX_VALUE;
                    maxDst = Double.MIN_VALUE;
                    for (double t : tresh) {
                        tot += t;
                        minDst = Double.min(t, minDst);
                        maxDst = Double.max(t, maxDst);
                    }
                    avg = tot / lenList.size();
                    tot = 0;
                    for (double t : tresh) {
                        tot += Math.pow(t - avg, 2);
                    }
                    double var = tot / lenList.size();

                    area[i] = PowerValue.GetArea( list, m++, lenList.size());

                    System.out.printf("Mes: %s, alpha: %.3f, gamma: %.3f, area tot: %.3f, min: %.3f, max: %.3f, var: %.3f\n",
                            measure, alphaValue, gamma, area[i], minDst, maxDst, var);
                } // foreach i in gamma
            } // foreach k
            saveJSonResults(list, measure, par.k, nullModel.name(), alternateModel.name(), area);
            System.out.println("** Saved ***\n");
        } // foreach measure
    }


    // analisi del dataset large da 200.000 a 10.000.000 step 200.000
    static void AnalyzeDataset5(double alphaValue, DatasetType nullModel, int firstMeasure, int lastMeasure) {

        String[] gValues = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma, minDst, maxDst, avg, area[] = new double[gValues.length];

        ExpParams par = null;
        String measure = "";

        int m,  kValues[] = {4, 6, 8, 10};

        ArrayList<ExpParams> lenList = getNullModelSeqLengths(nullModel, Paths.get( resultsPath, String.format("k=%d", kValues[0])));

        double[] tresh = new double[lenList.size()];
        System.out.printf("Distances for alternate model: %s\n", alternateModel.toString());

        ReadDistances ar = new ReadDistances(nullModel, alternateModel, lenList.get(0).filename, 0.); // solo per leggere la lista di misure disponibili
        lastMeasure = lastMeasure >= 0 ? lastMeasure : ar.measureNames.length;

        // per tutte le misure previste dall'esperimento
        for(int d = firstMeasure; d < lastMeasure; d++) {

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();
            m = 0;
            // per tutti i valori di k
            for(int k : kValues) {
                System.out.printf("Analyzing distance: %s for alpha = %.3f and k = %d\n", ar.measureNames[d], alphaValue, k);
                Path  kdir = Paths.get(resultsPath, String.format("k=%d", k));

                // per tutti i gamma presenti nelle properties (devono coincidere con quelli effettivamente utilizzati)
                for (int i = 0; i < gValues.length; i++) {

                    gamma = Double.parseDouble(gValues[i]);

                    // per tutte le lunghezze definite per il null model (cf. getNullModelSeqLengths())
                    for (int j = 0; j < lenList.size(); j++) {
                        par = lenList.get(j);

                        ar = new ReadDistances(nullModel, alternateModel, kdir.toString(), par.numPairs, k, par.seqLength, gamma);
                        // il costruttore legge automaticamente il nullModel
                        ar.readAlternateModel();
                        ar.readType1Check();

                        ar.setRefMeasure(d);        // ordina anche i risultati infunzione della misura scelta

                        ar.setThreshold(alphaValue);
                        tresh[j] = ar.getThreshold();
                        System.out.printf("Threshold(alpha = %.3f): %.3f\t", alphaValue, ar.getThreshold());

                        measure = ar.getCurrentMeasureName();
                        double power = ar.getPowerValue();

                        double t1 = ar.getType1Value();

                        System.out.printf("%s Power(k = %d, gamma = %.3f, n = %d) = %.4f, T1 = %.3f\n",
                                ar.getCurrentMeasureName(), k, gamma, par.seqLength, power, t1);

                        list.add(new PowerValue(k, par.seqLength, alphaValue, gamma, power, t1));
                    } // foreach len
                    double tot = 0;
                    minDst = Double.MAX_VALUE;
                    maxDst = Double.MIN_VALUE;
                    for (double t : tresh) {
                        tot += t;
                        minDst = Double.min(t, minDst);
                        maxDst = Double.max(t, maxDst);
                    }
                    avg = tot / lenList.size();
                    tot = 0;
                    for (double t : tresh) {
                        tot += Math.pow(t - avg, 2);
                    }
                    double var = tot / lenList.size();

                    area[i] = PowerValue.GetArea( list, m++, lenList.size());

                    System.out.printf("Mes: %s, alpha: %.3f, gamma: %.3f, area tot: %.3f, min: %.3f, max: %.3f, var: %.3f\n",
                            measure, alphaValue, gamma, area[i], minDst, maxDst, var);
                } // foreach i in gamma
            } // foreach k
            saveJSonResults(list, measure, par.k, nullModel.name(), alternateModel.name(), area);
            System.out.println("** Saved ***\n");
        } // foreach measure
    }


    static void AnalyzeDataset4(double alphaValue,DatasetType nullModel) {

        String[] gValues = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma, minDst, maxDst, avg, area[] = new double[gValues.length];

        ExpParams par = null;
        String measure = "";

        createNewCSVFile( alternateModel, alphaValue);

        int m,  kValues[] = {4, 6, 8, 10};

        ArrayList<ExpParams> lenList = getNullModelSeqLengths(nullModel, Paths.get( resultsPath, Integer.toString(kValues[0])));

        double[] tresh = new double[lenList.size()];
        System.out.printf("Distances for alternate model: %s\n", alternateModel.toString());

        ReadDistances ar = new ReadDistances(nullModel, alternateModel, lenList.get(0).filename, 0.); // solo per leggere la lista di misure disponibili

        // per tutte le misure previste dall'esperimento
        for(int d = 0; d < ar.measureNames.length; d++) {

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();
            m = 0;
            // per tutti i valori di k
            for(int k : kValues) {
                System.out.printf("Analyzing distance: %s for alpha = %.3f and k = %d\n", ar.measureNames[d], alphaValue, k);
                Path  kdir = Paths.get(resultsPath, Integer.toString(k));

                // per tutti i gamma presenti nelle properties (devono coincidere con quelli effettivamente utilizzati)
                for (int i = 0; i < gValues.length; i++) {

                    gamma = Double.parseDouble(gValues[i]);

                    // per tutte le lunghezze definite per il null model (cf. getNullModelSeqLengths())
                    for (int j = 0; j < lenList.size(); j++) {
                        par = lenList.get(j);

                        ar = new ReadDistances(nullModel, alternateModel, kdir.toString(), par.numPairs, k, par.seqLength, gamma);
                        // il costruttore legge automaticamente il nullModel
                        ar.readAlternateModel();
                        ar.readType1Check();

                        ar.setRefMeasure(d);        // ordina anche i risultati infunzione della misura scelta

                        ar.setThreshold(alphaValue);
                        tresh[j] = ar.getThreshold();
                        System.out.printf("Threshold(%.3f): %.3f\t", alphaValue, ar.getThreshold());

                        measure = ar.getCurrentMeasureName();
                        double power = ar.getPowerValue();

                        double t1 = ar.getType1Value();

                        System.out.printf("%s Power(k = %d, gamma = %.3f, n = %d) = %.4f, T1 = %.3f\n",
                                ar.getCurrentMeasureName(), k, gamma, par.seqLength, power, t1);

                        list.add(new PowerValue(k, par.seqLength, alphaValue, gamma, power, t1));
                    } // foreach len
                    double tot = 0;
                    minDst = Double.MAX_VALUE;
                    maxDst = Double.MIN_VALUE;
                    for (double t : tresh) {
                        tot += t;
                        minDst = Double.min(t, minDst);
                        maxDst = Double.max(t, maxDst);
                    }
                    avg = tot / lenList.size();
                    tot = 0;
                    for (double t : tresh) {
                        tot += Math.pow(t - avg, 2);
                    }
                    double var = tot / lenList.size();

                    area[i] = PowerValue.GetArea( list, m++, lenList.size());

                    System.out.printf("Mes: %s, alpha: %.3f, gamma: %.3f, area tot: %.3f, min: %.3f, max: %.3f, var: %.3f\n",
                            measure, alphaValue, gamma, area[i], minDst, maxDst, var);
                } // foreach i in gamma
                saveCSVResults(measure, nullModel.name(), alternateModel.name(), k, area);
            } // foreach k
            saveJSonResults(list, measure, par.k, nullModel.name(), alternateModel.name(), area);
            System.out.println("** Saved ***\n");
        } // foreach measure
    }

    static void AnalyzeDataset3(double alphaValue, DatasetType nullModel) {

        String[] gValues = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma, minDst, maxDst, avg, area[] = new double[gValues.length];
        ExpParams par = null;
        String measure = "";

        createNewCSVFile( alternateModel, alphaValue);

        ArrayList<ExpParams> lenList = getNullModelSeqLengths(nullModel, Paths.get(resultsPath));

        double[] tresh = new double[lenList.size()];
        System.out.printf("Distances for alternate model: %s\n", alternateModel.toString());

        ReadDistances ar = new ReadDistances(nullModel, alternateModel, lenList.get(0).filename, 0.); // solo per leggere la lista di misure disponibili

        // per tutte le misure previste dall'esperimento
        for(int d = 0; d < ar.measureNames.length; d++) {

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();

            System.out.printf("Analyzing distance: %s for alpha = %.3f\n", ar.measureNames[d], alphaValue);

            // per tutti i gamma presenti nelle properties (devono coincidere con quelli effettivamente utilizzati)
            for (int i = 0; i < gValues.length; i++) {

                gamma = Double.parseDouble(gValues[i]);

                // per tutte le lunghezze definite per il null model (cf. getNullModelSeqLengths())
                for (int j = 0; j < lenList.size(); j++) {
                    par = lenList.get(j);

                    ar = new ReadDistances(nullModel, alternateModel, resultsPath, par.numPairs, par.k, par.seqLength, gamma);
                    // il costruttore legge automaticamente il nullModel
                    ar.readAlternateModel();
                    ar.readType1Check();

                    ar.setRefMeasure(d);        // ordina anche i risultati infunzione della misura scelta

                    ar.setThreshold(alphaValue);
                    tresh[j] = ar.getThreshold();
                    System.out.printf("Threshold(%.3f): %.3f\t", alphaValue, ar.getThreshold());

                    measure = ar.getCurrentMeasureName();
                    double power = ar.getPowerValue();

                    double t1 = ar.getType1Value();

                    System.out.printf("%s Power(k = %d, gamma = %.3f, n = %d) = %.4f, T1 = %.3f\n",
                            ar.getCurrentMeasureName(), par.k, gamma, par.seqLength, power, t1);

                    list.add(new PowerValue(par.k, par.seqLength, alphaValue, gamma, power, t1));
                } // foreach len
                double tot = 0;
                minDst = Double.MAX_VALUE;
                maxDst = Double.MIN_VALUE;
                for(double t : tresh) {
                    tot += t;
                    minDst = Double.min(t, minDst);
                    maxDst = Double.max(t, maxDst);
                }
                avg = tot / lenList.size();
                tot = 0;
                for(double t : tresh) {
                    tot += Math.pow( t - avg, 2);
                }
                double var = tot / lenList.size();

                area[i] = PowerValue.GetArea( list, i, lenList.size());

                System.out.printf("Mes: %s, alpha: %.3f, gamma: %.3f, area tot: %.3f, min: %.3f, max: %.3f, var: %.3f\n",
                        measure, alphaValue, gamma, area[i], minDst, maxDst, var);
            } // foreach gamma
            saveCSVResults(measure, nullModel.name(), alternateModel.name(), 0, area);
            saveJSonResults(list, measure, par.k, nullModel.name(), alternateModel.name(), area);
            System.out.println("** Saved ***\n");
        }
    }


    static void AnalyzeDataset2( double alfa) {

        String[] gVal = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma;
        ExpParams par = null;
        String measure = "";
        ArrayList<ExpParams> lenList = getNullModelSeqLengths(nullModel, Paths.get(resultsPath));

        ReadDistances ar = new ReadDistances(nullModel, alternateModel, lenList.get(0).filename, 0.); // solo per leggere la lista di misure disponibili

        // per tutte le misure previste dall'esperimento
        for(int d = 0; d < ar.measureNames.length; d++) {

            System.out.printf("Analyzing distance: %s for alpha = %.3f\n", ar.measureNames[d], alfa);

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();
            // per tutti i gamma presenti nelle properties (devono coincidere con quelli effettivamente utilizzati)
            for (int i = 0; i < gVal.length; i++) {

                gamma = Double.parseDouble(gVal[i]);

                // per tutte le lunghezze definite per il null model (cf. getNullModelSeqLengths())
                for (int j = 0; j < lenList.size(); j++) {
                    par = lenList.get(j);

                    ar = new ReadDistances(nullModel, alternateModel, resultsPath, par.numPairs, par.k, par.seqLength, gamma);
                    // il costruttore legge automaticamente il nullModel
                    ar.readAlternateModel();
                    ar.readType1Check();

                    ar.setRefMeasure(d);        // ordina anche i risultati infunzione della misura scelta

                    ar.setThreshold(alfa);
                    System.out.printf("Threshold(%.3f): %.3f\t", alfa, ar.getThreshold());

                    measure = ar.getCurrentMeasureName();
                    double power = ar.getPowerValue();

                    double t1 = ar.getType1Value();

                    System.out.printf("%s Power(k = %d, gamma = %.3f, n = %d) = %.4f, T1 = %.3f\n",
                            ar.getCurrentMeasureName(), par.k, gamma, par.seqLength, power, t1);

                    list.add(new PowerValue(par.k, par.seqLength, alfa, gamma, power, t1));
                }
            }
            saveJSonResults(list, measure, par.k, nullModel.name(), alternateModel.name(), null);
            System.out.println("** Saved ***");
        }
    }


    static void AnalyzeDataset1( int k, double alfa) {

        String[] stVal = appProps.getAppProps().getProperty("powerstatistics.datasetBuilder.gammaProbabilities").split(",");
        double gamma;
        int seqLen;

        ReadDistances ar = new ReadDistances(nullModel, alternateModel,
                resultsPath, 10000, k,200, 0.001);

        ar.readNullModel();
        for(int d = 0; d < ar.measureNames.length; d++) {

            ar.setRefMeasure( d);
            System.out.printf("Analyzing distance: %s for alpha = %.3f\n", ar.measureNames[d], alfa);

            ArrayList<PowerValue> list = new ArrayList<PowerValue>();
            String measure = "";

            for (int i = 0; i < stVal.length; i++) {

                gamma = Double.parseDouble(stVal[i]);

                for (int j = 1; j <= 8; j++) {
                    seqLen = (int) Math.pow(2, j) * 100;

                    ar = new ReadDistances(nullModel, alternateModel, resultsPath, 10000, k, seqLen, gamma);
                    // il costruttore legge automaticamente il nullModel
                    ar.readAlternateModel();

                    ar.setRefMeasure(d);   // ordina anche i risultati in funzione della misura scelta

                    ar.setThreshold(alfa);
                    System.out.printf("Threshold(%.3f): %.3f\t", alfa, ar.getThreshold());

                    measure = ar.getCurrentMeasureName();
                    double power = ar.getPowerValue();
                    System.out.printf("%s Power(k = %d, gamma = %.3f, n =%6d) = %.4f\n",
                            ar.getCurrentMeasureName(), k, gamma, seqLen, power);

                    list.add(new PowerValue(k, seqLen, alfa, gamma, power, 0.));
                }
            }
            saveJSonResults(list, measure, k, nullModel.name(), alternateModel.name(), null);
            System.out.println("** Saved ***");
        }
    }


    static ArrayList<ExpParams>  getNullModelSeqLengths(DatasetType nullModel, Path dir) {
        String pattern = "";

        pattern = "*" + nullModel.toString() + "-[0-9]*.csv";

        ArrayList<ExpParams> params = new ArrayList<ExpParams>();
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(dir, pattern)) {
            for (Path entry: stream) {
                params.add(ReadDistances.parseFilename(entry.toString()));
            }
        }
        catch (IOException ex) {
            System.err.printf("Error reading folder: %s %s\n", dir.toString(), ex.getMessage());
            return null;
        }
        if (params.size() == 0) {
            System.out.printf("getNullModelSeqLengths: No csv file found for pattern: %s\n", pattern);
            System.exit(-1);
        }
        params.sort(new ExpParams.SortBySequenceLegth());
        return params;
    }

    static void createNewCSVFile(DatasetType model, double alpha) {

        Current_CSV_File = String.format("%s-%s-%.3f.csv", AREAS_CSV_FILE, model.toString(), alpha);
        try (
                BufferedWriter writer = Files.newBufferedWriter(Paths.get(resultsPath, Current_CSV_File),
                        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
                CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT
                        .withHeader("Measure", "Null Model", "Alternate Model", "K", "PS G=0.010", "PS G=0.050", "PS G=0.100"));
        ) {
            csvPrinter.flush();
        }
        catch (Exception ex) {
            System.err.println("Cannot create areas: " + ex.getMessage());
        }
    }


    static void saveCSVResults(String measure, String nullModelName, String alternateModelName, int k, double[] area) {
        try (
                BufferedWriter writer = Files.newBufferedWriter(Paths.get(resultsPath, Current_CSV_File),
                                                                StandardOpenOption.CREATE, StandardOpenOption.APPEND);
                CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT)
        ) {
            DecimalFormat df = new DecimalFormat("#.####");
            csvPrinter.printRecord(measure, nullModelName, alternateModelName, k,
                    df.format(area[0]), df.format(area[1]), df.format(area[2]));
            csvPrinter.flush();
            csvPrinter.close();
        }
        catch (Exception ex) {
            System.err.println("Cannot save areas: " + ex.getMessage());
        }
    }


    static void saveJSonResults(ArrayList<PowerValue> list, String distance, int k,
                                String nullModel, String alternateModel, double[] area) {

        // JSONParser jsonParser = new JSONParser();

        try {
            // Object run = jsonParser.parse(new FileReader(String.format("%s.json", distance)));
            // JSONArray jsonArray = (JSONArray)run;
            JSONArray areas = new JSONArray();
            if (area != null) {
                for (double x : area) {
                    // JSONObject a = new JSONObject();
                    // a.put
                    areas.add( new Double(x));
                }
            }

            JSONObject header = new JSONObject();
            header.put("distanceName", distance);
            header.put("nullModel", nullModel);
            header.put("alternateModel", alternateModel);
            header.put("areas", areas);

            JSONArray powerValues = new JSONArray();
            JSONObject run = null;

            for (PowerValue x : list) {
                run = x.toJsonObject();
                powerValues.add(run);
            }

            JSONObject expReport = new JSONObject();
            expReport.put("header", header);
            expReport.put("values", powerValues);

            String al = String.format("%.3f", list.get(0).alpha);
            String jsonDir = String.format("%s/json-%s", resultsPath, al.substring(2));
            File dstDir = new File(jsonDir);
            if (!dstDir.exists()) {
                dstDir.mkdirs();
            }
            String fname = String.format("%s/%s-%s-%s.json", jsonDir, distance,
                    alternateModel.toString(), nullModel.toString().substring(0, 2));
            FileWriter file = new FileWriter(fname);
            file.write(expReport.toJSONString());
            file.flush();
            file.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
