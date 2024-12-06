package it.unisa.di.bio.unitTest;

import java.io.FileInputStream;
import java.util.Hashtable;
import java.util.Scanner;

public class KMersStatistics {

    Hashtable<String, Integer> histogram;

    public KMersStatistics( String fname) {

        histogram = new Hashtable<String, Integer>(100000);
        int cnt = 0;
        try {
            FileInputStream fis = new FileInputStream(fname);
            Scanner sc = new Scanner(fis);

            while (sc.hasNextLine()) {
                String[] buf = sc.nextLine().split("\t");
                if (buf.length != 2) {
                    throw new Exception("Malformed input line");
                }
                histogram.put(buf[0], new Integer(buf[1]));
                cnt++;
            }
            System.out.printf("%d k-mers counters read\n", cnt);
        }
        catch (Exception ex) {
            System.err.println("Error reading input file: " + fname + " " + ex.getMessage());
            ex.printStackTrace();
        }
    }
}
