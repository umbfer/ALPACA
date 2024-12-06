package it.unisa.di.bio.powerstatistics;

import javax.json.JsonException;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.util.ArrayList;

import static java.lang.Math.log10;


public class PowerValue {

    double  t1;
    double  power;
    double  gamma;
    double  alpha;
    int     k;
    int     len;

    public PowerValue( int k, int len, double alpha, double gamma, double power, double t1) {
        this.t1     = t1;
        this.power  = power;
        this.gamma  = gamma;
        this.alpha  = alpha;
        this.k      = k;
        this.len    = len;
    }

    public JSONObject toJsonObject() {
        JSONObject obj = new JSONObject();
        try {
            obj.put("k", this.k);
            obj.put("len", this.len);
            obj.put("alpha", this.alpha);
            obj.put("gamma", this.gamma);
            obj.put("power", this.power);
            obj.put("t1", this.t1);
        } catch (
                JsonException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return obj;
    }

    public static double GetArea(ArrayList<PowerValue> list, int start, int len) {

        double maxlen = list.get((start+1)*len - 1).len; // l'ultimo dell'elenco e' la max length
        double areaTotale = 0.;

        for(int i = start * len; i < (start + 1) * len; i++) {
//            PowerValue v1 = list.get(i-1);
//            PowerValue v2 = list.get(i);
//            double h1 = v1.power;
//            double h2 = v2.power;
//            double x1 = log10((v1.len * 10) / maxlen); // normalizzato a 10 => log10 = 1
//            double x2 = log10((v2.len * 10) / maxlen);
//            areaTotale += (x2 - x1) * h1 + (((h2 - h1) * (x2 - x1)) / 2); // somma o sottrae l'area del triangolo a secondo h2 > h1 o no
            PowerValue v1 = list.get(i);
            double h1 = v1.power / 17.;     // max 1 anzichè 17
            areaTotale += h1;
        }
        return areaTotale;
    }
}
