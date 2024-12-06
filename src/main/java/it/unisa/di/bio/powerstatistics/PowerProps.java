package it.unisa.di.bio.powerstatistics;

import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Properties;
import java.util.ResourceBundle;


public class PowerProps {

    // private String fname= "PowerStatistics.properties";
    private String fname= "PowerStatistics.properties";

    private Properties appProps;

    private static PowerProps instance;

    private PowerProps(){
        try {

            appProps = new Properties();
            ClassLoader classloader = Thread.currentThread().getContextClassLoader();
            InputStream is = classloader.getResourceAsStream(fname);
            if (is != null)
                appProps.load(is);
            else
                throw new FileNotFoundException("resource file: " + fname);

            // appProps.list(System.out);
            is.close();
        } catch (Exception ex) {
            System.err.printf("%s Resource file not found: %s. Exiting.",
                    this.getClass().getCanonicalName(), ex.getMessage());
            System.exit(-1);
        }
    }

    public static PowerProps getInstance(){
        if(instance == null){
            instance = new PowerProps();
        }
        return instance;
    }

    public Properties getAppProps() {
        return appProps;
    }

}
