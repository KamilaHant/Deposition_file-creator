package metagenerator;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Scanner;

public class ConfigReader {
    static Hashtable<String, ArrayList<String>> config;
    static Hashtable<String, String> parametres;


    static void reader() throws FileNotFoundException {
        config = new Hashtable();
        parametres = new Hashtable();
        Scanner readData = new Scanner(new File(Main.configDir + "config"));

        while(readData.hasNext()) {
            String auxiliary = readData.nextLine();
            ArrayList<String> line = new ArrayList(Arrays.asList(auxiliary.split(" ")));
            String hash = (String)line.get(0);
            line.remove(0);
            config.put(hash, line);
        }

        readData = new Scanner(new File(Main.configDir+ "parametres"));

        while(readData.hasNext()) {
            String[] auxiliary = readData.nextLine().split(" ");
            parametres.put(auxiliary[0], auxiliary[1]);
        }

        System.out.print("---reading of config and parametres was succesfull---");
    }
}
