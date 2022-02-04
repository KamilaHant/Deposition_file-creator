package metagenerator;

import java.io.IOException;
import java.util.ArrayList;

public class Variables {
    static String Ehigh;
    static String Elow;
    static String R;
    static String x;
    static String T;
    static String newPath;
    static String key;


    static void getter(int step) {
        key = mutliKey();
        if (key == "Ehigh") {
            Ehigh = (String)((ArrayList)ConfigReader.config.get("Ehigh")).get(step);
            Elow = (String)((ArrayList)ConfigReader.config.get("Elow")).get(0);
            R = (String)((ArrayList)ConfigReader.config.get("R")).get(0);
            x = (String)((ArrayList)ConfigReader.config.get("x")).get(0);
            T = (String)((ArrayList)ConfigReader.config.get("T")).get(0);
        }

        if (key == "Elow") {
            Ehigh = (String)((ArrayList)ConfigReader.config.get("Ehigh")).get(0);
            Elow = (String)((ArrayList)ConfigReader.config.get("Elow")).get(step);
            R = (String)((ArrayList)ConfigReader.config.get("R")).get(0);
            x = (String)((ArrayList)ConfigReader.config.get("x")).get(0);
            T = (String)((ArrayList)ConfigReader.config.get("T")).get(0);
        }

        if (key == "R") {
            Ehigh = (String)((ArrayList)ConfigReader.config.get("Ehigh")).get(0);
            Elow = (String)((ArrayList)ConfigReader.config.get("Elow")).get(0);
            R = (String)((ArrayList)ConfigReader.config.get("R")).get(step);
            x = (String)((ArrayList)ConfigReader.config.get("x")).get(0);
            T = (String)((ArrayList)ConfigReader.config.get("T")).get(0);
        }

        if (key == "x") {
            Ehigh = (String)((ArrayList)ConfigReader.config.get("Ehigh")).get(0);
            Elow = (String)((ArrayList)ConfigReader.config.get("Elow")).get(0);
            R = (String)((ArrayList)ConfigReader.config.get("R")).get(0);
            x = (String)((ArrayList)ConfigReader.config.get("x")).get(step);
            T = (String)((ArrayList)ConfigReader.config.get("T")).get(0);
        }

        if (key == "T") {
            Ehigh = (String)((ArrayList)ConfigReader.config.get("Ehigh")).get(0);
            Elow = (String)((ArrayList)ConfigReader.config.get("Elow")).get(0);
            R = (String)((ArrayList)ConfigReader.config.get("R")).get(0);
            x = (String)((ArrayList)ConfigReader.config.get("x")).get(0);
            T = (String)((ArrayList)ConfigReader.config.get("T")).get(step);
        }

    }

    static void pathName() throws IOException {
        newPath = "D:\\generovani_input-new\\\\ZnO" + x + "_" + Ehigh + "-" + Elow + "_" + R + "_" + T + "K";
    }

    static String mutliKey() {
        if (((ArrayList)ConfigReader.config.get("Ehigh")).toArray().length > 1) {
            key = "Ehigh";
        }

        if (((ArrayList)ConfigReader.config.get("Elow")).toArray().length > 1) {
            key = "Elow";
        }

        if (((ArrayList)ConfigReader.config.get("R")).toArray().length > 1) {
            key = "R";
        }

        if (((ArrayList)ConfigReader.config.get("x")).toArray().length > 1) {
            key = "x";
        }

        if (((ArrayList)ConfigReader.config.get("T")).toArray().length > 1) {
            key = "T";
        }

        return key;
    }
}
