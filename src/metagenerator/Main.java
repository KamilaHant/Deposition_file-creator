package metagenerator;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class Main {

    static String dataDir = "D:\\generovani_input\\data\\";
    static String configDir = "D:\\generovani_input\\";



    public static void main(String[] args) throws IOException {
        try {
            ConfigReader.reader();
            String key = Variables.mutliKey();
            if (key != null) {
                for(int i = 0; i < ((ArrayList)ConfigReader.config.get(key)).toArray().length; ++i) {
                    Variables.getter(i);
                    Variables.pathName();
                    FileMaker.newfilemaker();
                    JobfileReplacer.replace();
                    DepositionReplacer.depositionReplacer();
                    ControlReplacer.controlNVE();
                    ControlReplacer.controlNVT();
                    RevconRescalingReplacer.revconRescaling();
                    ZpracovaniReplacer.zpracovaniReplacer();
                }
            }
        } catch (FileNotFoundException var3) {
            var3.printStackTrace();
        }

    }

}
