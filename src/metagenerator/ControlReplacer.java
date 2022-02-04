package metagenerator;


import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

public class ControlReplacer {

    static void controlNVT() {
        File textFile = new File(Variables.newPath, "\\CONTROL_nvt");

        try {
            String data = FileUtils.readFileToString(textFile);
            data = data.replace("$frozen", "frozen");
            data = data.replace("$name1", "name1");
            data = data.replace("$name2", "name2");
            data = data.replace("$T", Variables.T);
            data = data.replace("$damp", "damp");
            data = data.replace("$runNVT", "runNVT");
            data = data.replace("$thermo", "themo");
            data = data.replace("$timestep", "timestep");
            FileUtils.writeStringToFile(textFile, data);
        } catch (IOException var2) {
            var2.printStackTrace();
        }

    }

    static void controlNVE() {
        File textFile = new File(Variables.newPath, "\\CONTROL_nve");

        try {
            String data = FileUtils.readFileToString(textFile);
            data = data.replace("$frozen", "frozen");
            data = data.replace("$atoms", "atoms");
            data = data.replace("$name1", "name1");
            data = data.replace("$name2", "name2");
            data = data.replace("$T", Variables.T);
            data = data.replace("$damp", "damp");
            data = data.replace("$runNVE", "runNVE");
            data = data.replace("$thermo", "themo");
            data = data.replace("$timestep", "timestep");
            FileUtils.writeStringToFile(textFile, data);
        } catch (IOException var2) {
            var2.printStackTrace();
        }

    }
}