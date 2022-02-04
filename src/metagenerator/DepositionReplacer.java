package metagenerator;

import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

public class DepositionReplacer {

    static void depositionReplacer() {
        File textFile = new File(Variables.newPath, "\\Deposition_ZnO_LAMMPS.c");

        try {
            String data = FileUtils.readFileToString(textFile);
            data = data.replace("$frozen1", (CharSequence)ConfigReader.parametres.get("frozen1"));
            data = data.replace("$frozen2", (CharSequence)ConfigReader.parametres.get("frozen2"));
            data = data.replace("$unfrozen1", (CharSequence)ConfigReader.parametres.get("unfrozen1"));
            data = data.replace("$unfrozen1", (CharSequence)ConfigReader.parametres.get("unfrozen2"));
            data = data.replace("$defineAtoms", (CharSequence)ConfigReader.parametres.get("defineAtoms"));
            data = data.replace("$Ehigh", Variables.Ehigh);
            data = data.replace("$Elow", Variables.Elow);
            data = data.replace("$R", Variables.R);
            data = data.replace("$x", Variables.x);
            data = data.replace("$mass1", (CharSequence)ConfigReader.parametres.get("mass1"));
            data = data.replace("$mass2", (CharSequence)ConfigReader.parametres.get("mass2"));
            data = data.replace("$charge1", (CharSequence)ConfigReader.parametres.get("charge1"));
            data = data.replace("$charge2", (CharSequence)ConfigReader.parametres.get("charge2"));
            data = data.replace("$name1", (CharSequence)ConfigReader.parametres.get("name1"));
            data = data.replace("$name2", (CharSequence)ConfigReader.parametres.get("name2"));
            data = data.replace("$bond12", (CharSequence)ConfigReader.parametres.get("bond12"));
            data = data.replace("$cutoff11", (CharSequence)ConfigReader.parametres.get("cutoff11"));
            data = data.replace("$cutoff12", (CharSequence)ConfigReader.parametres.get("cutoff12"));
            data = data.replace("$cutoff22", (CharSequence)ConfigReader.parametres.get("cutoff22"));
            data = data.replace("$thickness", (CharSequence)ConfigReader.parametres.get("thickness"));
            data = data.replace("$cellx", (CharSequence)ConfigReader.parametres.get("cellx"));
            data = data.replace("$celly", (CharSequence)ConfigReader.parametres.get("celly"));
            FileUtils.writeStringToFile(textFile, data);
        } catch (IOException var2) {
            var2.printStackTrace();
        }

    }
}
