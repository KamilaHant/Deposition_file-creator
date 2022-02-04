package metagenerator;

import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

public class ZpracovaniReplacer {
    static void zpracovaniReplacer() {
        File textFile = new File(Variables.newPath + "\\" + (String)ConfigReader.parametres.get("zpracovani"));

        try {
            String data = FileUtils.readFileToString(textFile);
            data = data.replace("$frozen1", (CharSequence)ConfigReader.parametres.get("frozen1"));
            data = data.replace("$frozen2", (CharSequence)ConfigReader.parametres.get("frozen2"));
            data = data.replace("$unfrozen1", (CharSequence)ConfigReader.parametres.get("unfrozen1"));
            data = data.replace("$unfrozen1", (CharSequence)ConfigReader.parametres.get("unfrozen2"));
            data = data.replace("$name1", (CharSequence)ConfigReader.parametres.get("name1"));
            data = data.replace("$name2", (CharSequence)ConfigReader.parametres.get("name2"));
            data = data.replace("$cutoff12", (CharSequence)ConfigReader.parametres.get("cutoff12"));
            data = data.replace("$cellx", (CharSequence)ConfigReader.parametres.get("cellx"));
            data = data.replace("$celly", (CharSequence)ConfigReader.parametres.get("celly"));
            FileUtils.writeStringToFile(textFile, data);
        } catch (IOException var2) {
            var2.printStackTrace();
        }

    }
}
