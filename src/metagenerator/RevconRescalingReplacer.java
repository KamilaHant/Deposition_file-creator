package metagenerator;

import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

public class RevconRescalingReplacer {
    static void revconRescaling() {
        File textFile = new File(Variables.newPath + "\\REVCON_to_CONFIG_rescaling.c");

        try {
            String data = FileUtils.readFileToString(textFile);
            data = data.replace("$T", Variables.T);
            FileUtils.writeStringToFile(textFile, data);
        } catch (IOException var2) {
            var2.printStackTrace();
        }

    }
}
