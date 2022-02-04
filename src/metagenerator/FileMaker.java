package metagenerator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

public class FileMaker {

    static void newfilemaker() throws IOException {
        String sourceDirectoryLocation = Main.dataDir;
        String destinationDirectoryLocation = Variables.newPath;
        Files.walk(Paths.get(sourceDirectoryLocation)).forEach((source) -> {
            Path destination = Paths.get(destinationDirectoryLocation, source.toString().substring(sourceDirectoryLocation.length() - 1));

            try {
                Files.copy(source, destination, StandardCopyOption.REPLACE_EXISTING);
            } catch (IOException var5) {
                var5.printStackTrace();
            }

        });
    }
}
