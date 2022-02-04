package metagenerator;

import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;

public class JobfileReplacer {

    static void replace() throws IOException {
        File textFile = new File(Variables.newPath + "\\Jobfile");
        String replaceLine = "";

        try {
            String line = FileUtils.readFileToString(textFile);
            line = line.replace("$name", "ZnO" + Variables.x + "_" + Variables.Ehigh + "-" + Variables.Elow + "_" + Variables.R);
            line = line.replace("$walltime", (CharSequence)ConfigReader.parametres.get("walltime"));
            line = line.replace("$datadir", Variables.newPath);
            line = line.replace("$zpracovani", (CharSequence)ConfigReader.parametres.get("zpracovani"));
            if (line.contains("$lines")) {
                int steps = Integer.parseInt((String)ConfigReader.parametres.get("steps"));

                for(int i = 0; i <= steps; ++i) {
                    if (i >= 10 && i <= steps - 10) {
                        replaceLine = replaceLine + linesMakerShort(i) + "\n";
                    } else {
                        replaceLine = replaceLine + linesMakerLong(i) + "\n\n";
                    }
                }
            }

            line = line.replace("$lines", replaceLine);
            FileUtils.writeStringToFile(textFile, line);
        } catch (IOException var5) {
            var5.printStackTrace();
        }

    }

    static String linesMakerLong(int step) {
        String lines = "./Deposition_ZnO_LAMMPS.exe \n./CONFIG_to_SNAPSHOT.exe \n";
        lines = lines + "mv SNAPSHOT files/SNAPSHOT_" + step + "_I\n";
        lines = lines + "cp CONFIG CONFIG_" + step + "_I\n";
        lines = lines + "mv CONFIG_" + step + "_I files/\n";
        lines = lines + "cp CONTROL_nve CONTROL\n";
        lines = lines + "mv REVCON files/REVCON_" + (step - 1) + "_III\n";
        lines = lines + "mpiexec -np 1 lammps < CONTROL > OUTPUT\n";
        lines = lines + "mv OUTPUT files/OUTPUT_" + step + "_deposition\n";
        lines = lines + "./REVCON_to_CONFIG_rescaling.exe\n./CONFIG_to_SNAPSHOT.exe\n";
        lines = lines + "mv SNAPSHOT files/SNAPSHOT_" + step + "_II\n";
        lines = lines + "cp CONFIG CONFIG_" + step + "_II\n";
        lines = lines + "mv CONFIG_" + step + "_II files/\n";
        lines = lines + "cp CONTROL_nvt CONTROL\n";
        lines = lines + "mv REVCON files/REVCON_" + step + "_II\n";
        lines = lines + "mpiexec -np 1 lammps < CONTROL > OUTPUT\n";
        lines = lines + "mv OUTPUT files/OUTPUT_" + step + "_thermalization\n";
        lines = lines + "./REVCON_to_CONFIG_rescaling.exe\n./CONFIG_to_SNAPSHOT.exe\n";
        lines = lines + "mv SNAPSHOT files/SNAPSHOT_" + step + "_III\n";
        lines = lines + "cp CONFIG CONFIG_" + step + "_III\n";
        lines = lines + "mv CONFIG_" + step + "_III files/";
        return lines;
    }

    static String linesMakerShort(int step) {
        String lines = "./Deposition_ZnO_LAMMPS.exe \ncp CONTROL_nve CONTROL\n";
        lines = lines + "mpiexec -np 1 lammps < CONTROL > OUTPUT\n";
        lines = lines + "./REVCON_to_CONFIG_rescaling.exe\n";
        lines = lines + "cp CONTROL_nvt CONTROL\n";
        lines = lines + "mpiexec -np 1 lammps < CONTROL > OUTPUT\n";
        lines = lines + "./REVCON_to_CONFIG_rescaling.exe\n./CONFIG_to_SNAPSHOT.exe\n";
        lines = lines + "mv SNAPSHOT files/SNAPSHOT_" + step + "_III\n";
        return lines;
    }
}
