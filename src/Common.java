import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Common {

    public Map<String, String> parseFasta() throws IOException {
        Map<String, String> seqs = new HashMap<String, String>();
        FileInputStream fstream = new FileInputStream("/Users/goswama2/Documents/personal/123/out/production/MultiSeqAlign/0.fas");
        BufferedReader br = new BufferedReader(new InputStreamReader(fstream));

        String strLine;

        String currID = null;
        String currSeq = "";

        while ((strLine = br.readLine()) != null)   {
            String line = strLine.trim();
            System.out.println (line);
            if (line.length() == 0)
                continue;
            if (line.charAt(0) == '>') {
                if (currID != null) {
                    seqs.put(currID, currSeq);
                    currSeq = "";
                }
                currID = line.substring(1);
            }
            else {
                currSeq += line;
            }
        }
        seqs.put(currID, currSeq);

        Iterator it = seqs.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            System.out.println(pair.getKey() + " = " + pair.getValue());
        }

        fstream.close();
        return seqs;
    }
}
