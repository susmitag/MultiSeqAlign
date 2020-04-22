import java.util.*;

public class PairAlign {

    Map<String, Map <String, Float>> M;

    {
        M = new HashMap<String, Map <String, Float>>();
        Map <String, Float> mapA = new HashMap<String, Float>();
        mapA.put("A", (float)4);
        mapA.put("C", (float)0);
        M.put("A", mapA);
    }

    float pairAlign(Map<String, String> seqs, float gap) {
        List<String> IDs = new ArrayList<>(seqs.keySet());
        Collections.sort(IDs);
        String sID = IDs.get(0);
        String tID = IDs.get(1);
        String s = seqs.get(sID);
        String t = seqs.get(tID);
        int lenS = s.length();
        int lenT = t.length();
        float[][] S = new float[lenS+1][lenT+1];
        String[][] B = new String[lenS+1][lenT+1];

        for(int i=1; i<lenS+1; ++i){
            S[i][0] = i*gap;
            B[i][0] = "i";
        }

        for(int j=1; j<lenT+1; ++j){
            S[0][j] = j*gap;
            B[0][j] = "j";
        }

        for(int i=1; i<lenS+1; ++i){
            for(int j=1; j<lenT+1; ++j){
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("ij", S[i-1][j-1] + M.get(s.charAt(i-1)).get(t.charAt(j-1));
                options.put("i", S[i-1][j] + gap);
                options.put("j", S[i][j-1] + gap);
                S[i][j] = Collections.max(options.values());
                for (String o : options.keySet()) {
                    if (S[i][j] == options.get(o)) {
                        B[i][j] = o;
                        break;
                    }
                }
            }
        }

        int i = lenS;
        int j = lenT;
        Map<String, String> out = new HashMap<String, String>();
        out.put(sID, "");
        out.put(tID, "");

        while(i+j != 0) {
            String bt = B[i][j];
            if(bt.contains("i")){
                String v = out.get(sID);
                out.put(sID, v + s.charAt(i-1));
                i -= 1;
            } else {
                String v = out.get(sID);
                out.put(sID, v + "-");
            }
            if(bt.contains("j")){
                String v = out.get(tID);
                out.put(tID, v + s.charAt(j-1));
                j -= 1;
            } else {
                String v = out.get(tID);
                out.put(tID, v + "-");
            }
        }

        String rev = new StringBuilder(out.get(sID)).reverse().toString();
        out.put(sID, rev);
        rev = new StringBuilder(out.get(tID)).reverse().toString();
        out.put(tID, rev);

        return S[lenS][lenT];
    }
}
