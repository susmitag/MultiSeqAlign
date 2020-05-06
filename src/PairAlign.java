import java.util.*;

public class PairAlign {

    //BLOSUM62
    //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V

    Map<String, Integer> MIndirect = new HashMap<String, Integer>() {{
        put("A", 0);
        put("R", 1);
        put("N", 2);
        put("D", 3);
        put("C", 4);
        put("Q", 5);
        put("E", 6);
        put("G", 7);
        put("H", 8);
        put("I", 9);
        put("L", 10);
        put("K", 11);
        put("M", 12);
        put("F", 13);
        put("P", 14);
        put("S", 15);
        put("T", 16);
        put("W", 17);
        put("Y", 18);
        put("V", 19);
    }};

    int[][] M = {
            { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
            {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
            {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
            {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
            { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
            {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
            {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
            { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
            {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
            {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
            {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
            {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
            {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
            {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
            {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
            { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
            { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
            {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
            {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
            { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};

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

        // perform DP 1D
        for(int i=1; i<lenS+1; ++i){
            S[i][0] = i*gap;
            B[i][0] = "i";
        }

        for(int j=1; j<lenT+1; ++j){
            S[0][j] = j*gap;
            B[0][j] = "j";
        }

        // perform DP 2D
        for(int i=1; i<lenS+1; ++i){
            for(int j=1; j<lenT+1; ++j){
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("ij", S[i-1][j-1] + M[MIndirect.get(Character.toString(s.charAt(i-1)))][MIndirect.get(Character.toString(t.charAt(j-1)))]);
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

        // backtrack
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
                out.put(tID, v + t.charAt(j-1));
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
