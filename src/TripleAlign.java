import java.util.*;

public class TripleAlign {
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

    private static final int[][] M = {
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

    public Map<String, String> tripleAlign(Map<String, String> seqs, float g) {
        List<String> IDs = new ArrayList<>(seqs.keySet());
        String rID = IDs.get(0);
        String sID = IDs.get(1);
        String tID = IDs.get(2);
        String r = seqs.get(rID);
        String s = seqs.get(sID);
        String t = seqs.get(tID);
        int lenR = r.length();
        int lenS = s.length();
        int lenT = t.length();
        float[][][] S = new float[lenR+1][lenS+1][lenT+1];
        String[][][] B = new String[lenR+1][lenS+1][lenT+1];

        for(int x=1; x<lenR+1; ++x){
            S[x][0][0] = 2*x*g;
            B[x][0][0] = "x";
        }

        for(int y=1; y<lenS+1; ++y){
            S[0][y][0] = 2*y*g;
            B[0][y][0] = "y";
        }

        for(int z=1; z<lenT+1; ++z){
            S[0][0][z] = 2*z*g;
            B[0][0][z] = "z";
        }

        for(int x=1; x<lenR+1; ++x){
            for(int y=1; y<lenS+1; ++y){
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("xy", S[x-1][y-1][0] + M[MIndirect.get(r.charAt(x-1))][MIndirect.get(s.charAt(y-1))] + 2*g);
                options.put("x", S[x-1][y][0] + 2*g);
                options.put("y", S[x][y-1][0] + 2*g);
                S[x][y][0] = Collections.max(options.values());
                for (String o : options.keySet()) {
                    if (S[x][y][0] == options.get(o)) {
                        B[x][y][0] = o;
                        break;
                    }
                }
            }
        }

        for(int x=1; x<lenR+1; ++x){
            for(int z=1; z<lenT+1; ++z){
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("xz", S[x-1][0][z-1] + M[MIndirect.get(r.charAt(x-1))][MIndirect.get(t.charAt(z-1))] + 2*g);
                options.put("x", S[x-1][0][z] + 2*g);
                options.put("z", S[x][0][z-1] + 2*g);
                S[x][0][z] = Collections.max(options.values());
                for (String o : options.keySet()) {
                    if (S[x][0][z] == options.get(o)) {
                        B[x][0][z] = o;
                        break;
                    }
                }
            }
        }

        for(int y=1; y<lenS+1; ++y){
            for(int z=1; z<lenT+1; ++z){
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("yz", S[0][y-1][z-1] + M[MIndirect.get(s.charAt(y-1))][MIndirect.get(t.charAt(z-1))] + 2*g);
                options.put("y", S[0][y-1][z] + 2*g);
                options.put("z", S[0][y][z-1] + 2*g);
                S[0][y][z] = Collections.max(options.values());
                for (String o : options.keySet()) {
                    if (S[0][y][z] == options.get(o)) {
                        B[0][y][z] = o;
                        break;
                    }
                }
            }
        }

        for(int x=1; x<lenR+1; ++x){
            for(int y=1; y<lenS+1; ++y) {
                for (int z = 1; z < lenT + 1; ++z) {
                    Map<String, Float> options = new HashMap<String, Float>();
                    options.put("xyz", S[x - 1][y - 1][z - 1] + M[MIndirect.get(r.charAt(x - 1))][MIndirect.get(s.charAt(y - 1))]
                            + M[MIndirect.get(r.charAt(x - 1))][MIndirect.get(t.charAt(z - 1))]
                            + M[MIndirect.get(s.charAt(y - 1))][MIndirect.get(t.charAt(z - 1))]);
                    options.put("xy", S[x - 1][y - 1][z] + M[MIndirect.get(r.charAt(x - 1))][MIndirect.get(s.charAt(y - 1))] + 2 * g);
                    options.put("xz", S[x - 1][y][z - 1] + M[MIndirect.get(r.charAt(x - 1))][MIndirect.get(t.charAt(z - 1))] + 2 * g);
                    options.put("yz", S[x][y - 1][z - 1] + M[MIndirect.get(s.charAt(y - 1))][MIndirect.get(t.charAt(z - 1))] + 2 * g);
                    options.put("x", S[x - 1][y][z] + 2 * g);
                    options.put("y", S[x][y - 1][z] + 2 * g);
                    options.put("z", S[x][y][z - 1] + 2 * g);
                    S[x][y][z] = Collections.max(options.values());
                    for (String o : options.keySet()) {
                        if (S[x][y][z] == options.get(o)) {
                            B[x][y][z] = o;
                            break;
                        }
                    }
                }
            }
        }


        int x = lenR;
        int y = lenS;
        int z = lenT;
        Map<String, String> out = new HashMap<String, String>();
        out.put(rID, "");
        out.put(sID, "");
        out.put(tID, "");



        while(x+y+z != 0) {
            String bt = B[x][y][z];
            if(bt.contains("x")){
                String v = out.get(rID);
                out.put(rID, v + r.charAt(x-1));
                x -= 1;
            } else {
                String v = out.get(rID);
                out.put(rID, v + "-");
            }
            if(bt.contains("y")){
                String v = out.get(sID);
                out.put(tID, v + s.charAt(y-1));
                y -= 1;
            } else {
                String v = out.get(sID);
                out.put(sID, v + "-");
            }
            if(bt.contains("z")){
                String v = out.get(tID);
                out.put(tID, v + t.charAt(z-1));
                z -= 1;
            } else {
                String v = out.get(tID);
                out.put(tID, v + "-");
            }
        }

        for (String ID : out.keySet()) {
            String rev = new StringBuilder(out.get(ID)).reverse().toString();
            out.put(ID, rev);
        }

        return out;
    }
}
