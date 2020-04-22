import java.util.*;

public class MultiAlign {

    Map<String, Map <String, Float>> M;

    {
        M = new HashMap<String, Map <String, Float>>();
        Map <String, Float> mapA = new HashMap<String, Float>();
        mapA.put("A", (float)4);
        mapA.put("C", (float)0);
        M.put("A", mapA);
    }

    public float mergeAlign(Map<String, String> aln, String ID, String seq, float gap) {
        List<String> keys = new ArrayList<>(aln.keySet());
        Collections.sort(keys);

        int n = aln.get(keys.get(0)).length();
        int m = seq.length();
        float[][] S = new float[n+1][m+1];
        String[][] B = new String[n+1][m+1];

        float[] sop = new float[n];
        for (int iKey = 0; iKey < keys.size() - 1; ++iKey) {
            String s = aln.get(keys.get(iKey));
            for (int jKey = iKey + 1; jKey < keys.size(); ++jKey) {
                String t = aln.get(keys.get(iKey));
                for (int i = 0; i < n; ++i) {
                    if (((s.charAt(i) == '-') && (t.charAt(i) != '-')) || ((s.charAt(i) != '-') && (t.charAt(i) == '-'))) {
                        sop[i] += gap;
                    } else if ((s.charAt(i) != '-') && (t.charAt(i) != '-')) {
                        sop[i] += M.get(s.charAt(i)).get(t.charAt(i));
                    }
                }
            }
        }

        for (int i=1; i<n+1; ++i) {
            S[i][0] = S[i-1][0] + keys.size()*gap;
            B[i][0] = "i";
        }
        for (int j=1; j<m+1; ++j) {
            S[0][j] = S[0][j-1] + keys.size()*gap;
            B[0][j] = "j";
        }
        for (int i=1; i<n+1; ++i) {
            for (int j=1; j<m+1; ++j) {
                Map<String, Float> options = new HashMap<String, Float>();
                options.put("i", S[i-1][j] + sop[i-1] + keys.size()*gap);
                options.put("j", S[i][j-1] + keys.size()*gap);
                options.put("ij", S[i-1][j-1] + sop[i-1]);
                for (String s : aln.values()) {
                    Float v = options.get("ij");
                    if (s.charAt(i-1) == '-') {
                        options.put("ij", v+gap);
                    } else {
                        options.put("ij", v+M.get(s.charAt(i-1)).get(seq.charAt(j-1)));
                    }
                }
                S[i][j] = Collections.max(options.values());
                for (String o : options.keySet()) {
                    if (S[i][j] == options.get(o)) {
                        B[i][j] = o;
                        break;
                    }
                }
            }
        }
        int i = n;
        int j = m;
        Map<String, String> out = new HashMap<String, String>();
        for (String key : keys) {
            out.put(key, "");
        }
        out.put(ID, "");
        float outScore = S[i][j];
        while (i+j != 0) {
            String bt = B[i][j];
            if(bt.contains("i")){
                for (String key : keys) {
                    String v = out.get(key);
                    out.put(key, v + aln.get(key).charAt(i-1));
                }
                i -= 1;
            } else {
                for (String key : keys) {
                    String v = out.get(key);
                    out.put(key, v + "-");
                }
            }
            if(bt.contains("j")){
                String v = out.get(ID);
                out.put(ID, v + seq.charAt(j-1));
                j -= 1;
            } else {
                String v = out.get(ID);
                out.put(ID, v + "-");
            }
        }
        for(String key : out.keySet()) {
            String rev = new StringBuilder(out.get(key)).reverse().toString();
            aln.put(key, rev);
        }
        return outScore;
    }

    public class Entry implements Comparable<Entry> {
        private int key;
        private String value;

        public Entry(int key, String value) {
            this.key = key;
            this.value = value;
        }

        // getters
        public int getKey() {
            return key;
        }

        public String getValue() {
            return value;
        }

        @Override
        public int compareTo(Entry other) {
            if (this.getKey() < other.getKey())
                return 1;
            else if (this.getKey() > other.getKey())
                return -1;
            return 0;
        }
    }

    public Map<String, String> multiAlign(Map<String, String> seqs, float gap) {
        List<String> IDs = new ArrayList<>(seqs.keySet());
        Collections.sort(IDs);
        Map<String, Map<String, Integer>> dm = new HashMap<String, Map<String, Integer>>();
        for(String ID : IDs) {
            dm.put(ID, new HashMap<String, Integer>());
        }

        PairAlign pA = new PairAlign();
        for (int i = 0; i < IDs.size() - 1; ++i) {
            for (int j = i + 1; j < IDs.size(); ++i) {
                String iID = IDs.get(i);
                String jID = IDs.get(j);
                Map<String, String> in = new HashMap<String, String>();
                in.put(iID, seqs.get(iID));
                in.put(jID, seqs.get(jID));
                float d = -1 * pA.pairAlign(in, gap);
                Map<String, Integer> dmiID = dm.get(iID);
                dmiID.put(jID, d);
                dm.put(iID, dmiID);
                Map<String, Integer> dmjID = dm.get(jID);
                dmjID.put(iID, d);
                dm.put(jID, dmjID);
            }
        }

        String bestPairiID = null;
        String bestPairjID = null;
        float bestPaird = Float.MAX_VALUE;
        for (int i = 0; i < IDs.size() - 1; ++i) {
            String iID = IDs.get(i);
            for (int j = i + 1; j < IDs.size(); ++i) {
                String jID = IDs.get(j);
                float d = dm.get(iID).get(jID);
                if (((bestPairiID == null) && (bestPairjID == null)) || (d < bestPaird)) {
                    bestPairiID = iID;
                    bestPairjID = jID;
                    bestPaird = d;
                }
            }
        }
        String bestThirdID = null;
        float bestThirdd = Float.MAX_VALUE;
        for (String ID : IDs) {
            float d = Float.MAX_VALUE;
            if(!bestPairiID.equals(ID) && !bestPairjID.equals(ID)) {
                d = dm.get(bestPairiID).get(ID) + dm.get(bestPairjID).get(ID);
            }
            if((bestThirdID == null) || (d < bestThirdd)) {
                bestThirdID = ID;
                bestThirdd = d;
            }
        }
        Map<String, String> in = new HashMap<String, String>();
        in.put(bestPairiID, seqs.get(bestPairiID));
        in.put(bestPairjID, seqs.get(bestPairjID));
        in.put(bestThirdID, seqs.get(bestThirdID));
        TripleAlign tA = new TripleAlign();
        Map<String, String> out = tA.tripleAlign(in, gap);

        PriorityQueue<Entry> pq = new PriorityQueue<>();

        for (String ID : IDs) {
            if (!out.containsKey(ID)) {
                pq.add(new Entry(dm.get(bestPairiID).get(ID), ID));
                pq.add(new Entry(dm.get(bestPairjID).get(ID), ID));
                pq.add(new Entry(dm.get(bestThirdID).get(ID), ID));
            }
        }

        while (out.keySet().size() < IDs.size()) {
            Entry dID = pq.poll();
            String ID = dID.getValue();
            if (!out.containsKey(ID)) {
                continue;
            }
            mergeAlign(out, ID, seqs.get(ID), gap);
            for (String newID : IDs) {
                if (!out.containsKey(newID)) {
                    pq.add(new Entry(dm.get(newID).get(ID), newID));
                }
            }
        }

        return out;
    }
}