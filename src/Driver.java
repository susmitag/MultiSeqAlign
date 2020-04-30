import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class Driver {
    public static void main(String args[]) throws IOException {
        Common c = new Common();
        Map<String, String> seqs = c.parseFasta();
        MultiAlign m = new MultiAlign();
        float gap = (float)-3.32192809489;
        Map<String, String> aln = m.multiAlign(seqs, gap);
        List<String> IDs = new ArrayList<>(aln.keySet());
        Collections.sort(IDs);
        for (String ID : IDs) {
            System.out.println(">");
            System.out.println(ID);
            System.out.println();
            System.out.println(aln.get(ID));
            System.out.println();
        }
    }
}
