import java.io.IOException;
import java.util.Map;

public class Driver {
    public static void main(String args[]) throws IOException {
        Common c = new Common();
        Map<String, String> seqs = c.parseFasta();
    }
}
