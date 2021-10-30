package qut;
import jdk.internal.access.JavaIOFileDescriptorAccess;
import org.junit.jupiter.api.BeforeEach;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import static org.junit.jupiter.api.Assertions.*;

class CompareResult {

    private HashMap<String, String> defaultConsensus;

    @BeforeEach
    void declareDefaultConsensus() {
        defaultConsensus = new HashMap<>();
        defaultConsensus.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        defaultConsensus.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        defaultConsensus.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        defaultConsensus.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        defaultConsensus.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        defaultConsensus.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        defaultConsensus.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        defaultConsensus.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        defaultConsensus.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");
    }

    /**
     * Compare Parallel consensus with consensus from analyzing default dataset
     */
    @org.junit.jupiter.api.Test
    void parallel_defaultResult() throws IOException, ExecutionException, InterruptedException {
        Parallel.main(new String[0]);
        for (Map.Entry<String, Sigma70Consensus> entry : Parallel.getConsensus().entrySet()) {
            assertEquals(defaultConsensus.get(entry.getKey()), entry.getValue().toString());
        }
    }

    /**
     * Compare ExplicitThreading consensus with consensus from analyzing default dataset
     */
    @org.junit.jupiter.api.Test
    void explicit_defaultResult() throws InterruptedException {
        ExplicitThreading.main(null);
        for (Map.Entry<String, Sigma70Consensus> entry : ExplicitThreading.getConsensus().entrySet()) {
            assertEquals(defaultConsensus.get(entry.getKey()), entry.getValue().toString());
        }
    }


    /**
     * Compare Parallel consensus with consensus from analyzing default dataset
     * Note: This test should only be run when input data set has been changed
     * as Sequential takes much time to run
     */
    @org.junit.jupiter.api.Test
    void parallel_sequentialResult() throws IOException, ExecutionException, InterruptedException {
        Parallel.main(null);
        Sequential.main(null);
        HashMap<String, Sigma70Consensus> sequentialConsensus = Sequential.getConsensus();
        for (Map.Entry<String, Sigma70Consensus> entry : Parallel.getConsensus().entrySet()) {
            assertEquals(sequentialConsensus.get(entry.getKey()).toString(), entry.getValue().toString());
        }
    }

    /**
     * Compare ExplicitThreading consensus with consensus from analyzing default dataset.
     * Note: This test should only be run when input data set has been changed
     * as Sequential takes much time to run
     */
    @org.junit.jupiter.api.Test
    void explicitThreading_sequentialResult() throws IOException, InterruptedException {
        ExplicitThreading.main(null);
        Sequential.main(null);
        HashMap<String, Sigma70Consensus> sequentialConsensus = Sequential.getConsensus();
        for (Map.Entry<String, Sigma70Consensus> entry : ExplicitThreading.getConsensus().entrySet()) {
            assertEquals(sequentialConsensus.get(entry.getKey()).toString(), entry.getValue().toString());
        }
    }
}





//            System.out.println("\ndefault: " + defaultConsensus.get(entry.getKey()));
//            System.out.println("actual: " + entry.getValue());
//            System.out.println("default: " + defaultConsensus.get(entry.getKey()).getClass().getSimpleName());
//            System.out.println("actual: " + entry.getValue().getClass().getSimpleName());
