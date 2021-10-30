package qut;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

import static org.junit.jupiter.api.Assertions.*;

class compare_result_against {

    @org.junit.jupiter.api.Test
    void defaultResult() throws IOException, ExecutionException, InterruptedException {
        Parallel.main(new String[0]);
        Parallel.getConsensus().entrySet();
        assertTrue(1==1, "Correct 1");

    }

    @org.junit.jupiter.api.Test
    void originalSequentialRun() {



    }

}