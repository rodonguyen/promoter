package qut;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class test {
    // https://www.javacodemonk.com/java-8-parallel-stream-custom-threadpool-48643a91
    // Try parallelStream with this example
    public static void testParallelOperation() {
        long start = System.currentTimeMillis();
//        IntStream s = IntStream.range(0, 20);
        List<Integer> s = Arrays.asList(1, 2, 3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "20");
        s.parallelStream().forEach(i -> {
            try {
                Thread.sleep(100);
            } catch (Exception ignore) {}
            System.out.print((System.currentTimeMillis() - start) + "ms ");
        });
        System.out.println("\nOverall time consumed: "+ (System.currentTimeMillis() - start)+" ms");
    }


    public static void main(String[] args) {

        testParallelOperation();
    }
}
