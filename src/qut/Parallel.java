package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.locks.ReentrantLock;

public class Parallel
{
    // Changed them to final
    private static final HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
//    private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    private static final ThreadLocal<Series> sigma70_pattern =  ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    // https://www.baeldung.com/java-threadlocal
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];
    static ReentrantLock lock = new ReentrantLock();
    // https://www.linkedin.com/learning/parallel-and-concurrent-programming-with-java-1/livelock-java-demo?u=57080313

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }

    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true)
        {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    private static boolean Homologous(PeptideSequence A, PeptideSequence B)
    {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene)
    {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
           upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion)
    {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir)
    {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir)
    {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }


    public synchronized void addConsensus(String name, Match prediction) {
        consensus.get(name).addMatch(prediction);
        consensus.get("all").addMatch(prediction);
    }
 /* **==========================================================================================**
    ||                                   PARALLEL STREAM CODE                                   ||
    ||                                          below                                           || */

   /**
     * Run app by utilizing parallelStream method with an initial step of preparing data
     */
    public void runParallelStream(String referenceFile, String dir) throws IOException {
        // Preparing Data and store in List<TaskHandler>
        List<TaskHandler> taskHandlers = new ArrayList<>();
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        for (Gene referenceGene : referenceGenes) {
            for (String filename : ListGenbankFiles(dir)) {
                GenbankRecord record = Parse(filename);
                for (Gene gene : record.genes) {
                    taskHandlers.add(new TaskHandler(referenceGene, gene, record));
                }
            }
        }
        System.out.println("parallelStream - Preparing data finished!\nComputing...\n");

        // Initialize Threads with ParallelStream() and compute
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(Runtime.getRuntime().availableProcessors()));
        taskHandlers.parallelStream()
            .filter(task -> Homologous(task.getGene().sequence, task.getReferenceGene().sequence))
            .forEach(task -> {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(task.getRecord().nucleotides, task.getGene());
                Match prediction = PredictPromoter(upStreamRegion);
                if (prediction != null) addConsensus(task.getReferenceGene().name, prediction);
            });

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

/* ||                                                                                          ||
   ||                                                                                          ||
   **==========================================================================================**
   ||                                  EXECUTOR SERVICE CODE                                   ||
   ||                                          below                                           ||*/

    /**
     * Run app by utilizing executorService method
     */
    public void runExecutorService(String referenceFile, String dir) throws IOException, ExecutionException, InterruptedException {
        // Set number of threads equal to number of threads available on machine
        ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<Future> futureTasks = new ArrayList<>();

        for (Gene referenceGene : referenceGenes) {
            for (String filename : ListGenbankFiles(dir)) {
                GenbankRecord record = Parse(filename);
//                List<Gene> genes = record.genes;
//                genes.parallelStream().forEach(gene -> {
                for (Gene gene : record.genes) {
                    Future futureTask = executorService.submit(new RunnableTask(referenceGene, gene, record));
                    futureTasks.add(futureTask);
                }
            }
        }
        System.out.println("executorService - Preparing data finished!\nComputing...\n");

        // Shut down executorService to stop accepting new tasks and to close threads as they finished
        executorService.shutdown();
        for (Future futureTask : futureTasks) futureTask.get();

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }
    /**
     * Runnable class used in 'runExecutorService'
     */
    public class RunnableTask implements Runnable {
        private final Gene referenceGene;
        private final Gene gene;
        private final GenbankRecord record;

        public RunnableTask(Gene referenceGene, Gene gene, GenbankRecord record) {
            this.referenceGene = referenceGene;
            this.gene = gene;
            this.record = record;
        }

        public void run() {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                Match prediction = PredictPromoter(upStreamRegion);
                if (prediction != null) {
                    addConsensus(referenceGene.name, prediction);
                }
            }
        }
    }
/* ||                                                                                          ||
   ||                                                                                          ||
   **==========================================================================================**
   ||                                  SEQUENTIAL CODE                                         ||
   ||                                       below                                              ||*/

    public static void run(String referenceFile, String dir) throws IOException
    {
        // Set fixed number of threads (thread pool)
        // https://www.baeldung.com/java-8-parallel-streams-custom-threadpool
        //   Avoid Thread pool being used by other components within the app
        //   memory leak -> remember to shutdown
        // https://stackoverflow.com/questions/21163108/custom-thread-pool-in-java-8-parallel-stream
        //   2nd answer
        // https://www.javacodemonk.com/java-8-parallel-stream-custom-threadpool-48643a91
        // System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "12");

        // Get referenceGene (which Ecoli genes will be compared to)
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        long startTime = System.currentTimeMillis();

        //  ************************* 1st For Loop *******************************
         for (Gene referenceGene : referenceGenes) {
            System.out.println(referenceGene.name);
            // ************************ 2nd For Loop ************************
            // Get Ecoli file in 'dir'
            List<String> filenames = ListGenbankFiles(dir);
//            filenames.parallelStream().forEach(filename -> {
             for (String filename : ListGenbankFiles(dir)) {
                System.out.println(filename);
                GenbankRecord record = null;
                try {
                    record = Parse(filename);
                } catch (IOException e) {
                    e.printStackTrace();
                }

                // ************************** 3rd For Loop *************************
		        // For each gene in the genes record
//                 List<Gene> genes = record.genes;
                 GenbankRecord finalRecord = record;
//                  genes.parallelStream().forEach(gene -> {
                 for (Gene gene : record.genes) {
                     // Compare gene from reference with gene from records
                     // Homologous: determine if 2 genes serve the same purpose,
                     //      using SmithWatermanGototh algorithm (expensive)
                     if (Homologous(gene.sequence, referenceGene.sequence)) {
                         // Extract upstreamRegion
                         NucleotideSequence upStreamRegion = GetUpstreamRegion(finalRecord.nucleotides, gene);
                         // Predict whether if it is a promoter
                         Match prediction = PredictPromoter(upStreamRegion);
                         if (prediction != null) {
                             // Store result in 'concensus'
                             consensus.get(referenceGene.name).addMatch(prediction);
                             consensus.get("all").addMatch(prediction);
//                             System.out.println("Finish a gene prediction.");
                         }
                     }
                 }
            }
        }

        // Print result from 'concensus'
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
           System.out.println(entry.getKey() + " " + entry.getValue());

        System.out.println("Execution time: " + (System.currentTimeMillis() - startTime) + " ms");
    }

    public static void main(String[] args) throws FileNotFoundException, IOException, ExecutionException, InterruptedException {
        int iteration = 1;
        long[] durations = new long[iteration];
        for (int i=0; i<iteration; i++) {
            long start = System.currentTimeMillis();
            new Parallel().runParallelStream("referenceGenes.list", "Ecoli");
//            new Sequential2Parallel().runParallelStream("referenceGenes.list", "Ecoli");
//            run("src/referenceGenes.list", "src/Ecoli");
            durations[i] = System.currentTimeMillis() - start;
        }


        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());

        System.out.println("Average execution time is " + Arrays.stream(durations).sum()/iteration/1000 + " s");
        System.out.println(
                "\n-------------------------------" +
                "\nSmall progress is still progress." +
                "\nYou will get there Rodo! Get the 7" +
                "\n-------------------------------");
    }

}
