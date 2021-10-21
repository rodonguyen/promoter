package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.locks.ReentrantLock;

class Sequential_for_explicit_threading implements Runnable {
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];


    // =================== New ====================
    private String ecoliFilename;
    private static String referenceFile;
    static ReentrantLock lock1 = new ReentrantLock();

    public Sequential_for_explicit_threading(String ecoliFilename, String referenceFile) {
        this.ecoliFilename = ecoliFilename;
        this.referenceFile = referenceFile;
    }
    // ============================================

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
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
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
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    // TODO:
    // - create 4 threads for 4 Ecoli
    // - tag ID: 1 2 3 4
    // - print the ID tagged with the dir
    // - process and print <ID> - <promoter found>
    // - parallelStream lib

    public void run() {
        // Get referenceGene (which Ecoli genes will be compared to)
        List<Gene> referenceGenes = null;
        try {
            referenceGenes = ParseReferenceGenes(referenceFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

//        // Get Ecoli file in 'dir'
//        for (String filename : ListGenbankFiles(dir))
//        {
            System.out.println(ecoliFilename);

        GenbankRecord record = null;
        try {
            record = Parse(ecoliFilename);
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (Gene referenceGene : referenceGenes)
            {
		        System.out.println(referenceGene.name);

		        // For each gene in the record
                for (Gene gene : record.genes)
                    // Compare gene from reference with gene from records
                    // Homologous: determine if 2 genes serve the same purpose,
                    //      using SmithWatermanGototh algorithm (expensive)
                    if (Homologous(gene.sequence, referenceGene.sequence))
                    {
                        lock1.lock();
                        // Extract upstreamRegion
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
	                    // Predict whether if it is a promoter
                        Match prediction = PredictPromoter(upStreamRegion);
                        if (prediction != null)
                        {
                            // Store result in 'concensus'
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                        lock1.unlock();
                    }
            }
//        }
        // Print result from 'concensus'
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
           System.out.println(entry.getKey() + " " + entry.getValue());
    }

//    public static void main(String[] args) throws FileNotFoundException, IOException
//    {
////        run("src/referenceGenes.list", "src/Ecoli");
//        //////////////////
//
//        // Get Ecoli filename in 'dir'
//        for (String ecoliFilename : ListGenbankFiles("src/Ecoli")) {
//            new Sequential(ecoliFilename).start();
//        }
//    }
}

class ExplicitThreading {
    static ReentrantLock lock1 = new ReentrantLock();

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

    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
        // Get Ecoli filename in 'dir'
        List<String> ecoliFilenames = ListGenbankFiles("src/Ecoli");

        long startTime = System.nanoTime();

        List<Thread> threads = new ArrayList<Thread>();
        List<String> listGenBankFiles = ListGenbankFiles("src/Ecoli");
        for ( int i=0; i < listGenBankFiles.size(); i++) {
            String ecoliFilename = listGenBankFiles.get(i);
            Runnable task = new Sequential_for_explicit_threading(ecoliFilename, "src/referenceGenes.list");
            Thread thread = new Thread(task);
            thread.setName(String.valueOf(i));
            thread.start();

            threads.add(thread);
        }

//        https://www.vogella.com/tutorials/JavaConcurrency/article.html
//        int running = 0;
//        do {
//            running = 0;
//            for (Thread thread : threads) {
//                if (thread.isAlive()) running++;
//            }
//            System.out.println("We have " + running + " running threads. ");
//            Thread.sleep(1000);
//        } while (running > 0);

        for (Thread thread : threads)  thread.join();

        long timeElapsed = System.nanoTime() - startTime ;
        System.out.println("Execution time: " + timeElapsed/1000000 + " ms");
    }
}
