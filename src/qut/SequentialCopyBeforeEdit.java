package qut;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import jaligner.BLOSUM62;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

import java.io.*;
import java.util.*;
import java.util.concurrent.locks.ReentrantLock;

public class SequentialCopyBeforeEdit
{
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

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
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
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

    // TODO:
    // - create 4 threads for 4 Ecoli
    // - tag ID: 1 2 3 4
    // - print the ID tagged with the dir
    // - process and print <ID> - <promoter found>

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException
    {
        // Get referenceGene (which Ecoli genes will be compared to)
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);

        // Get Ecoli file in 'dir'
        for (String filename : ListGenbankFiles(dir))
        {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
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
                    }
            }
        }

        // Print result from 'concensus'
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
           System.out.println(entry.getKey() + " " + entry.getValue());
    }

    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        run("src/referenceGenes.list", "src/Ecoli");
    }
}
