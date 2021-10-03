package qut;

import java.io.*;
import java.util.*;

public class Gene
{
    public int strand;
    public int location;
    public String name;
    public PeptideSequence sequence;

    public Gene()
    {
    }

    public Gene(String name, int strand, int location, String sequence)
    {
        this.name = name;
        this.strand = strand;
        this.location = location;
        this.sequence = new PeptideSequence(sequence);
    }

    public static List<Gene> ParseGenes(BufferedReader reader) throws IOException
    {
        List<Gene> genes = new ArrayList<Gene>();   
        
        // loop through all features(genes) until we come to ORIGIN marking start of DNA
        while (true)
        {
            String line = reader.readLine();
            if (line == null)
                return null;

            String label = line.substring(0, Math.min(20, line.length())).trim();

            if ("CDS".equals(label))
            {
                int location = -1;
                int strand = 1;
                String GeneName = null, GeneID = null, translation = null;
                    
                String loc = line.substring(21);
                if (loc.startsWith("complement("))
                {
                    loc = loc.substring(11, loc.length()-1);
                    strand = -1;
                }

                if (loc.startsWith("join("))
                {
                    loc = loc.substring(5, loc.length()-1);
                    List<Integer> starts = new ArrayList<Integer>();
                    for (String part : loc.split(","))
                    {
                        if (part.startsWith("<"))
                        {
                            starts.add(-1);
                            break;
                        }
                        int separator = part.indexOf('.');
                        if (separator <= 0)
                            starts.add(-1);
                        else
                        {
                            String start = part.substring(0, separator);
                            int i = Integer.parseInt(start);
                            starts.add(i);
                        }
                    }
                    Collections.sort(starts);
                    location = starts.get(0);
                }
                else if (loc.startsWith("<"))
                    location = -1;
                else
                {
                    int separator = loc.indexOf('.');
                    String start = loc.substring(0, separator);
                    location = Integer.parseInt(start);
                }


                // loop through attributes of CDS (/db_xref, /translation, etc)
                while (true)
                {
                    line = reader.readLine();
                    label = line.substring(0,Math.min(line.length(), 20)).trim();
            
                    line = line.trim();
                    if ("".equals(label))
                    {
                        if (line.startsWith("/gene=\""))
                            GeneName = line.substring(7, line.length()-1);
                        else if (line.startsWith("/db_xref=\"GeneID:"))
                            GeneID = line.substring(10, line.length()-1);                                   
                        else if (line.startsWith("/translation="))
                        {
                            StringBuilder builder = new StringBuilder();
                            boolean end = line.endsWith("\"");
                            if (end)    
                                builder.append(line.substring(14, line.length()-1));
                            else
                            {
                                builder.append(line.substring(14));
                                while (!end)
                                {
                                    line= reader.readLine().trim();
                                    end = line.endsWith("\"");
                                    if (end)
                                        builder.append(line.substring(0, line.length()-1));
                                    else
                                        builder.append(line);
                                }
                            }
                            translation = builder.toString();
                        }
                    }
                    else // we've come to the end of the CDS
                        break;
                }
                if (GeneID != null && location > 0 && translation != null)
                    genes.add(new Gene(GeneID + "(" + GeneName + ")", strand, location, translation));
            }
            if ("ORIGIN".equals(label))
                return genes;
        }
    }
}
