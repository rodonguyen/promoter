package qut;

import java.io.*;


public class NucleotideSequence
{
    public byte[] bytes;

    private static byte[] complement = new byte['z'];

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }
    
    public NucleotideSequence()
    {
    }

    public NucleotideSequence(byte[] sequence)
    {
        this.bytes = sequence;
    }    
    
    public NucleotideSequence(String string)
    {
        bytes = string.getBytes();
    }
    
    public NucleotideSequence GetUpstreamRegion(Gene gene)
    {               
        int upStreamDistance = 250;
        if (gene.location <= upStreamDistance)
           upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }    
    
    @Override
    public String toString()
    {
        return new String(bytes);
    }
}
