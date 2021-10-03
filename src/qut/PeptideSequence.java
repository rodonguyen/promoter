package qut;

import jaligner.*;
import jaligner.matrix.*;
import java.io.*;


public class PeptideSequence 
{
    public byte[] bytes;

    private static final Matrix BLOSUM_62 = BLOSUM62.Load();

    public PeptideSequence()
    {
    }
    
    public PeptideSequence(String string)
    {
        bytes = string.getBytes();
    }
    
    public static double Similarity(PeptideSequence A, PeptideSequence B)
    {  
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore();       
    }
    
    @Override
    public String toString()
    {
        return new String(bytes);
    }
}
