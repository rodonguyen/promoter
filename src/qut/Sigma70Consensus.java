package qut;

import java.io.*;
import edu.au.jacobi.pattern.*;


public class Sigma70Consensus
{
    private int[][] _35Count = new int[6][4];
    private int[][] _10Count= new int[6][4];
    private int predictions;
    private int gapTotal;

    private int[] map = new int[128];
    private char[] alphabet = new char[4];

    public Sigma70Consensus()
    {
        gapTotal = 0;
        predictions = 0;
        
        map['A'] = 0; map['a'] = 0; alphabet[0] = 'A'; 
        map['C'] = 1; map['c'] = 1; alphabet[1] = 'C'; 
        map['G'] = 2; map['g'] = 2; alphabet[2] = 'G'; 
        map['T'] = 3; map['t'] = 3; alphabet[3] = 'T'; 
    }

    private double getGapAverage()
    {
        return (double) gapTotal / predictions;
    }

    public void addConsensus(Sigma70Consensus s)
    {
        predictions += s.predictions;
        gapTotal += s.gapTotal;
        for (int i=0; i<6; i++)
            for (int j=0; j<4; j++)
	    {
                _35Count[i][j] += s._35Count[i][j];
                _10Count[i][j] += s._10Count[i][j];
	    }
    }

    public void addMatch(Match match)
    {
        predictions += 1;
        gapTotal += match.getSubMatch(1).calcLength();
        add_10Sequence(match.getSubMatch(2).letters().getBytes());
        add_35Sequence(match.getSubMatch(0).letters().getBytes());
    }    
    
    public void set(Match match)
    {
        predictions = 0;
        gapTotal = 0;
        for (int i=0; i<6; i++)
            for (int j=0; j<4; j++)
	    {
                _35Count[i][j] = 0;
                _10Count[i][j] = 0;
	    }        
        
        
        addMatch(match);
    }
    
    private void add_35Sequence(byte[] sequence)
    {
        for (int pos = 0; pos < sequence.length; pos++)
            _35Count[pos][map[sequence[pos]]]++;
    }

    private void add_10Sequence(byte[] sequence)
    {
        for (int pos = 0; pos < sequence.length; pos++)
            _10Count[pos][map[sequence[pos]]]++;
    }
    
    private void Display(StringBuffer sb, String name, int[][] count)
    {
        sb.append(name);
        for (int i=0; i<count.length; i++)
        {
            int max = 0;
            int maxValue = count[i][0];
            for (int j=1; j<count[i].length; j++)
                if (count[i][j] > maxValue)
                {
                    max = j;
                    maxValue = count[i][j];
                }
            sb.append(alphabet[max]).append(' ');         
        }
    }
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer(" Consensus: ");
        Display(sb, "-35: ", _35Count);
        sb.append("gap: ");
        sb.append(String.format("%.1f", getGapAverage()));
        Display(sb, " -10: ", _10Count);        
        sb.append(" (" + predictions + " matches)");
        return sb.toString();
    }
}
