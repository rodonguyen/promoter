package qut;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import edu.au.jacobi.sequence.Sequence;
import edu.au.jacobi.sequence.list.FeatureList;
import java.util.Comparator;

public class BioPatterns 
{
    public static Match getBestMatch(Series pattern, String inputString)
    {
        return new Sequence("DNA",inputString).searchBest(1, 1, pattern);
    }	        
}
