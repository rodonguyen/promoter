package qut;

import edu.au.jacobi.alphabet.*;
import edu.au.jacobi.pattern.*;

public class Sigma70Definition 
{
    private static PWM getMinus10Pwm()
    {
        Alphabet alphabet = AlphabetDNA.instance();
        edu.au.jacobi.pattern.PWM pwm = new PWM(alphabet, 0.0);
        pwm.add('a', "-2.357, +1.744, -0.440, +1.344, +1.335, -2.480");
        pwm.add('c', "-1.496, -2.357, -0.960, -0.916, -0.752, -1.562");
        pwm.add('g', "-1.863, -2.762, -0.792, -1.101, -1.006, -2.762");
        pwm.add('t', "+1.667, -1.705, +1.119, -1.101, -1.372, +1.738");
        return pwm;
    }

    private static PWM getMinus35Pwm()
    {
        Alphabet alphabet = AlphabetDNA.instance();
        edu.au.jacobi.pattern.PWM pwm = new PWM(alphabet, 0.0);
        pwm.add('a', "-1.949, -1.782, -1.562, +1.129, -0.184, +1.191");
        pwm.add('c', "-1.705, -2.614, -1.372, -0.440, +1.009, -1.372");
        pwm.add('g', "-2.244, -2.357, +1.362, -1.632, -1.632, -0.605");
        pwm.add('t', "+1.688, +1.744, -0.504, -0.409, -0.349, -0.571");
        return pwm;
    }        

    public static Series getSeriesAll_Unanchored(double threshold)
    {
        PWM pwmM10 = getMinus10Pwm();
        PWM pwmM35 = getMinus35Pwm();
        Gap spacer = new Gap(14, 20, new double[] {0.00,0.10,0.15,0.36,1.00,0.31,0.15,0.10,0.00});

        pwmM10.setImpact(1.0);
        pwmM35.setImpact(0.9086);  //set weight = fractional range difference        
        spacer.setImpact(0.15);

        Series series = new SeriesAll(threshold);
        series.add(pwmM35);
        series.add(spacer);
        series.add(pwmM10);
        return series;
    }   	    
}
