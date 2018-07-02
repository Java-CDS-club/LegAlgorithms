
import java.util.ArrayList;
import java.util.Iterator;

/**
 * SCStatistics - utility methods of mathematical statistics. Various tools for extracting data patterns 
 * and making predictions about the vessels' positions.
 */
public class SCStatistics {
	
    //-------------------------------------------------------------------------
    /**
     * Creates a double array of speed values from an ArrayList of totes
     * @param  totes given ArrayList of Tote objects
     * @return double[] - array of speeds extracted form totes
     */
    public static double[] getSpeeds(ArrayList<Tote> totes) {
        int index = 0;
        double[] da = new double[totes.size()];
        Iterator<Tote> iter = totes.iterator();
        while (iter.hasNext())
            da[index++] = iter.next().dspeed;

        return da;
    }

    //-------------------------------------------------------------------------
    /**
     * Creates a double array of speed values from an ArrayList of totes within the given indexes
     * @param  totes given ArrayList of Tote objects
     * @param  istart (inclusive) lower index of the totes-array segment extracted for calculation
     * @param  iend (exclusive) upper index of the totes-array segment extracted for calculation
     * @return double[] - array of speeds extracted form totes
     */
    public static double[] getSpeeds(ArrayList<Tote> totes, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getSpeeds - 'istart' index (" + istart + ") is less than 0.");

        if (totes.size() < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getSpeeds - 'iend' index (" + iend + ") is greater than array length.");
        
        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.getSpeeds - 'istart' (" + istart + ") is less than 'iend' (" + iend + ").");
        
        int index = 0;
        double[] da = new double[iend - istart];
        for (int ii=istart; ii<iend; ii++)
            da[index++] = totes.get(ii).dspeed;

        return da;
    }

    //-------------------------------------------------------------------------
    /**
     * Creates a double array of heading values from an ArrayList of totes
     * @param  totes given ArrayList of Tote objects
     * @return double[] - array of headings extracted form totes
     */
    public static double[] getHeadings(ArrayList<Tote> totes) {
        double[] dl = new double[totes.size()];
        int index = 0;
        Iterator<Tote> iter = totes.iterator();
        while (iter.hasNext())
            dl[index++] = iter.next().dheading;

        return dl;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Creates a double array of sheading values from an ArrayList of totes within the given indexes
     * @param  totes given ArrayList of Tote objects
     * @param  istart (inclusive) lower index of the totes-array segment extracted for calculation
     * @param  iend (exclusive) upper index of the totes-array segment extracted for calculation
     * @return double[] - array of headings extracted form totes
     */
    public static double[] getHeadings(ArrayList<Tote> totes, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getHeadings - 'istart' index (" + istart + ") is less than 0.");

        if (totes.size() < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getHeadings - 'iend' index (" + iend + ") is greater than array length.");
        
        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.getHeadings - 'istart' (" + istart + ") is less than 'iend' (" + iend + ").");
        
        int index = 0;
        double[] da = new double[iend - istart];
        for (int ii=istart; ii<iend; ii++)
            da[index++] = totes.get(ii).dheading;

        return da;
    }

    //-------------------------------------------------------------------------
    /**
     * Returns maximal value from an array of doubles
     * @param  ardoubles given array
     * @return double - max value within the array
     */
    public static double max(double[] ardoubles) {
        double dmax = Double.MIN_VALUE;
        for (int ii=0; ii<ardoubles.length; ii++)
            dmax = Double.max(dmax, ardoubles[ii]);

        return dmax;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Returns minimal value from an array of doubles
     * @param  ardoubles given array
     * @return double - min value within the array
     */
    public static double min(double[] ardoubles) {
        double dmin = Double.MAX_VALUE;
        for (int ii=0; ii<ardoubles.length; ii++)
                dmin = Double.min(dmin, ardoubles[ii]);

        return dmin;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Returns mean from an array of normally distributed values
     * @param  ardoubles given array
     * @return double - mean value (calculated as simple arithmetic mean)
     */
    public static double mean(double[] ardoubles) {
        double dsum = 0.0;
        for (int ii=0; ii<ardoubles.length; ii++)
            dsum += ardoubles[ii];

        return dsum / ardoubles.length;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Returns mean from an array of values with different 'weight' (e.g. probability)
     * @param  ardoubles given array of values
     * @param  arweights given array of weights
     * @return double - mean value (calculated as weighted arithmetic mean)
     */
    public static double weighted_mean(double[] ardoubles, double[] arweights){
        if (ardoubles.length != arweights.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.weighted_mean - arrays of different length");
        
        double dsum = 0.0;
        double dsumweights = 0.0;
        for (int ii=0; ii<ardoubles.length; ii++) {
            dsum += ardoubles[ii] * arweights[ii];
            dsumweights += arweights[ii];
        }

        return dsum / dsumweights;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Returns standard deviation an array of normally distributed values (expected mean value is unknown)
     * @param  ardoubles given array
     * @return double - standard deviation (calculated as simple arithmetic mean)
     */
    public static double stdev(double[] ardoubles) {
        double dmean = mean(ardoubles);
        double dsum_squares = 0.0;
        for (int ii=0; ii<ardoubles.length; ii++)
            dsum_squares += ardoubles[ii]*ardoubles[ii];
        
        return Math.sqrt((dsum_squares - dmean*dmean * ardoubles.length) / (ardoubles.length - 1));
    }
    
    //-------------------------------------------------------------------------
    /**
     * Returns standard deviation an array of normally distributed values (please note: expected 'mean' value is known)
     * @param  ardoubles given array
     * @param  dmean the expected value of stochastic variable which is represented by array values
     * @return double - standard deviation (calculated as simple arithmetic mean)
     */
    public static double stdev_mn(double[] ardoubles, double dmean) {
        double dsum_squares = 0.0;
        for (int ii=0; ii<ardoubles.length; ii++)
            dsum_squares += ardoubles[ii]*ardoubles[ii];
        
        return Math.sqrt((dsum_squares - dmean*dmean * ardoubles.length) / ardoubles.length);
    }
    
    //-------------------------------------------------------------------------
    /**
     * Returns standard deviation from an array of values with different 'weight' (e.g. probability)
     * @param  ardoubles given array
     * @param  arweights given array of weights
     * @return double - standard deviation (calculated as simple arithmetic mean)
     */
    public static double weighted_stdev(double[] ardoubles, double[] arweights) {
        if (ardoubles.length != arweights.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.weighted_mean - arrays of different length");
	    
        double dmean = weighted_mean(ardoubles, arweights);
        double dsum_squares = 0.0;
        double dsumweights = 0.0;
        for (int ii=0; ii<ardoubles.length; ii++) {
            dsum_squares += ardoubles[ii]*ardoubles[ii] * arweights[ii];
            dsumweights += arweights[ii];
        }
        
        return Math.sqrt((dsum_squares - dmean*dmean * dsumweights) / dsumweights);
    }

    //-------------------------------------------------------------------------
    /**
     * Returns two sided quantile of Student's distribution with the 95% confidence interval
     * @param  inum number of measured data in dataset
     * @return double - quantile of Student's distribution
     * @see https://en.wikipedia.org/wiki/Student%27s_t-distribution
     * @todo confidence interval could be parameterized
     */
    static double get95StudentQuantil(int inum) {
        if (inum < 1) 
            throw new IllegalArgumentException("Calling ScStatistics.get95StudentQuantil - number of elements leass than 1");
        
        final double[] quantiles_1_30 = {
        12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 
        2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
        }; 
        
        final double[] quantiles_32_50 = { // increment = 2
        2.037, 2.032, 2.028, 2.024, 2.021, 2.018, 2.015, 2.013, 2.011, 2.009
        };

        final double[] quantiles_55_70 = { // increment = 5
        2.004, 2.000, 1.997, 1.994 
        };
        
        final double[] quantiles_80_100 = { // increment = 10
        1.990, 1.987, 1.984 
        };
        
        final double quantil_120 = 1.980;
        
        final double[] quantiles_150_300 = { // increment = 50
        1.976, 1.972, 1.970, 1.968 
        };
        
        final double[] quantiles_400_500 = { // increment = 100
        1.966, 1.965 
        };
        
        final double quantil_infinity = 1.960;
        
        if (inum < 30) return quantiles_1_30[inum-1];
        else if (inum < 32) return quantiles_1_30[29]; // 2.042
        else if (inum <= 50) return quantiles_32_50[(inum-32) / 2];
        else if (inum < 55) return quantiles_32_50[9]; // 2.009
        else if (inum <= 70) return quantiles_55_70[(inum-55) / 5];
        else if (inum < 80) return quantiles_55_70[3]; // 1.994
        else if (inum <= 100) return quantiles_80_100[(inum-80) / 10];
        else if (inum < 120) return quantiles_80_100[2]; // 1.984
        else if (inum < 150) return quantil_120; // 1.980
        else if (inum <= 300) return quantiles_150_300[(inum-150) / 50];
        else if (inum < 400) return quantiles_150_300[3]; // 1.968
        else if (inum < 500) return quantiles_400_500[0]; // 1.966
        else if (inum < 700) return quantiles_400_500[1]; // 1.965
        else
            return quantil_infinity; // 1.96
    }
    
    //-------------------------------------------------------------------------
    /**
     * Equal mean values statistical test. Returns true if two mean values 
     * can be consider stochastically equal (confidence interval - 95%).
     * The assumption is that the both datasets are homogenous in the sense of 
     * the measurement accuracy. (e.g. measured with the same sensor)
     * @param  dmean1 first mean
     * @param  dstdev1 standard deviation of measurements (values) in the first dataset
     * @param  inum1 number of elements in the first dataset
     * @param  dmean2 second  mean
     * @param  dstdev2 standard deviation of measurements (values) in the second dataset
     * @param  inum2 number of elements in the second dataset
     * @return  boolean - true if means can be considered equal (in the stochastic sense)
     */
    public static boolean areMeansEqual(double dmean1, double dstdev1, int inum1, 
                                        double dmean2, double dstdev2, int inum2) {
        
        final double quantil95 = get95StudentQuantil(Integer.min(inum1, inum2));
        
        double dstdev =  Math.sqrt(dstdev1*dstdev1 / inum1 + dstdev2*dstdev2 / inum2);
        double ddiff = Math.abs(dmean1 - dmean2);
        double dtest_value = ddiff / dstdev;
        
        return dtest_value <= quantil95;
    }
}	
