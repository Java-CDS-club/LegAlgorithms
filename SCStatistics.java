
import java.util.ArrayList;
import java.util.Iterator;

/**
 * SCStatistics - utility methods of mathematical statistics. Various tools for extracting data patterns 
 * and making predictions about the vessels positions.
 */
public class SCStatistics {
	
    public static class RegrResults {
        int istart;
        int iend;
        
        double taumax; // maximal normalized deviation
        
        double sigma0; // standard deviation
        
        // estimated parameters of linear model: Yi = a*Xi + b
        double a;
        double b;
        
        // ... and their estimated standard error (i.e reciprocal accuracy)
        double ma;
        double mb;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Creates a double array of time values from an ArrayList of totes
     * @param  totes given ArrayList of Tote objects
     * @return double[] - array of speeds extracted from totes
     */
    public static double[] getTimes(ArrayList<Tote> totes) {
        return getTimes(totes, 0, totes.size());
    }

    //-------------------------------------------------------------------------
    /**
     * Creates a double array of absolute time values from an ArrayList of totes within the given indexes
     * @param  totes given ArrayList of Tote objects
     * @param  istart (inclusive) lower index of the totes-array segment extracted for calculation
     * @param  iend (exclusive) upper index of the totes-array segment extracted for calculation
     * @return double[] - array of absolute times extracted from totes
     */
    public static double[] getTimes(ArrayList<Tote> totes, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getTimes - 'istart' index (" + istart + ") is less than 0.");

        if (totes.size() < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getTimes - 'iend' index (" + iend + ") is greater than array length.");
        
        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.getTimes - 'istart' (" + istart + ") is less than 'iend' (" + iend + ").");
        
        int index = 0;
        double[] da = new double[iend - istart];
        for (int ii=istart; ii<iend; ii++)
            da[index++] = totes.get(ii).dabsolute_time;

        return da;
    }

    //-------------------------------------------------------------------------
    /**
     * Creates a double array of relative time values from an ArrayList of totes within the given indexes. The values are relative to the first time in list i.e the elapsed times from totes[istart] are recorded.
     * @param  totes given ArrayList of Tote objects
     * @return double[] - array of speeds extracted from totes
     */
    public static double[] getRelativeTimes(ArrayList<Tote> totes) {
        return getRelativeTimes(totes, 0, totes.size());
    }

    //-------------------------------------------------------------------------
    /**
     * Creates a double array of relative time values from an ArrayList of totes within the given indexes. The values are relative to the first time in list i.e the elapsed times from totes[istart] are recorded.
     * @param  totes given ArrayList of Tote objects
     * @param  istart (inclusive) lower index of the totes-array segment extracted for calculation
     * @param  iend (exclusive) upper index of the totes-array segment extracted for calculation
     * @return double[] - array of absolute times extracted from totes
     */
    public static double[] getRelativeTimes(ArrayList<Tote> totes, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getRelativeTimes - 'istart' index (" + istart + ") is less than 0.");

        if (totes.size() < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.getRelativeTimes - 'iend' index (" + iend + ") is greater than array length.");
        
        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.getRelativeTimes - 'istart' (" + istart + ") is less than 'iend' (" + iend + ").");
        
        int index = 0;
        double dreduce = totes.get(istart).dabsolute_time;
        double[] da = new double[iend - istart];
        for (int ii=istart; ii<iend; ii++)
            da[index++] = totes.get(ii).dabsolute_time - dreduce;

        return da;
    }

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
     * Returns maximal value from a sub-array within an array of doubles
     * @param  ardoubles given array
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return double - max value within the array
     */
    public static double max(double[] ardoubles, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.max - 'istart' index (" + istart + ") is less than 0.");

        if (ardoubles.length < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.max - 'iend' index (" + iend + ") is greater than array length.");

        double dmax = Double.MIN_VALUE;
        for (int ii=istart; ii<iend; ii++)
            dmax = Double.max(dmax, ardoubles[ii]);

        return dmax;
    }
	
    //-------------------------------------------------------------------------
    /**
     * Returns minimal value from a sub-array within an array of doubles
     * @param  ardoubles given array
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return double - min value within the array
     */
    public static double min(double[] ardoubles, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.min - 'istart' index (" + istart + ") is less than 0.");

        if (ardoubles.length < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.min - 'iend' index (" + iend + ") is greater than array length.");

        double dmin = Double.MAX_VALUE;
        for (int ii=istart; ii<iend; ii++)
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
     * Returns mean of an sub-array from an array of normally distributed values
     * @param  ardoubles given array
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return double - mean value (calculated as simple arithmetic mean)
     */
    public static double mean(double[] ardoubles, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.mean - 'istart' index (" + istart + ") is less than 0.");

        if (ardoubles.length < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.mean - 'iend' index (" + iend + ") is greater than array length.");

        double dsum = 0.0;
        for (int ii=istart; ii<iend; ii++)
            dsum += ardoubles[ii];

        return dsum / (iend - istart);
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
     * Returns standard deviation from a sub-array of normally distributed values (expected mean value is unknown)
     * @param  ardoubles given array
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return double - standard deviation (calculated as simple arithmetic mean)
     */
    public static double stdev(double[] ardoubles, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.stdev - 'istart' index (" + istart + ") is less than 0.");

        if (ardoubles.length < iend)
            throw new IndexOutOfBoundsException("Calling ScStatistics.stdev - 'iend' index (" + iend + ") is greater than array length.");

        int numelemets = iend - istart;
        double dmean = mean(ardoubles, istart, iend);
        double dsum_squares = 0.0;
        for (int ii=istart; ii<iend; ii++)
            dsum_squares += ardoubles[ii]*ardoubles[ii];
        
        return Math.sqrt((dsum_squares - dmean*dmean * numelemets) / (numelemets - 1));
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
     * @param  inum number of observed/measured data in dataset
     * @return double - quantile of Student's distribution
     * @see https://en.wikipedia.org/wiki/Student%27s_t-distribution
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
     * Returns two sided quantile of Student's distribution with the 99% confidence interval
     * @param  inum number of observed/measured data in dataset
     * @return double - quantile of Student's distribution
     */
    static double get99StudentQuantil(int inum) {
        if (inum < 1) 
            throw new IllegalArgumentException("Calling ScStatistics.get99StudentQuantil - number of elements leass than 1");
        
        final double[] quantiles_1_30 = {
        63.66, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
        3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
        2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
        }; 
        
        final double[] quantiles_32_50 = { // increment = 2
        2.739, 2.728, 2.720, 2.712, 2.705, 2.698, 2.692, 2.687, 2.682, 2.678
        };

        final double[] quantiles_55_70 = { // increment = 5
        2.668, 2.660, 2.654, 2.648
        };
        
        final double[] quantiles_80_100 = { // increment = 10
        2.639, 2.632, 2.626
        };
        
        final double quantil_120 = 2.617;
        
        final double[] quantiles_150_300 = { // increment = 50 
        2.609, 2.601, 2.596, 2.592
        };
        
        final double[] quantiles_400_500 = { // increment = 100
        2.588, 2.586 
        };
        
        final double quantil_infinity = 2.576;
        
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
     * Returns two sided quantile of Student's distribution with the 99.9% confidence interval
     * @param  inum number of observed/measured data in dataset
     * @return double - quantile of Student's distribution
     */
    static double get999StudentQuantil(int inum) {
        if (inum < 1) 
            throw new IllegalArgumentException("Calling ScStatistics.get999StudentQuantil - number of elements leass than 1");
        
        final double[] quantiles_1_30 = {
        636.6, 31.96, 12.92, 8.610, 6.869, 5.959, 5.408, 5.041, 4.718, 4.587,    
        4.437, 4.318, 4.221, 4.140, 4.073, 4.015, 3.965, 3.922, 3.883, 3.850,
        3.819, 3.792, 3.767, 3.745, 3.725, 3.707, 3.690, 3.674, 3.659, 3.646
        }; 
        
        final double[] quantiles_32_50 = { // increment = 2
        3.622, 3.601, 3.582, 3.566, 3.551, 3.538, 3.526, 3.515, 3.505, 3.496
        };

        final double[] quantiles_55_70 = { // increment = 5
        3.476, 3.460, 3.447, 3.435
        };
        
        final double[] quantiles_80_100 = { // increment = 10
        3.416, 3.402, 3.390
        };
        
        final double quantil_120 = 3.374;
        
        final double[] quantiles_150_300 = { // increment = 50 
        3.357, 3.340, 3.330, 3.323
        };
        
        final double[] quantiles_400_500 = { // increment = 100
        3.315, 3.310 
        };
        
        final double quantil_infinity = 3.291;
        
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
    
    //-------------------------------------------------------------------------
    /**
     * Returns the set of resulting parameters to analyze linear regression
     * @param  x given domain array
     * @param  y given value array
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return double - max value within the array
     */
    public static RegrResults lregression(double[] x, double[] y, int istart, int iend) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.lregression - 'istart' index (" + istart + ") is less than 0.");

        if (x.length < iend - istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.lregression - 'iend' index (" + iend + ") is greater than array length.");

        if (x.length != y.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.lregression - input arrays of different size.");

        RegrResults result = new RegrResults();
        
        result.istart = istart;
        result.iend = iend;
        
        // find out the estimated values of basic linear regression parameters (a, b)
        int num = iend - istart;	// number of elements in the arrays
        double N11 = 0.0;			// sum(x*x)
        double N12 = 0.0;                       // sum(x)
        double N22 = num;                       // 
        
        double n1 = 0.0;                        // -sum(x*y)
        double n2 = 0.0;                        // -sum(y)
        
        for(int ii=istart; ii<iend; ii++) {
            N11 += x[ii] * x[ii];
            N12 += x[ii];
            n1 -= x[ii] * y[ii];
            n2 -= y[ii];
        }
        
        double det = N11 * N22 - N12*N12; 
        if (det < 0) // N is positive definite matrix and 'det' should be always positive
            throw new IllegalStateException("Calling ScStatistics.lregression - negative determinant of positive definite matrix. It should not be possible.");
        
        double Qx11 =  N22 / det;
        double Qx12 = -N12 / det;
        double Qx22 =  N11 / det;
        		
        result.a = -(Qx11*n1 + Qx12*n2);        // -( sumxy * num   + sumx  * sumy ) / det;
        result.b = -(Qx12*n1 + Qx22*n2);        // -( sumx  * sumxy - sumxx * sumy ) / det;
        
        // accuracy, standard deviations, allowed errors...
        double sumvv = 0.0;
        double vmax = Double.MIN_VALUE; // redundant; just for tests
        for(int jj=istart; jj<iend; jj++) {
        	double yest = x[jj]*result.a + result.b;	// estimated (expected) value
        	double v = y[jj] - yest;					// residuals - measured value minus expected value
        	vmax = Double.max(vmax, Math.abs(v));
        	sumvv += v*v;
        }
        
        result.sigma0 = Math.sqrt(sumvv / (num - 2));
        
        result.ma = result.sigma0 * Math.sqrt(Qx11);
        result.mb = result.sigma0 * Math.sqrt(Qx22);
        
        
        result.taumax = Double.MIN_VALUE;
        for(int ii=0; ii<num; ii++) {
        	double v = y[ii] - (x[ii]*result.a + result.b);
        	double Qv = 1. - (Qx11*x[ii]*x[ii] + 2*Qx12*x[ii] + Qx22);
        	double tau = Math.abs(v) / (result.sigma0 * Math.sqrt(Qv));
        	result.taumax = Double.max(result.taumax, tau);
        }

        return result;
    }
}

/*
 * Linear regression in matrix notation (method of least squares):
 * 
 * Mathematical model: Y = a*X + b
 * Stochastical model: The accuracy of all measured Y is the same
 * 
 * 
 * Matrix of coefficients (Jacobian - Ai1 = ∂Yi/∂a, Ai2 = ∂Yi/∂b) and vector of residuals:
 *
 *     | X1    1  |        | -Y1   |
 *     | X2    1  |        | -Y2   |
 * A = | ..    .. |    f = |  ...  |
 *     | ..    .. |        |  ...  |
 *     | Xnum  1  |        | -Ynum |
 * 
 *
 * Matrix of normal equations:
 * 
 *                  | sumxx  sumx |
 * N = trans(A)*A = |             |
 *                  | sumx   num  |
 *
 *                                    
 * detN = num*sumxx - sumx*sumx
 * 
 *
 *          |  num   -sumx  |
 * inv(N) = |               | / detN
 *          | -sumx   sumxx |
 * 
 *
 *
 * Vector of cofactors:
 * 
 *                  | -sumxy |
 * n = trans(A)*f = |        |
 *                  | -sumy  |
 *
 *
 *
 * The estimation of the linear regression parameters:
 * 
 *      | a |               | sumxy*num  - sumx*sumy  |
 * x' = |   | = -inv(N)*n = |                         | / detN
 *      | b |               | sumxx*sumy - sumx*sumxy |
 *
 *
 *
 * The covariance matrix:
 * 
 *               |  num  -sumx  |
 * Qx = inv(N) = |              | / detN
 *               | -sumx  sumxx |
 *
 *
 * Standard deviation of estimated parameters:
 * 
 * ma = sigma0 * sqrt(Qx[11])
 * mb = sigma0 * sqrt(Qx[22])
 *
 *
 * Corrections:
 * 
 * v = A*x' + f
 * 
 * The covariance matrix of corrections:
 * 
 * Qv = E - A*Qx*trans(A)
 * 
 * Diagonal members of Qv:
 * 
 * Qv[ii] = 1 - ( X[i]^2 * Qx[11]  +  2*X[i]*Qx[12] + Qx[22] )
 * 
 * mv[i] = sigma0 * sqrt(Qv[ii])
 * 
 * tau[i] = v[i] / mv[i]  (should be less then corresponding quantile)
 * 
 */

