import java.util.ArrayList;

public class SCAlgorithms {
    
	// if all the values in interval do not differ more, we can consider it the steady course/speed interval
	static double COURSE_STEADY_RANGE = 0.5;
	static double SPEED_STEADY_RANGE = 0.1;
    
	// if standard deviation is not greater, we can consider it the steady course/speed interval
	static double COURSE_STEADY_STDEV = 0.05;
	static double SPEED_STEADY_STDEV = 0.01;
    
    // internal utility data structures
    public static class SpanPair{
        public SpanPair (int _first, int _second) {
            first = _first;
            second = _second;
        }
        
        public int first;
        public int second;
    }
    
    public class KeyValuePair{
        public KeyValuePair (int _key, int _value) {
            key = _key;
            value = _value;
        }
        
        public int key;
        public int value;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady values (course or speed) intervals from an array of input values.
     * @param  input_array the array of input values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @param maxspan the maximal difference between the start value and end value
     * @return double[] - array of speeds extracted form totes
     */
    public static ArrayList<SpanPair> fifo_span(double[] input_array, int minelements, double maxspan) {
        
        ArrayList<SpanPair> periods = new ArrayList<>();
        
        int istart=0, iend=0;
        while (istart <= input_array.length - minelements && iend < input_array.length) {
            if(input_array.length == ++iend || Math.abs(input_array[iend] - input_array[istart]) > maxspan) {
                if (iend - istart > minelements) {
                    periods.add(new SpanPair(istart, iend));
                    istart = iend;
                }
                else
                    iend = ++istart;
            }
        }
                    
        return periods;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady values (course or speed) intervals from an 
     * array of input values. Based on statistically determined maximal allowed span value.
     * @param  input_array the array of input values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @return double[] - array of speeds extracted form totes
     */
    public static ArrayList<SpanPair> fifo_span_st(double[] input_array, int minelements) {
        
        ArrayList<SpanPair> periods = new ArrayList<>();
        
        int istart=0, iend=0;
        while (istart <= input_array.length - minelements && iend < input_array.length) {
            boolean bcondition = false;
            if (iend++ == istart) bcondition = Math.abs(input_array[iend] - input_array[istart]) < 2.5; // max jump, 10.0 for course, 2.5 for speed
            else if (input_array.length == iend) bcondition = true;
            else {
                double dstdev = SCStatistics.stdev(input_array, istart, iend);
                double drange = Math.abs(input_array[iend] - input_array[istart]);
                bcondition = drange * Math.sqrt(2.0) / dstdev <= SCStatistics.get99StudentQuantil(iend - istart);
            }
            if(false == bcondition) {
                if (iend - istart > minelements) {
                    periods.add(new SpanPair(istart, iend));
                    istart = iend;
                }
                else iend = ++istart;
            }
        }

        return periods;
    }
    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady value (course or speed) intervals from an 
     * array of input values. Based on statistics and calculated limits 
     * for deviation from the mean value.
     * @param  input_array the array of input values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @return double[] - array of speeds extracted form totes
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> fifo_mean_st(double[] input_array, int minelements) {
        
        ArrayList<SpanPair> periods = new ArrayList<>();
        
        int istart=0, iend=1;
        while ((istart <= input_array.length - minelements) && (iend < input_array.length)) {
            boolean bcondition = false;
            if (++iend < input_array.length) {
                int numelements = iend - istart;
                final double dmean = SCStatistics.mean(input_array, istart, iend);
                final double dstdev = SCStatistics.stdev(input_array, istart, iend);
                if (dstdev < 0.01) // Preventing devzerro errors
                    bcondition = true;
                else { // The 2nd test statistic is compared with the quantile of Student's S(0,1) distribution with inumelements-2 degrees of freedom
                    final double dval = input_array[iend];
                    final double ddeviation = Math.abs(dval - dmean);
                    final double drange = Math.abs(input_array[iend] - input_array[istart]);
                    final double dtest1 = drange * Math.sqrt(2.0) / dstdev;
                    final double dtest2 = ddeviation / dstdev;
                    final double dquantile1 = SCStatistics.get999StudentQuantil(numelements - 1);
                    final double dquantile2 = numelements > 2 ? SCStatistics.get99StudentQuantil(numelements - 2) : dquantile1;
                    bcondition = dtest1 <= dquantile1;
                    bcondition &= dtest2 <= dquantile2;
                }
            }

            if(false == bcondition) {
                if (iend - istart > minelements) {
                    periods.add(new SpanPair(istart, iend));
                    istart = iend;
                    iend = istart + 1;
                }
                else iend = ++istart + 1;
            }
        }

        return periods;
    }

    //-------------------------------------------------------------------------
    /**
     * Returns true if we all the deviations from mean value are within the statistically allowed limits 
     * (hypothesis |x(i)-x'| > 0 is tested; alternative hypothesis is |x(i)-x'| = 0
     * test value (x(i)-x')/sigma(x)*sqrt(n-1)/sqrt(n) has Student(0,1) distribution with n-1 degrees of freedom)
     * @param  times the array of times
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @param  steady_range the regression line will be considered horizontal if its grow is less then this predefined argument
     * @return boolean - true if we can consider regression line horizontal
     */
    static boolean areDeviationsInAllowedLimits(double[] values, int istart, int iend) {
        final double dstdev = SCStatistics.stdev(values, istart, iend);
        if (dstdev < SPEED_STEADY_STDEV) // Preventing devzero errors
            return true;
        
        final double dmax = SCStatistics.max(values, istart, iend);
        final double dmin = SCStatistics.min(values, istart, iend);
        final int numelements = iend - istart;
        final double dmean = SCStatistics.mean(values, istart, iend);
        final double ddeviation = Double.max(dmax - dmean, dmean - dmin);
        final double dtest = ddeviation / dstdev * Math.sqrt(numelements / (numelements-1));
        final double dquantile = numelements > 2 ? SCStatistics.get99StudentQuantil(numelements - 2) : SCStatistics.get99StudentQuantil(numelements - 1);
        
        return dtest <= dquantile;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Returns true if we can consider regression line horizontal 
     * (regression analysis: hypothesis |a| > 0 is tested; alternative hypothesis is a = 0
     * test value (a-0)/sigma(a) has Student(0,1) distribution with n-2 degrees of freedom)
     * @param  times the array of times
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @param  steady_range the regression line will be considered horizontal if it cannot grow more than this predefined argument
     * @return boolean - true if we can consider regression line horizontal
     */
    static boolean isRegressionLineHorizontal(double[] times, double[] values, int istart, int iend, double steady_range) {
        SCStatistics.RegrResults res = SCStatistics.lregression(times, values, istart, iend);
        final int numelements = iend - istart;
        final double dstdev = SCStatistics.stdev(values, istart, iend);
        final double regression_grow = Math.abs(res.a * (times[iend-1] - times[istart]));
        boolean cond1 = regression_grow <= Double.max(steady_range, dstdev);
        
        // if regression line cannot increase significntly we can consider the line horizontal 
        if (cond1) return true;
        
        // ... otherwise we have to compare test statistics
        boolean cond2 = true;
        if(res.ma > 1.e-8) { // Preventing devzero errors
            final double dtest = Math.abs(res.a / res.ma);
            final double dquantile = numelements > 2 ? SCStatistics.get999StudentQuantil(numelements - 2) : SCStatistics.get999StudentQuantil(1) ;
            cond2 = dtest <= dquantile;
        }
        
        return cond2;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady value (course or speed) intervals from an 
     * array of input values. Based on statistics and calculated limits 
     * for deviation from the mean value. The maximal deviation is compared to the test statistics..
     * @param  times the array of times
     * @param  values the array of values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @return ArrayList< SpanPair > - array of steady intervals extracted form totes
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> fifo_mean_st_maxdev(double[] times, double[] values, int minelements, boolean bRegressionAnalysis) {
        if (times.length != values.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.fifo_mean_st_maxdev - input arrays of different length");

        ArrayList<SpanPair> periods = new ArrayList<>();
        
        int istart=0;
        int iend = istart + minelements;
                        
        while ((istart <= values.length - minelements) && (iend <= values.length)) {
            boolean bcondition;
            final double dmax = SCStatistics.max(values, istart, iend);
            final double dmin = SCStatistics.min(values, istart, iend);
            final double dstdev = SCStatistics.stdev(values, istart, iend);
            
            // if all the values are in the predefined small range or 
            // if they deviate within the numbers precision (2 decimals, here)
            // we can drop out calculations and consider it to be steady course/speed interval
            if (SPEED_STEADY_RANGE >= (dmax - dmin) || SPEED_STEADY_STDEV >= dstdev) 
                bcondition = true;
            else {
                bcondition = areDeviationsInAllowedLimits(values, istart, iend);
                if(bcondition && bRegressionAnalysis)
                    bcondition &= isRegressionLineHorizontal(times, values, istart, iend, SPEED_STEADY_RANGE);
            }

            if(!bcondition || values.length == iend) {
                if (iend - istart != minelements) { // i.e. iend - istart > minelements
                    SpanPair sp = new SpanPair(istart, (iend != values.length ? --iend : iend));
                    periods.add(sp);
                    
                    istart = iend;
                    iend += minelements;
                }
                else
                    iend = ++istart + minelements;
            }
            else
                iend++;
        }

        return periods;
    }
    
    /**
     * Merge intervals which can be considered same (the same mean, in statistical sense)
     * @param  values the array of values
     * @param  periods_in the array of steady course/speed intervals
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> merge_intervals(double[] values, ArrayList<SpanPair> periods_in) {
        // todo: throw here...
        ArrayList<SpanPair> periods_out = new ArrayList<>();
        
        SpanPair in = periods_in.get(0);
        SpanPair out = new SpanPair(in.first, in.second);
        int num1 = out.second - out.first;
        
        
        for(int jj=1; jj<periods_in.size(); jj++) {
            SpanPair merging = periods_in.get(jj);
            double mean1 = SCStatistics.mean(values, out.first, out.second);
            double stdev1 = SCStatistics.stdev(values, out.first, out.second);
            
            int num2 = merging.second - merging.first;
            double mean2 = SCStatistics.mean(values, merging.first, merging.second);
            double stdev2 = SCStatistics.stdev(values, merging.first, merging.second);
            
            // todo: Fischer's test of dispersion equality
            
            // big-interval standard deviation
            double stdev_big = SCStatistics.stdev(values, out.first, merging.second);

            // standard deviation of the difference of two mean values
            double stedev_diff = Math.sqrt(stdev1*stdev1/num1 + stdev2*stdev2/num2);
            
            double dquantile_diff = SCStatistics.get99StudentQuantil(num1 + num2 - 2);
                
            boolean cond = Math.abs(mean1 - mean2) / stedev_diff <= dquantile_diff;
            
            if(cond) { // merge intervals
                int numelements = merging.second - out.first;
                double dquantile = SCStatistics.get99StudentQuantil(numelements - 2);
                
                int numdev = 0; // number of corrections beyond limits (d/md > quantile)
                for(int kk=out.first; kk<merging.second; kk++) {
                    double dval = values[kk];
                    double ddeviation = Math.abs(dval - mean2);
                    double dtest = ddeviation / stdev_big * Math.sqrt(numelements / (numelements-1));
                    if((dtest > dquantile) && (++numdev > 1 + numelements/100)) {
                        cond = false;
                        break;
                    }
                }
                
                if(cond) {
                    out.second = merging.second;
                    if(periods_in.size()-1 == jj) // last steady interval in the input arraylist
                        periods_out.add(out);
                    
                    continue;
                }
            }
            
            periods_out.add(out);
            out = new SpanPair(merging.first, merging.second);
            if(periods_in.size()-1 == jj) // last steady interval in the input arraylist
                periods_out.add(out);
        }
        
        return periods_out;
    }
}
