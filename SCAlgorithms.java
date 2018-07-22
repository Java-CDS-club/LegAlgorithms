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
     * for deviation from the mean value. The maximal deviation is compared to the test statistics.
     * Inclination of regression line is analyzed afterwards.
     * @param  times the array of times
     * @param  values the array of values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @param bRegressionAnalysis perform regression analysis if true
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
            if (SCConstants.SPEED_STEADY_RANGE >= (dmax - dmin) || SCConstants.SPEED_STEADY_STDEV >= dstdev) 
                bcondition = true;
            else {
                bcondition = areDeviationsInAllowedLimits(values, istart, iend);
            }

            if(!bcondition || values.length == iend) {

                if(iend - istart == minelements) // if no success; if no steady interval at the very beginning - then, iterate forward.
                    iend = ++istart + minelements;

                else {
                    boolean cond_regression = bRegressionAnalysis ? true : false;

                    // Is horizontal line - linear regressiona analysis. If not, iterate backward until finding horizontal line.
                    if(bRegressionAnalysis) {
                        do {
                            cond_regression = isRegressionLineHorizontal(times, values, istart, iend, SCConstants.SPEED_STEADY_RANGE);
                            } while (!cond_regression && iend-- > istart + minelements);
                        }

                    if(!bRegressionAnalysis || cond_regression) {
                        int ishift = shiftIntervalRight(times, values, istart, iend, bRegressionAnalysis);
                        istart += ishift;
                        iend += ishift;

                        SpanPair sp = new SpanPair(istart, (iend != values.length ? --iend : iend));
                        periods.add(sp);
                        
                        istart = iend;
                        iend += minelements;
                    }
                    else // if no interval passed regression test - iterate forward
                        iend = ++istart + minelements;
                }

            }
            else
                iend++;
        }

        return periods;
    }
    
    //-------------------------------------------------------------------------
    /**
     * Returns a number of possible steady-interval shifting (to the right) so that the new interval should be statistically better.
     * (Example: the initial interval [5,15) is statistically good; [5,16) is bad; [6,16) is good and even better
     * than [5,15) - it is recommendable to shift the initial steady interval to the right for 1)
     * @param  times the array of times
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @return int - Return a number of possible (and recommended) shifting to the right
     */
    public static int shiftIntervalRight(double[] times, double[] values, int istart, int iend, boolean bRegressionAnalysis) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.shiftIntervalRight - 'istart' index (" + istart + ") is less than 0.");

        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.shiftIntervalRight - 'istart' (" + istart + ") is bigger than 'iend' (" + iend + ").");

        if (times.length != values.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.fifo_mean_st_maxdev - input arrays of different length");

        if (iend > values.length - 1)
            return 0;

        // next interval is not proper if it doesn't fulfill the basic condition
        if(!areDeviationsInAllowedLimits(values, istart+1, iend+1))
            return 0;

        int numelements = iend - istart;
        double stdev_old = SCStatistics.stdev(values, istart, iend);

        // shifted aray
        double stdev_new = SCStatistics.stdev(values, istart+1, iend+1);
        if(stdev_new >= stdev_old)
            return 0;

        if(bRegressionAnalysis) {
            double dmax = SCStatistics.max(values, istart, iend);
            double dmin = SCStatistics.min(values, istart, iend);
            if (dmax - dmin > 0.1) { // no need to analyze regression of the next interval if all the values in array are within 0.01 range
                SCStatistics.RegrResults res = SCStatistics.lregression(times, values, istart+1, iend+1);
                double dtest = Math.abs(res.a / res.ma);
                double dquantile = numelements > 2 ? SCStatistics.get999StudentQuantil(numelements - 2) : SCStatistics.get999StudentQuantil(1) ;
                double regression_grow = Math.abs(res.a * (times[iend] - times[istart+1]));
                boolean bcond = (regression_grow <= stdev_new || dtest <= dquantile);
                if(!bcond) { // no shift if the next interval regression line is more inclined
                    SCStatistics.RegrResults res0 = SCStatistics.lregression(times, values, istart, iend);
                    if(Math.abs(res0.a) < Math.abs(res.a))
                        return 0;
                }
            }
        }

        // the steady interval should be shifted for 1 or even more (recursive call)
        return 1 + shiftIntervalRight(times, values, istart+1, iend+1, bRegressionAnalysis);
    }

    //-------------------------------------------------------------------------
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
    
    //-------------------------------------------------------------------------
    /**
     * Intersects two lists of intervals (e.g. steady course and steady speed intervals).
     * @param  first the first list of intervals (the list has to be ordered)
     * @param  second the second list of intervals (the list has to be ordered)
     * @return rrayList<SpanPair> - Returns new list of a combined/intersected intervals.
     */
    public static ArrayList<SpanPair> intersectLists(ArrayList<SpanPair> first, ArrayList<SpanPair> second) {
        ArrayList<SpanPair> intersection = new ArrayList<>();

        // examine and treat properly all the topological possibilities

        int index_second = 0; // index of a SpanPair element from the second list
        Iterator<SpanPair> iter = first.iterator();
        while (iter.hasNext()) { // iterate through the first list
            SpanPair pair1 = iter.next();
            for(int ii=index_second; ii<second.size(); ii++) { // iterate through the second list, starting at index_second
                SpanPair pair2 = second.get(ii);
                if(pair1.first >= pair2.first && pair1.first < pair2.second) { // the beginning of pair1 falls between pair2 interval
                    intersection.add( new SpanPair(pair1.first, Integer.min(pair1.second, pair2.second)) );
                    index_second = pair1.second < pair2.second ? ii : ii+1;
                    if(pair1.second <= pair2.second)
                    break;
                }
                else if(pair1.second > pair2.first && pair1.second <= pair2.second) {// the end of pair1 falls between pair2 interval
                    intersection.add( new SpanPair(Integer.max(pair1.first, pair2.first), pair1.second) );
                    index_second = ii;
                    break;
                }
                else if(pair1.first <= pair2.first && pair1.second >= pair2.second) // if pair1 envelopes pair2
                    intersection.add( new SpanPair(pair2.first, pair2.second) );
                else if(pair1.second <= pair2.first) // pair1 and pair2 are disjunct intervals (the first pair precedes the seconds)
                    break;
               // pair1 and pair2 are disjunct intervals (the second pair precedes the first) - nothing to be done
            }
        }

        return intersection;
    }
}
