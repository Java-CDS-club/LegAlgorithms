import java.util.ArrayList;

public class SCAlgorithms {
    
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
     * Creates an array of steady value (course or speed) intervals from an 
     * array of input values. Based on statistics and calculated limits 
     * for deviation from the mean value. The maximal deviation is compared.
     * @param  input_array the array of input values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @return double[] - array of speeds extracted form totes
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> fifo_mean_st_maxdev(double[] input_array, int minelements) {
        
        ArrayList<SpanPair> periods = new ArrayList<>();
        
        int istart=0, iend=1;
        while ((istart <= input_array.length - minelements) && (iend < input_array.length)) {
            boolean bcondition = false;
            if (++iend < input_array.length) {
                int numelements = iend - istart;
                final double dmax = SCStatistics.max(input_array, istart, iend);
                final double dmin = SCStatistics.min(input_array, istart, iend);
                if (dmax - dmin < 0.01)
                    bcondition = true;
                else { // The 2nd test statistic is compared to the quantile of Student's S(0,1) distribution with inumelements-2 degrees of freedom
                    final double dstdev = SCStatistics.stdev(input_array, istart, iend);
                    if (dstdev < 0.01) // Preventing devzerro errors
                        bcondition = true;
                    else {
                        final double dmean = SCStatistics.mean(input_array, istart, iend);
                        final double ddeviation = Math.max(dmax - dmean, dmean - dmin);
                        final double drange = Math.abs(input_array[iend-1] - input_array[istart]);
                        final double dtest1 = drange * Math.sqrt(2.0) / dstdev;
                        final double dtest2 = ddeviation / dstdev;
                        final double dquantile1 = SCStatistics.get99StudentQuantil(numelements - 1);
                        final double dquantile2 = numelements > 2 ? SCStatistics.get999StudentQuantil(numelements - 2) : dquantile1;
                        bcondition = dtest1 <= dquantile1 && dtest2 <= dquantile2;
                    }
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
