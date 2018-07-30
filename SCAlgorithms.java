
import java.util.Collections;
import java.util.ArrayList;
import java.util.Iterator;

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

    //-------------------------------------------------------------------------
    /**
     * Returns true if we all the deviations from mean value are within the statistically allowed limits
     * (hypothesis |x(i)-x'| > 0 is tested; alternative hypothesis is |x(i)-x'| = 0
     * test value (x(i)-x')/sigma(x)*sqrt(n-1)/sqrt(n) has Student(0,1) distribution with n-1 degrees of freedom)
     * @param  times the array of times
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     * @return boolean - true if we can consider regression line horizontal
     */
    static boolean areDeviationsInAllowedLimits(double[] values, int istart, int iend, double steady_stdev) {

        final double dstdev = SCStatistics.stdev(values, istart, iend);
        if (dstdev < steady_stdev) // Preventing devzero errors
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
     * Returns true if the maximal range in an interval is within the statistically allowed limits.
     * Test value (max(x)-min(x))/sigma(range) should be less then corresponding quantile
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     * @return boolean - true if maximal range is in the allowed limits
     */
    static boolean isMaxRangeInAllowedLimits(double[] values, int istart, int iend, double steady_stdev) {

        final double dstdev = SCStatistics.stdev(values, istart, iend);
        if (dstdev < steady_stdev) // Preventing devzero errors
            return true;

        final int numelements = iend - istart;
        final double dmax = SCStatistics.max(values, istart, iend);
        final double dmin = SCStatistics.min(values, istart, iend);
        final double dtest = (dmax - dmin) / Math.sqrt(dstdev);
        final double dquantile = SCStatistics.get99RangeQuantile(numelements);

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
     * @param  steady_range it will be considered that interval is steady if its max and min valueas are in this predefined range
     * @param  steady_stdev it will be considered that interval is steady if standard deviation is less then this predefined argument
     * @return ArrayList< SpanPair > - array of steady intervals extracted form totes
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> fifo_mean_st_maxdev(double[] times, double[] values, int minelements, boolean bRegressionAnalysis, double steady_range, double steady_stdev) {

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
            if (steady_range >= (dmax - dmin) || steady_stdev >= dstdev)
                bcondition = true;
            else {
                bcondition = areDeviationsInAllowedLimits(values, istart, iend, steady_stdev);
            }

            if(!bcondition || values.length == iend) {

                if(iend - istart == minelements) // if no success; if no steady interval at the very beginning - then, iterate forward
                    iend = ++istart + minelements;

                else {
                    boolean cond_regression = bRegressionAnalysis ? true : false;

                    // Is horizontal line - linear regressiona analysis. If not, iterate backward until finding horizontal line.
                    if(bRegressionAnalysis) {
                        do {
                            cond_regression = isRegressionLineHorizontal(times, values, istart, iend-1, steady_range);
                            } while (!cond_regression && iend-- > istart + minelements);
                        }

                    if(!bRegressionAnalysis || cond_regression) {
                        if(shiftDevIntervalRight(times, values, istart, iend-1, bRegressionAnalysis, steady_stdev)) {
                            istart++;
                            continue;
                        }

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
     * Creates an array of steady value (course or speed) intervals from an
     * array of input values. Based on statistics and calculated allowed ranges. The maximal range is compared to the test statistics.
     * Optionally - the inclination of regression line is analyzed afterwards.
     * @param  times the array of times
     * @param  values the array of values
     * @param minelements the minimal number of elements in the 'steady' interval
     * @param bRegressionAnalysis perform regression analysis if true
     * @param  steady_range it will be considered that interval is steady if its max and min valueas are in this predefined range
     * @param  steady_stdev it will be considered that interval is steady if standard deviation is less then this predefined argument
     * @return ArrayList< SpanPair > - array of steady intervals extracted from totes
     * todo min_elapsedtime instead of minelements
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> fifo_maxrange(double[] times, double[] values, int minelements, boolean bRegressionAnalysis, double steady_range, double steady_stdev) {
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
            if (steady_range >= (dmax - dmin) || steady_stdev >= dstdev)
                bcondition = true;
            else {
                bcondition = isMaxRangeInAllowedLimits(values, istart, iend, steady_stdev);
            }

            if(!bcondition || values.length == iend) {

                if(iend - istart == minelements) // if no success; if no steady interval at the very beginning - then, iterate forward
                    iend = ++istart + minelements;

                else {
                    boolean cond_regression = bRegressionAnalysis ? true : false;

                    // Is horizontal line - linear regression analysis. If not, iterate backward until finding horizontal line.
                    if(bRegressionAnalysis) {
                        do {
                            cond_regression = isRegressionLineHorizontal(times, values, istart, iend, steady_range);
                           } while (!cond_regression && iend-- > istart + minelements);
                        }

                    if(!bRegressionAnalysis || cond_regression) {
                        boolean bshift = shiftRangeIntervalRight(times, values, istart, iend, bRegressionAnalysis, steady_stdev);
                        if(bshift) {
                            System.out.println("dev-shifted: " + istart + " - " + iend + "  (" + (iend-istart) + ")");
                            istart++;
                            continue;
                        }

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
     * @param bRegressionAnalysis perform regression analysis if true
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     * @return boolean - true if (deviations) interval should be shifted
     */
    public static boolean shiftDevIntervalRight(double[] times, double[] values, int istart, int iend, boolean bRegressionAnalysis, double steady_stdev) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.shiftIntervalRight - 'istart' index (" + istart + ") is less than 0.");

        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.shiftIntervalRight - 'istart' (" + istart + ") is bigger than 'iend' (" + iend + ").");

        if (times.length != values.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.fifo_mean_st_maxdev - input arrays of different length");

        if (iend >= values.length - 1)
            return false;

        // start,end of right shifted interval
        int istartnew = istart + 1;
        int iendnew = iend + 1;

        // next interval is not proper if it doesn't fulfill the basic condition
        if(!areDeviationsInAllowedLimits(values, istartnew, iendnew, steady_stdev))
            return false;

        int numelements = iend - istart;
        double stdev_old = SCStatistics.stdev(values, istart, iend);

        // shifted aray
        double stdev_new = SCStatistics.stdev(values, istartnew, iendnew);
        if(stdev_new >= stdev_old)
            return false;

        if(bRegressionAnalysis) {
            double dmax = SCStatistics.max(values, istartnew, iendnew);
            double dmin = SCStatistics.min(values, istartnew, iendnew);
            if (dmax - dmin > 0.1) { // no need to analyze regression of the next interval if all the values in array are within 0.01 range
                SCStatistics.RegrResults res = SCStatistics.lregression(times, values, istartnew, iendnew);
                double dtest = Math.abs(res.a / res.ma);
                double dquantile = numelements > 2 ? SCStatistics.get999StudentQuantil(numelements - 2) : SCStatistics.get999StudentQuantil(1) ;
                double regression_grow = Math.abs(res.a * (times[iend] - times[istart+1])); // iend (not iendnew). The end of intervals are always exclusive.
                boolean bcond = (regression_grow <= stdev_new || dtest <= dquantile);
                if(bcond) { // no shift if the next interval regression line is more inclined
                    SCStatistics.RegrResults res0 = SCStatistics.lregression(times, values, istart, iend);
                    if(Math.abs(res0.a) < Math.abs(res.a))
                        return false;
                }
            }
        }

        // If we came here, the interval should be shifted
        return true;
    }

    //-------------------------------------------------------------------------
    /**
     * Returns a number of possible steady-interval shifting (to the right) so that the new interval should be statistically better (satisfying maximal range condition).
     * (Example: the initial interval [5,15) is statistically good; [5,16) is bad; [6,16) is good and even better
     * than [5,15) - it is recommendable to shift the initial steady interval to the right for 1)
     * @param  times the array of times
     * @param  values the array of values
     * @param  istart (inclusive) lower index of the sub-array
     * @param  iend (exclusive) upper index of the sub-array
     * @param bRegressionAnalysis perform regression analysis if true
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     * @return int - Return a number of possible (and recommended) shifting to the right
     */
    public static boolean shiftRangeIntervalRight(double[] times, double[] values, int istart, int iend, boolean bRegressionAnalysis, double steady_stdev) {
        if (0 > istart)
            throw new IndexOutOfBoundsException("Calling ScStatistics.shiftIntervalRight - 'istart' index (" + istart + ") is less than 0.");

        if (istart > iend)
            throw new NegativeArraySizeException("Calling ScStatistics.shiftIntervalRight - 'istart' (" + istart + ") is bigger than 'iend' (" + iend + ").");

        if (times.length != values.length)
            throw new IndexOutOfBoundsException("Calling ScStatistics.fifo_mean_st_maxdev - input arrays of different length");

        if (iend > values.length - 1)
            return false;

        // start,end of right shifted interval
        int istartnew = istart + 1;
        int iendnew = iend + 1;

        // next interval is not proper if it doesn't fulfill the basic condition
        if(!isMaxRangeInAllowedLimits(values, istartnew, iendnew, steady_stdev))
            return false;

        int numelements = iend - istart;
        double stdev_old = SCStatistics.stdev(values, istart, iend);

        // shifted aray
        double stdev_new = SCStatistics.stdev(values, istartnew, iendnew);
        if(stdev_new >= stdev_old)
            return false;

        if(bRegressionAnalysis) {
            double dmax = SCStatistics.max(values, istartnew, iendnew);
            double dmin = SCStatistics.min(values, istartnew, iendnew);
            if (dmax - dmin > 0.1) { // no need to analyze regression of the next interval if all the values in array are within 0.01 range
                SCStatistics.RegrResults res = SCStatistics.lregression(times, values, istartnew, iendnew);
                double dtest = Math.abs(res.a / res.ma);
                double dquantile = numelements > 2 ? SCStatistics.get999StudentQuantil(numelements - 2) : SCStatistics.get999StudentQuantil(1) ;
                double regression_grow = Math.abs(res.a * (times[iendnew] - times[istartnew]));
                boolean bcond = (regression_grow <= stdev_new || dtest <= dquantile);
                if(bcond) { // no shift if the next interval regression line is more inclined
                    SCStatistics.RegrResults res0 = SCStatistics.lregression(times, values, istart, iend);
                    if(Math.abs(res0.a) < Math.abs(res.a))
                        return false;
                }
            }
        }

        // If we came here, the interval should be shifted
        return true;
    }

    //-------------------------------------------------------------------------
    /**
     * Merge intervals which can be considered same (the same mean, in statistical sense)
     * Criterion that has to be fulfilled - all the deviations (in the merged interval) have to be within allowed limits
     * @param  values the array of values
     * @param  periods_in the array of steady course/speed intervals
     * todo it should be optimized. It can run significantly faster.
     */
    public static ArrayList<SpanPair> mergeDevIntervals(double[] values, ArrayList<SpanPair> periods_in) {

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
     * Merge intervals which can be considered same (the same mean, in statistical sense)
     * Criterion that has to be fulfilled - the maximal range (in merged interval) has to be within allowed limits
     * @param  values the array of values
     * @param  periods_in the array of steady course/speed intervals
     * todo it should be optimized. It can run significantly faster
     */
    public static ArrayList<SpanPair> mergeRangeIntervals(double[] values, ArrayList<SpanPair> periods_in) {
        // todo: throw here...
        ArrayList<SpanPair> periods_out = new ArrayList<>();

        SpanPair in = periods_in.get(0);
        SpanPair out = new SpanPair(in.first, in.second);
        //int num1 = out.second - out.first;

        for(int jj=1; jj<periods_in.size(); jj++) {
            SpanPair merging = periods_in.get(jj);

            // todo: Fischer's test of dispersion equality

            // big-interval standard deviation
            final double stdev_big = SCStatistics.stdev(values, out.first, merging.second);
            final double max_big = SCStatistics.max(values, out.first, merging.second);
            final double min_big = SCStatistics.min(values, out.first, merging.second);

            final double dtest = (max_big - min_big) / Math.sqrt(stdev_big);
            final double dquantile = SCStatistics.get99RangeQuantile(merging.second - merging.first);

            boolean cond = (dtest <= dquantile);

            if(cond) {
                out.second = merging.second;
                if(periods_in.size()-1 == jj) // last steady interval in the input arraylist
                    periods_out.add(out);

                continue;
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

    //-------------------------------------------------------------------------
    /**
     * Adjusts touching intervals so that the sum of squares of deviations should have a minimal possible value.
     * (Optionally, it can be checked whether the regression line is horizontal)
     * Criteria that must be fulfilled - the deviations have to be within allowed limits
     * @param times the times in the examined period
     * @param values the values being examined (headings or speeds)
     * @param periods the array of steady course/speed intervals
     * @param minelements the minimal number of elements in a 'steady' interval
     * @param bRegressionAnalysis perform regression analysis if true
     * @param  steady_range it will be considered that interval is steady if its max and min values are in this predefined range
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     */
    public static void adjustDevTouchingIntervals(double[] times, double[] values, ArrayList<SpanPair> periods, int minelements, boolean bRegressionAnalysis, double steady_range, double steady_stdev) {

        // Iterate through the list, having in focus only neighboring intervals
        for(int ii=0, jj=1; jj<periods.size(); ii++, jj++) {
            SpanPair left = periods.get(ii);
            SpanPair right = periods.get(jj);

            // Treat only touching intervals
            if(left.second != right.first)
                continue;

            int touching_final = left.second;
            double minsumdev2 = Double.MAX_VALUE; // sum of squares of deviations - the crucial parameter

            // Test all possible touching points between...
            for(int touching = left.first + minelements; touching <= right.second - minelements; touching++) {

                // First, check whether deviations are in allowed limits
                boolean cond_first_ok = areDeviationsInAllowedLimits(values, left.first, touching, steady_stdev);
                if(!cond_first_ok) continue;
                boolean cond_second_ok = areDeviationsInAllowedLimits(values, touching, right.second, steady_stdev);
                if(!cond_second_ok) continue;

                // Optionally, verify that the regression lines are horizontal
                if(bRegressionAnalysis) {
                    cond_first_ok &= isRegressionLineHorizontal(times, values, left.first, touching, steady_range);
                    if(!cond_first_ok) continue;
                    cond_second_ok &= isRegressionLineHorizontal(times, values, touching, right.second, steady_range);
                    if(!cond_second_ok) continue;
                }

                // Is this is the optimal case, so far - remember it
                double stdev_left = SCStatistics.stdev(values, left.first, touching);
                double stdev_right = SCStatistics.stdev(values, touching, right.second);
                double sumdev2 = stdev_left*(touching - left.first - 1) + stdev_right*(right.second - touching - 1);

                if(sumdev2 < minsumdev2) {
                    minsumdev2 = sumdev2;
                    touching_final = touching;
                }
            }

            left.second = right.first = touching_final;
        }
    }

    //-------------------------------------------------------------------------
    /**
     * Adjusts touching intervals so that the sum of squares of deviations should have a minimal possible value.
     * (Optionally, it can be checked whether the regression line is horizontal)
     * Criteria that must be fulfilled - the max range has to be within allowed limits
     * @param times the times in the examined period
     * @param values the values being examined (headings or speeds)
     * @param periods the array of steady course/speed intervals
     * @param minelements the minimal number of elements in a 'steady' interval
     * @param bRegressionAnalysis perform regression analysis if true
     * @param  steady_range it will be considered that interval is steady if its max and min values are in this predefined range
     * @param  steady_stdev it will be considered that deviations are
     * in allowed limits if standard deviation is less then this predefined argument
     */
    //public static ArrayList<SpanPair> merge_intervals(double[] values, ArrayList<SpanPair> periods_in) {
    public static void adjustRangeTouchingIntervals(double[] times, double[] values, ArrayList<SpanPair> periods, int minelements, boolean bRegressionAnalysis, double steady_range, double steady_stdev) {

        // Iterate through the list, having in focus only neighboring intervals
        for(int ii=0, jj=1; jj<periods.size(); ii++, jj++) {
            SpanPair left = periods.get(ii);
            SpanPair right = periods.get(jj);

            // Treat only touching intervals
            if(left.second != right.first)
                continue;

            int touching_final = left.second;
            double minsumdev2 = Double.MAX_VALUE; // sum of squares of deviations - the crucial parameter

            // Test all possible touching points between...
            for(int touching = left.first + minelements; touching <= right.second - minelements; touching++) {

                // First, check whether deviations are in allowed limits
                boolean cond_first_ok = isMaxRangeInAllowedLimits(values, left.first, touching, steady_stdev);
                if(!cond_first_ok) continue;
                boolean cond_second_ok = isMaxRangeInAllowedLimits(values, touching, right.second, steady_stdev);
                if(!cond_second_ok) continue;

                // Optionally, verify that the regression lines are horizontal
                if(bRegressionAnalysis) {
                    cond_first_ok &= isRegressionLineHorizontal(times, values, left.first, touching, steady_range);
                    if(!cond_first_ok) continue;
                    cond_second_ok &= isRegressionLineHorizontal(times, values, touching, right.second, steady_range);
                    if(!cond_second_ok) continue;
                }

                // Is this is the optimal case, so far - remember it
                double stdev_left = SCStatistics.stdev(values, left.first, touching);
                double stdev_right = SCStatistics.stdev(values, touching, right.second);
                double sumdev2 = stdev_left*(touching - left.first - 1) + stdev_right*(right.second - touching - 1);

                if(sumdev2 < minsumdev2) {
                    minsumdev2 = sumdev2;
                    touching_final = touching;
                }
            }

            left.second = right.first = touching_final;
        }
    }

    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady-speed intervals from an input array of totes.
     * @param  totes given ArrayList of Tote objects
     * @return ArrayList< SpanPair > - array of steady-speed intervals extracted from totes
     */
    public static ArrayList<SpanPair> extractStaedySpeeds(ArrayList<Tote> totes) {

        double[] times = SCStatistics.getRelativeTimes(totes);
        double[] values = SCStatistics.getSpeeds(totes);

        int minelements = 35;

        // Apply max-deviation + regression check algorithm
        ArrayList<SCAlgorithms.SpanPair> speed_intervals0 = fifo_mean_st_maxdev(times, values, minelements, true, SCConstants.SPEED_STEADY_RANGE, SCConstants.SPEED_STEADY_STDEV);

        // Remove redundant intervals (peaks, holes). They all have non-homogeneous variance.

        // Find them...
        ArrayList<SCStatistics.Variance> list = new ArrayList<>();
        for(int ii=0; ii<speed_intervals0.size(); ii++) {
            int index_start = speed_intervals0.get(ii).first;
            int index_end = speed_intervals0.get(ii).second;
            int numelements = index_end - index_start;
            double dstdev = SCStatistics.stdev(values, index_start, index_end);
            list.add(new SCStatistics.Variance(ii, numelements-1,  dstdev*dstdev));
        }

        ArrayList<Integer> _4remove = SCStatistics.isolateNonHomogeneous(list, SCConstants.SPEED_STEADY_STDEV);

        // ... and remove them
        Collections.sort(_4remove);
        for(int jj = _4remove.size()-1; jj >= 0; jj--) {
            Integer myint = _4remove.get(jj);
            speed_intervals0.remove(myint.intValue());
        }

        // Adjust touching steady-course intervals (criteria: sum of squares of deviations = min)
        adjustTouchingIntervals(times, values, speed_intervals0, minelements, true, SCConstants.SPEED_STEADY_RANGE, SCConstants.SPEED_STEADY_STDEV);

        // Merge neighboring steady-speed intervals (those that pass statistical equal-means test)
        ArrayList<SCAlgorithms.SpanPair> speed_intervals2 = merge_intervals(values, speed_intervals0);

        return speed_intervals2;
    }

    //-------------------------------------------------------------------------
    /**
     * Creates an array of steady-course intervals from an input array of totes.
     * @param  totes given ArrayList of Tote objects
     * @return ArrayList< SpanPair > - array of steady-course intervals extracted from totes
     */
    public static ArrayList<SpanPair> extractStaedyHeadings(ArrayList<Tote> totes) {

        double[] times = SCStatistics.getRelativeTimes(totes);
        double[] values = SCStatistics.getHeadings(totes);

        int minelements = 35;

        // Apply max-deviation + regression check algorithm
        ArrayList<SCAlgorithms.SpanPair> course_intervals0 = fifo_mean_st_maxdev(times, values, minelements, true, SCConstants.COURSE_STEADY_RANGE, SCConstants.COURSE_STEADY_STDEV);

        // Adjust touching steady-course intervals (criteria: sum of squares of deviations = min)
        adjustTouchingIntervals(times, values, course_intervals0, minelements, true, SCConstants.COURSE_STEADY_RANGE, SCConstants.COURSE_STEADY_STDEV);

        // Merge neighboring steady-course intervals (those that pass statistical equal-means test)
        ArrayList<SCAlgorithms.SpanPair> course_intervals1 = merge_intervals(values, course_intervals0);

        return course_intervals1;
    }

}
