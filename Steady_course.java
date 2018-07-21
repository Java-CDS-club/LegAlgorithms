
import java.util.ArrayList;
import java.util.Iterator;

public class Steady_course {

    //-------------------------------------------------------------------------
    // internal utility data structure
    public static class SCPair {
        public SCPair (Tote _first, Tote _second) {
            first = _first;
            second = _second;
        }
        
        public Tote first;
        public Tote second;
    }
      
    //-------------------------------------------------------------------------
    // First try to extract steady-coarse periods
    // It's a very rough determination (without any statistics). Just first try; just the beginning...
    // (todo: min time interval as argument)
    public static ArrayList<SCPair> steady_coarse_periods(ArrayList<Tote> totes, int minelements, double djump) {
        
        ArrayList<SCPair> periods = new ArrayList<>();
        
        Tote first = totes.get(0);
        
        Iterator<Tote> iter = totes.iterator();
        iter.next(); // start iterating from the second element in array
        while (iter.hasNext()) {
            Tote this_tote = iter.next();
            if (Math.abs(this_tote.dheading - first.dheading) > djump) {
                if (this_tote.index - first.index >= minelements) periods.add(new SCPair(first, totes.get(this_tote.index-1)));

                first = this_tote;
            }
        }

        return periods;
    }
    
    //-------------------------------------------------------------------------
    // First try to extract steady-speed periods
    // It's a very rough determination (rather useless and without any 
    // statistics, at the moment). Just first try; just the beginning...
    // (todo: min time interval as argument)
    public static ArrayList<SCPair> steady_speed_periods(ArrayList<Tote> totes, int minelements, double djump) {
        
        ArrayList<SCPair> periods = new ArrayList<>();
        
        Tote first = totes.get(0);
        
        Iterator<Tote> iter = totes.iterator();
        iter.next(); // start iterating from the second element in array
        while (iter.hasNext()) {
            Tote this_tote = iter.next();
            if (Math.abs(this_tote.dspeed - first.dspeed) > djump) {
                if (this_tote.index - first.index >= minelements) periods.add(new SCPair(first, totes.get(this_tote.index-1)));

                first = this_tote;
            }
        }

        return periods;
    }

    //-------------------------------------------------------------------------
    // Test printing of totes. Elements between istart (inclusive)
    // and iend (exclusive) will be printed.
    static void print(ArrayList<Tote> totes, int istart, int iend) {
        
        // todo - validate istart and iend arguments
        double delapsed_max = totes.get(iend-1).dabsolute_time - totes.get(0).dabsolute_time;
        String str_format = delapsed_max < 10   ? "%3.1f" :
                            delapsed_max < 100  ? "%4.1f" :
                            delapsed_max < 1000 ? "%5.1f" :
                            delapsed_max < 1000 ? "%6.1f" :
                                                  "%7.1f" ;

        for(int ii=istart; ii<iend; ii++) { 
            Tote tote = totes.get(ii);
            String str_space = tote.index < 10 ? "    " :
                               tote.index < 100 ? "   " :
                               tote.index < 1000 ? "  " :
                               tote.index < 10000 ? " " :
                                                     "" ;

            double delapsed = tote.dabsolute_time - totes.get(0).dabsolute_time;

            System.out.println(tote.index + ")" + str_space + "time: "    + String.format("%.2f", tote.dabsolute_time) + "  " +
                                                              "elapsed: " + String.format(str_format, delapsed) + "  "  +
                                                              "course: "  + String.format("%.2f", tote.dheading) + "  " +
                                                              "speed: "   + String.format("%.2f", tote.dspeed) );
        }
    }

    
    //=========================================================================
    public static void main(String[] args) {
        
    	try {
            SCFileReader myreader = new SCFileReader();
            ArrayList<Tote> totes = myreader.process("OtherOwnship_Trimmed.txt");

            // printing test (only first 10 totes)
            // print(totes, 0, 10);
            
            // printing (roughly determined) steady-speed periods
            System.out.println();
            double[] times = SCStatistics.getRelativeTimes(totes);
            double[] values = SCStatistics.getSpeeds(totes);
            ArrayList<SCAlgorithms.SpanPair> speed_intervals0 = SCAlgorithms.fifo_mean_st_maxdev(times, values, 35, true);
            ArrayList<SCAlgorithms.SpanPair> speed_intervals = SCAlgorithms.merge_intervals(values, speed_intervals0);
            System.out.println("There are " + speed_intervals.size() + " steady speed intervals.");
            for(int ii=0; ii<speed_intervals.size(); ii++) {
                int index_start = speed_intervals.get(ii).first;
                int index_end = speed_intervals.get(ii).second;
                int numelements = index_end - index_start;
                double[] speeds = SCStatistics.getSpeeds(totes, index_start, index_end); 
                double dstartval = speeds[0];
                double dendval = speeds[index_end - index_start - 1];
                double dmean = SCStatistics.mean(speeds);
                double dmax = SCStatistics.max(speeds);
                double dmin = SCStatistics.min(speeds);
                double dstdev = SCStatistics.stdev(speeds);
                double dmaxdev = Math.max(dmax - dmean, dmean - dmin);
                double dtau = dmaxdev / dstdev * Math.sqrt(numelements / (numelements-1));
                double dquantile = SCStatistics.get99StudentQuantil(numelements - 2);

                double[] timesreg = SCStatistics.getRelativeTimes(totes, index_start, index_end);
                double[] speedsreg = SCStatistics.getSpeeds(totes, index_start, index_end);
                SCStatistics.RegrResults res = SCStatistics.lregression(timesreg, speedsreg, 0, numelements);
                double dtest3 = Math.abs(res.a / res.ma);
                double dquantile3 = SCStatistics.get999StudentQuantil(index_end - index_start - 2);
                //double grow = Math.abs(res.a * (timesreg[timesreg.length - 1] - timesreg[0]));
                double delapsedtime = timesreg[timesreg.length - 1] - timesreg[0];
                double grow = Math.abs(res.a * delapsedtime);

                if (true /*|| dtest3 <= dquantile3*/) {
                    String spath = "speed_" + index_start + "_" + index_end + ".txt";
                    SCFileReader.writeSpeeds2File(totes, index_start, index_end, spath, true, true);

                    System.out.println("index: " + index_start + " - " + index_end + 
                                       "  (" + (index_end-index_start) + ")" +
                                       "  mean=" + String.format("%.2f", dmean) + 
                                       "  stdev=" + String.format("%.3f", dstdev) +
                                       "  maxdev=" + String.format("%.2f", dmaxdev) +
                                       "  tau=" + String.format("%.2f", dtau) +
                                       "  quantile=" + String.format("%.2f", dquantile) +
                                       "  max=" + String.format("%.2f", dmax) + 
                                       "  min=" + String.format("%.2f", dmin) +
                                       "  start_val=" + String.format("%.2f", dstartval) + 
                                       "  end_val=" + String.format("%.2f", dendval) + 
                                       "    regression: a=" + String.format("%.8f", res.a) + 
                                       "  ma=" + String.format("%.8f", res.ma) + 
                                       "  test=" + String.format("%.3f", dtest3) + 
                                       "  quantil=" + String.format("%.3f", dquantile3) + 
                                       "  grow=" + String.format("%.3f", grow) + 
                                       (grow > Double.max(0.1, dstdev) && dtest3 > dquantile3 ? " ATTENTION!!" : ""));
                }
            }
    	}

    	catch(Exception e) {
            System.out.println(e);
    	}
    }
}


