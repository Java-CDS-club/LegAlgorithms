
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
            print(totes, 0, 10);
            
            // printing the first 3 (roughly determined) steady-corse periods
            System.out.println();
            ArrayList<SCPair> course_periods = steady_coarse_periods(totes, 10, 0.5);
            for(int ii=0; ii<3; ii++) {
                System.out.println("steady_coarse_period: index1=" + course_periods.get(ii).first.index + "  index2=" + course_periods.get(ii).second.index);
                System.out.println("      in other words: time1=" + course_periods.get(ii).first.stime + "  time2=" + course_periods.get(ii).second.stime);
            }
            
            // printing the first 3 (roughly determined) steady-speed periods
            System.out.println();
            ArrayList<SCPair> speed_periods = steady_speed_periods(totes, 10, 0.1);
            for(int ii=0; ii<3; ii++) {
                System.out.println("steady_coarse_period: index1=" + speed_periods.get(ii).first.index + "  index2=" + speed_periods.get(ii).second.index);
                System.out.println("      in other words: time1=" + speed_periods.get(ii).first.stime + "  time2=" + speed_periods.get(ii).second.stime);
            }
    	}
    	
    	catch(Exception e) {
            System.out.println(e);
    	}
    }
}


/*
            // small trivial unit test
            System.out.println();

            double[] dspeeds = SCStatistics.getSpeeds(totes, 0, 10);
            System.out.println("max speed = " + SCStatistics.max(dspeeds));
            System.out.println("min speed = " + SCStatistics.min(dspeeds));
            System.out.println("speed mean = " + SCStatistics.mean(dspeeds));
            System.out.println("speed stdev = " + SCStatistics.stdev(dspeeds));

            double[] weights0 = new double[dspeeds.length];
            for (int jj=0; jj<weights0.length; jj++) weights0[jj]=1.0;
            System.out.println("speed weighted_stdev  = " + SCStatistics.weighted_stdev(dspeeds, weights0));

            System.out.println();

            double[] dheadings = SCStatistics.getHeadings(totes, 0, 10);
            System.out.println("max course = " + SCStatistics.max(dheadings));
            System.out.println("min course = " + SCStatistics.min(dheadings));
            System.out.println("course mean = " + SCStatistics.mean(dheadings));
            System.out.println("course stdev = " + SCStatistics.stdev(dheadings));
            System.out.println("speed weighted_stdev  = " + SCStatistics.weighted_stdev(dheadings, weights0));

            System.out.println();

            double[] dheadings26 = SCStatistics.getHeadings(totes, 2, 6);
            double[] weights = new double[dheadings26.length];
            weights[0] = 1.0; weights[1] = 0.5; weights[2] = 0.5; weights[3] = 1.0;
            System.out.println("mean course2 = " + SCStatistics.mean(dheadings26));
            System.out.println("weighte mean course2 = " + SCStatistics.weighted_mean(dheadings26, weights));
            System.out.println("course2 stdev = " + SCStatistics.stdev(dheadings26));
            System.out.println("course2 weighted stdev = " + SCStatistics.weighted_stdev(dheadings26, weights));

            System.out.println();
            
            double[] dheadings2 = SCStatistics.getHeadings(totes, 4, 8);
            double dmean2 = SCStatistics.min(dheadings2);
            double dstdev2 = SCStatistics.stdev(dheadings2);

            double[] dheadings3 = SCStatistics.getHeadings(totes, 22, 29);
            double dmean3 = SCStatistics.min(dheadings3);
            double dstdev3 = SCStatistics.stdev(dheadings3);
            
            boolean beq = SCStatistics.areMeansEqual(dmean3, dstdev3, dheadings3.length, dmean2, dstdev2, dheadings2.length);
            System.out.println(beq ? "equal means" : "means are not equal");

*/
