
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


