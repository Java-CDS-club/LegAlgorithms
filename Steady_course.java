
import java.util.ArrayList;

public class Steady_course {

    //=========================================================================
    public static void main(String[] args) {
        
    	try {
            SCFileReader myreader = new SCFileReader();
            ArrayList<Tote> totes = myreader.process("OtherOwnship_Trimmed.txt");

            // Finding out steady-course periods
            ArrayList<SCAlgorithms.SpanPair> course0_intervals = SCAlgorithms.extractStaedyHeadings(totes);

            // Finding out steady-speed periods
            ArrayList<SCAlgorithms.SpanPair> speed0_intervals = SCAlgorithms.extractStaedySpeeds(totes);
            
            // Combine them
            ArrayList<SCAlgorithms.SpanPair> steady_CourseAndSpeed_intervals = SCAlgorithms.intersectLists(course0_intervals, speed0_intervals);

            // Trace the results
            //System.out.println("There are " + (any list from above).size() + " steady course intervals.");
            //for(int ii=0; ii<(any list from above).size(); ii++) {
            //    int index_start = (any list from above).get(ii).first;
            //    int index_end = (any list from above).get(ii).second;
            //    int numelements = index_end - index_start;

            //    double[] timesreg = SCStatistics.getRelativeTimes(totes, index_start, index_end);
            //    double delapsedtime = timesreg[timesreg.length - 1] - timesreg[0];

            //    System.out.println("index: " + index_start + " - " + index_end + 
            //                       "  numelements = " + numelements + 
            //                       "  elapsed time = " + delapsedtime);
            //}
    	}

    	catch(Exception e) {
            System.out.println(e);
    	}
    }
}


