
import java.util.ArrayList;
import java.util.Iterator;

public class Steady_course {

    //-------------------------------------------------------------------------
	static void printIntervals(ArrayList<SCAlgorithms.SpanPair> intervals) /*throws IOException*/ {
		Iterator<SCAlgorithms.SpanPair> iter = intervals.iterator();
		while(iter.hasNext()) {
			SCAlgorithms.SpanPair item = iter.next();
			System.out.println(String.format("%5d", item.first) + " " + String.format("%5d", item.second));
		}
	}

    //=========================================================================
    public static void main(String[] args) {
        
    	try {
            SCFileReader myreader = new SCFileReader();
            ArrayList<Tote> totes = myreader.process("OtherOwnship_Trimmed.txt");

            double mintime = 300.0; // 5min = 300sec

            // Finding out steady-course periods
            ArrayList<SCAlgorithms.SpanPair> course0_intervals = SCAlgorithms.extractSteadyHeadings(totes, mintime);
            System.out.println("\nSteady course intervals:");
            printIntervals(course0_intervals);

            // Finding out steady-speed periods
            ArrayList<SCAlgorithms.SpanPair> speed0_intervals = SCAlgorithms.extractSteadySpeeds(totes, mintime);
            System.out.println("\nSteady speed intervals:");
            printIntervals(speed0_intervals);
            
            // Combine them
            ArrayList<SCAlgorithms.SpanPair> steady_CourseAndSpeed_intervals = SCAlgorithms.intersectLists(course0_intervals, speed0_intervals);
            System.out.println("\nSteady course-speed combined intervals:");
            printIntervals(steady_CourseAndSpeed_intervals);

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


