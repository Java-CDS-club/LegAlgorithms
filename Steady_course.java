
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
            System.out.println("There are " + steady_CourseAndSpeed_intervals.size() + " steady course intervals.");
            for(int ii=0; ii<steady_CourseAndSpeed_intervals.size(); ii++) {
                int index_start = steady_CourseAndSpeed_intervals.get(ii).first;
                int index_end = steady_CourseAndSpeed_intervals.get(ii).second;
                int numelements = index_end - index_start;

                double[] timesreg = SCStatistics.getRelativeTimes(totes, index_start, index_end);
                double delapsedtime = timesreg[timesreg.length - 1] - timesreg[0];

                String spath = "course_" + index_start + "_" + index_end + ".txt";
                SCFileReader.writeHeadings2File(totes, index_start, index_end, spath, true, true);

                System.out.println("index: " + index_start + " - " + index_end + 
                                   "  numelements = " + numelements + 
                                   "  elapsed time = " + delapsedtime);
            }
    	}

    	catch(Exception e) {
            System.out.println(e);
    	}
    }
}


