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
}
