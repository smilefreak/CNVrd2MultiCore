#include <Rcpp.h>

using namespace Rcpp;

/**
 * Function stretches a matrix to the window size from a segment_scores
 * list.
 * 
 * @author James Boocock
 * @date 17/06/2013
 * @return a 
 */

// [[Rcpp::export]]

NumericMatrix getSubRegionMatrixFromSegmentScores(List segmentResults,int window_size,int st, int en){
    // calculate the matrix size
    int list_length = segmentResults.size(); 
    int mat_size = (en/window_size);
    NumericMatrix result(list_length,mat_size);
    int row_index = 0;
    int loc_index = 0;
    ////Rcout << " segment_length " << list_length << "\n";
    ////Rcout << " window_size " << window_size << "\n";
    ////Rcout << " st " << st << "\n";
    //Rcout << " en " << en << "\n";
    for(int list_index=0; list_index < list_length; list_index++){
        DataFrame temp_dat = as<DataFrame>(segmentResults[list_index]);
        NumericVector loc_start = as<NumericVector>(temp_dat["loc.start"]);
        NumericVector loc_end  =  as<NumericVector>(temp_dat["loc.end"]);
        NumericVector seg_mean = as<NumericVector>(temp_dat["seg.mean"]);
        //Rcout << "Can we extract these three items. \n";
        int end_pos = 0;
        int start_window = window_size;
        // Variables to work out what the distance between 2 points
        int to = 0;
        int from =0;
        int len_start = loc_start.size();
        for(int location_index= 0; location_index < len_start ; location_index++){
            end_pos = loc_end[location_index];
            while(start_window <= end_pos){
                start_window += window_size;
                to++;
            }
            for(int j = from; j < to; j++){
                result(list_index,j) = seg_mean(location_index);
            }
            from = to;
        }
        //result(list_index,to) = seg_mean(len_start - 1);
    }
   return result; 
}
