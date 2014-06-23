#include <Rcpp.h>

using namespace Rcpp;
/**
 * function stretches each row of a matrix to the window size.
 *
 * @author James Boocock
 * @data 10/06/2013
 *
*/


// [[Rcpp::export]]
NumericMatrix stretchToWindowSize(NumericMatrix subRegionMatrix,NumericVector start, NumericVector end, int window_size){
    int i,j,m;
    int len_start = start.size();
    int len_end = end.size();
  //  Rcout << len_start << " end =: " << len_end << "\n";
    int mat_row = subRegionMatrix.nrow();
    int mat_col = subRegionMatrix.ncol();
   // Rcout << "row = " << mat_row << " col = " << mat_col << "\n";
    int start_window = window_size;
    int max_item = end[len_end - 1];
//    Rcout << "Max position: " <<  max_item << "\n"; 
    NumericMatrix result(mat_row,max_item / window_size);
    int from = 0;
    // define end_pos which keeps track of the end position we are tracking up to
    int end_pos;
    int to =  0;
    for(i = 0; i < len_start; i++){
        end_pos = end[i];
    //    Rcout << "End pos " << end_pos << "\n";
        while(start_window <= end_pos ){
            start_window += window_size;
            to++;
        }
        for(j = from; j < to; j ++){
            for(m = 0; m < mat_row; m++){
                //#Rcout << m << "\n";
                result(m,j) = subRegionMatrix(m,i);
            }
       }
       //Rcout << "From = " << from << " To = " << to << "\n";
       from = to; 
       //Rcout << "Start_window = " << start_window << " end position = " << end_pos << "\n";
    }
    return result;
}
