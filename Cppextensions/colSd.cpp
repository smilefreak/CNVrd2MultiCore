#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
/**
 * Function takes a matrix and returns the standard deviation for each
 * col.
 *
 * @author James Boocock
 * @data 26/05/1989
 *
*/


// [[Rcpp::export]]
NumericVector colSd(NumericMatrix subRegionMatrix){
    int i,j;
    int col_len = subRegionMatrix.ncol();
    //Rcout << "Column length : " << col_len;
    int row_len = subRegionMatrix.nrow();
    //Rcout << "Row lencgth : " << row_len;
    double standard_dev, col_mean;
    NumericVector result(col_len);
    for(i = 0; i < col_len; i++){
        standard_dev = 0.0;
        col_mean = 0.0;
        for(j = 0; j < row_len;j++){
            col_mean += subRegionMatrix(j,i);
        }
        // calculate colMeans.
        col_mean /= row_len;
        //Rcout << "column mean = " << col_mean << "\n";
        
        for(j = 0; j < row_len;j++){
            standard_dev += pow((col_mean - subRegionMatrix(j,i)),2);
        }
        //Rcout << "standard_dev = " << standard_dev << "\n";
        // standard deviation calculation
        standard_dev = standard_dev / (row_len - 1);
        result[i] = sqrt(standard_dev); 
    }
   return result;
}
