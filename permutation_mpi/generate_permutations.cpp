#include <Rcpp.h>

using namespace Rcpp;

/**
 * Function takes a matrix with both a subregion and a start and end vector and returns a full matrix
 *
 * @author James Boocock
 * @date 18/05/2014
 *
 **/
// [[Rcpp::export]]

NumericMatrix permuteMatrix(NumericMatrix toPermute, NumericVector reorded_array){
    int i, j, m;
    NumericMatrix result(toPermute.nrow(),toPermute.ncol());
    

}
