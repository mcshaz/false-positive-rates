#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#define fisherPIndx 0 // "p.value" - could use string keys, but aiming for micro optimisation
#define fisherCIIndx 1 // "conf.int"
#define fisherORIndx 2 // "estimate"

// [[Rcpp::export]]
DataFrame monteCarloFisher(int n, IntegerMatrix outcomes){
  // lots of ugly code follows - i should instantiate a struct with odds ratio etc 
  Function f("fisher.test");
  int aloc = n / 2;
  std::vector<int> lu(aloc * (aloc + 3) / 2 + 1);
  
  IntegerVector v(4);
  // Set the number of rows and columns to attribute dim of the vector object.
  v.attr("dim") = Dimension(2, 2);
  
  NumericMatrix p(outcomes.nrow(), outcomes.ncol() / 2);
  NumericMatrix ors(outcomes.nrow(), outcomes.ncol() / 2);
  NumericMatrix ci_lb(outcomes.nrow(), outcomes.ncol() / 2);
  NumericMatrix ci_ub(outcomes.nrow(), outcomes.ncol() / 2);
  for(int i = 0; i < outcomes.nrow(); ++i){
    for(int j = 0; j < outcomes.ncol(); j += 2){
      int a = outcomes[j * outcomes.nrow() + i];
      int b = outcomes[(j + 1) * outcomes.nrow() + i];
      int min;
      int max;
      if (a > b) {
        min = b;
        max = a;
      } else {
        min = a;
        max = b;
      }
      int outIndx = j * outcomes.nrow() / 2 + i;
      int luIndx = max * (max + 1)/2 + min;
      if (lu[luIndx]) {
        int copyIndx = lu[luIndx] - 1;
        p[outIndx] = p[copyIndx];
        if ((ors[copyIndx] < 1.0) != (a < b)) {
          
          ors[outIndx] = 1.0 / ors[copyIndx];
          ci_lb[outIndx] = 1.0 / ci_ub[copyIndx];
          ci_ub[outIndx] = 1.0 / ci_lb[copyIndx];
        } else {
          ors[outIndx] = ors[copyIndx];
          ci_lb[outIndx] = ci_lb[copyIndx];
          ci_ub[outIndx] = ci_ub[copyIndx];
        }
      } else {
        v[0] = a;
        v[1] = b;
        v[2] = aloc - a;
        v[3] = aloc - b;
        
        List fe = f(v);
        NumericVector pValue = fe[fisherPIndx];
        p[outIndx] = pValue[0];
        NumericVector odr = fe[fisherORIndx];
        ors[outIndx] = odr[0];
        NumericVector cis = fe[fisherCIIndx];
        ci_lb[outIndx] = cis[0];
        ci_ub[outIndx] = cis[1];
        lu[luIndx] = outIndx + 1;
      }
    }
  }
  return DataFrame::create( Named("or") = ors, Named("p") = p, Named("ci_lb") = ci_lb, Named("ci_ub") = ci_ub );
}
