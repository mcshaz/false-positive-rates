#include <Rcpp.h>
#include <utility>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#define fisherPIndx 0 // "p.value" - could use string keys, but aiming for micro optimisation
#define fisherCIIndx 1 // "conf.int"
#define fisherORIndx 2 // "estimate"

inline static int binomialCombintnCnt(const int &alloc){
  // the grids can be the standard number of combinations of 2 numbers n * (n+1) /2, but 0 is a valid digit and also each identical pair i.e. 0,0 1,1 ...
  return alloc * (alloc + 3) / 2 + 1;
}

inline static int hashMinMax(const int &min, const int &max){
  return max * (max + 1)/2 + min;
}

inline static std::pair<int, int> findMinmax(IntegerVector &outcomes, const int &allocs) {
  IntegerVector::iterator it = outcomes.begin();
  int min;
  int max = min = *it;
  for(IntegerVector::iterator end = outcomes.end(); it != end; ++it)
  {
    if (*it < min) {
      min = *it;
    } else if (*it > max) {
      max = *it;
    } else {
      continue;
    }
    if (min == 0 && max == allocs) {
      break;
    }
  }
  return std::make_pair(min, max);
}

class VectorLookup
{
private:
  std::vector<int> _lookup;
  int _min;
  int lookupIndx(const int& a, const int& b);
  
public:
  VectorLookup(const int& min, const int& max) : _lookup(binomialCombintnCnt(max - min), -1), _min(min){};
  ~VectorLookup(){}
  int get(const int& a, const int& b);
  void set(const int& a, const int& b, const int& value);
};

int VectorLookup::lookupIndx(const int& a, const int& b) {
  return a > b 
    ? hashMinMax(b - _min, a - _min)
    : hashMinMax(a - _min, b - _min);
  
}

int VectorLookup::get(const int& a, const int& b)
{
  return _lookup[lookupIndx(a, b)];
}

void VectorLookup::set(const int& a, const int& b, const int& value)
{
  _lookup[lookupIndx(a, b)] = value;
}

// [[Rcpp::export]]
DataFrame monteCarloFisher(int alloc, IntegerMatrix outcomes){
  // lots of ugly code follows - i should instantiate a struct with odds ratio etc 
  Function f("fisher.test");
  auto mm = findMinmax(outcomes, alloc);
  VectorLookup lu(mm.first, mm.second);
  
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
      int outIndx = j * outcomes.nrow() / 2 + i;
      int copyIndx = lu.get(a, b);
      if (copyIndx != -1) {
        p[outIndx] = p[copyIndx];
        if ((ors[copyIndx] < 1.0) != (a <= b)) {
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
        v[2] = alloc - a;
        v[3] = alloc - b;
        
        List fe = f(v);
        NumericVector pValue = fe[fisherPIndx];
        p[outIndx] = pValue[0];
        NumericVector odr = fe[fisherORIndx];
        ors[outIndx] = odr[0];
        NumericVector cis = fe[fisherCIIndx];
        ci_lb[outIndx] = cis[0];
        ci_ub[outIndx] = cis[1];
        lu.set(a, b, outIndx);
      }
    }
  }
  return DataFrame::create( Named("or") = ors, Named("p") = p, Named("ci_lb") = ci_lb, Named("ci_ub") = ci_ub );
}
