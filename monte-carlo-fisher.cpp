#include <Rcpp.h>
#include <algorithm>
#include <chrono>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#define fisherPIndx 0 // "p.value" - could use string keys, but aiming for micro optimisation
#define fisherCIIndx 1 // "conf.int"
#define fisherORIndx 2 // "estimate"

inline static int binomialCombinations(const int &alloc){
  // the grids can be the standard number of combinations of 2 numbers n * (n+1) /2, but 0 is a valid digit and also each identical pair i.e. 0,0 1,1 ...
  return alloc * (alloc + 3) / 2 + 1;
}

inline static int lookupIndex(const int &a, const int &b){
  int min;
  int max;
  if (a > b) {
    min = b;
    max = a;
  } else {
    min = a;
    max = b;
  }
  return max * (max + 1)/2 + min;
}

inline static std::pair<const int&,const int&> minmax(IntegerVector &outcomes, const int &allocs) {
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

class AbstractLookup
{
public:
  AbstractLookup(){}
  virtual ~AbstractLookup(){}
  virtual int get(const int& lookupIndx) = 0;
  virtual void set(const int& lookupIndx, const int& value) = 0;
};

class VectorLookup : public AbstractLookup
{
private:
  std::vector<int> _lookup;
  
public:
  VectorLookup(size_t s) : _lookup(s){}
  ~VectorLookup(){}
  int get(const int& lookupIndx);
  void set(const int& lookupIndx, const int& value);
};

int VectorLookup::get(const int& lookupIndx)
{
  return _lookup[lookupIndx];
}
void VectorLookup::set(const int& lookupIndx, const int& value)
{
  _lookup[lookupIndx] = value;
}

class MapLookup : public AbstractLookup
{
private:
  std::unordered_map<int, int> _lookup;
  
public:
  MapLookup(size_t s){
    _lookup.reserve(s);
  }
  ~MapLookup(){}
  int get(const int& lookupIndx);
  void set(const int& lookupIndx, const int& value);
};

int MapLookup::get(const int& lookupIndx)
{
  return _lookup[lookupIndx];
}
void MapLookup::set(const int& lookupIndx, const int& value)
{
  _lookup[lookupIndx] = value;
}

class CalcStore
{
private:
  AbstractLookup *_lookup;
public:
  CalcStore(std::pair<const int&, const int&> &mm);
  int get(const int& lookupIndx);
  void set(const int& lookupIndx, const int& value);
  ~CalcStore() {
    delete _lookup;
  }
};

CalcStore::CalcStore(std::pair<const int&, const int&> &mm) {
  if (mm.first < 16L) {
    int memAlloc = binomialCombinations(mm.second);
    _lookup = new VectorLookup(memAlloc);
  } else {
    int memAlloc = binomialCombinations(mm.second - mm.first);
    _lookup = new MapLookup(memAlloc);
  }
}

int CalcStore::get(const int &lookupIndx) {
  return _lookup->get(lookupIndx) - 1;
}
void CalcStore::set(const int& lookupIndx, const int& value){
  _lookup->set(lookupIndx, value + 1);
}


// [[Rcpp::export]]
DataFrame monteCarloFisher(int alloc, IntegerMatrix outcomes){
  // lots of ugly code follows - i should instantiate a struct with odds ratio etc 
  Function f("fisher.test");
  auto start = std::chrono::high_resolution_clock::now();
  auto mm = minmax(outcomes, alloc);
  auto finish = std::chrono::high_resolution_clock::now();        
  std::chrono::duration<double> elapsed = finish - start;
  Rcout << "Elapsed time: " << elapsed.count() * 10000000 << " s\n";
  CalcStore lu(mm);
  
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
      int luIndx = lookupIndex(a, b);
      int copyIndx = lu.get(luIndx);
      if (copyIndx != -1) {
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
        lu.set(luIndx, outIndx);
      }
    }
  }
  return DataFrame::create( Named("or") = ors, Named("p") = p, Named("ci_lb") = ci_lb, Named("ci_ub") = ci_ub );
}
