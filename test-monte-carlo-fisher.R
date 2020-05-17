library(microbenchmark)

library(Rcpp)
setwd("~/GitHub/false-positive-rates")
if(!exists("monteCarloFisher", mode="function")) 
  sourceCpp("monte-carlo-fisher.cpp")

# todo shouldn't test with randomly generated values, however tests must not be unique, and values should not be uniformly distributed
# maybe use dput and paste

outcomes = cbind(rbinom(2000, 100, 0.5), # control group
                 rbinom(2000, 100, 0.6))


applyRFisher <- function(){
  prow <- function(r) {
    f <- fisher.test(matrix(c(r[1], 100 - r[1], r[2], 100 - r[2]), 2))
    return(c(f$p.value, unname(f$estimate), f$conf.int[1], f$conf.int[2]))
  }
  r <- t(apply(outcomes,1,prow)) # vapply to be fair
  return(r)
}

applyMCF <- function() {
  df <- monteCarloFisher(n = 200, outcomes = outcomes)
  r <- cbind(df$p, df$or, df$ci_lb, df$ci_ub)
  return(r)
}

mbm <- microbenchmark(applyRFisher(), applyMCF(), times = 2L, check='equivalent')
print(mbm)
