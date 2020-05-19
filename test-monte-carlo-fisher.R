library(microbenchmark)
library(Rcpp)
rm(list = ls()) # clear all vars from the current workspace

setwd("~/GitHub/false-positive-rates")
sourceCpp("monte-carlo-fisher.cpp")

outcomes <- as.matrix(expand.grid(0L:10L, 0L:10L)) # all permutations. R also has a t(combn(0:10, 2)) & then add cbind(0:10, 0:10)
outcomes <- rbind(outcomes, outcomes)
alloc <- 10L

applyRFisher <- function(){
  prow <- function(r) {
    f <- fisher.test(matrix(c(r[1], alloc - r[1], r[2], alloc - r[2]), 2))
    return(c(f$p.value, unname(f$estimate), f$conf.int[1], f$conf.int[2]))
  }
  r <- t(apply(outcomes,1,prow)) 
  return(r)
}

applyMCF <- function() {
  df <- monteCarloFisher(alloc = alloc, outcomes = outcomes)
  r <- cbind(df$p, df$or, df$ci_lb, df$ci_ub)
  return(r)
}

#note tolerance below - 2.792825e-05 4,6 and 4,8
stopifnot(isTRUE(all.equal(applyRFisher(), applyMCF(), tolerance=3e-5)))

stopifnot(isTRUE(all.equal(applyRFisher(), applyMCF())))

alloc <- 100L
outcomes <- cbind(rbinom(2000L, alloc, 0.5), # control group
                 rbinom(2000L, alloc, 0.6))


mbm <- microbenchmark(applyRFisher(), applyMCF(), times = 2L, check='equivalent')
print(mbm)
