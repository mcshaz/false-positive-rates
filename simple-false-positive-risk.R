rm(list = ls()) # clear all vars from the current workspace

setwd("~/GitHub/false-positive-rates")
# if(!exists("monteCarloFisher", mode="function")) 
    sourceCpp("monte-carlo-fisher.cpp")

a <- c(0.1)
simTrials <- function (
    monteCarloSims = 10000L,
    participantsPerArm = 100L,
    baselineRisk = 0.5,
    absoluteRR = 0.1) {
    
    ctrlOutcomes <- rbinom(monteCarloSims, participantsPerArm, baselineRisk)
    rxRisk <- baselineRisk - absoluteRR
    
    df <- monteCarloFisher(alloc = participantsPerArm,
                           outcomes = cbind(
                               rbinom(monteCarloSims, participantsPerArm, baselineRisk), # no dif interv
                               ctrlOutcomes, # control group
                               rbinom(monteCarloSims, participantsPerArm, rxRisk), # working Rx
                               ctrlOutcomes))
                                            
    
    sub <- paste0("2 arms both risk:", baselineRisk)
    hist(df$p.1, sub=sub)
    hist(df$p.1[df$p.1 < 0.1], sub=sub)
    
    sub <- paste0("control:", baselineRisk, "; intervention:", rxRisk)
    hist(df$p.2, sub=sub)
    hist(df$p.2[df$p.2 < 0.1], sub=sub)
    
    type1_tail1 <- sum(df$or.1 < 1.0 & df$p.1 < 0.05)
    type1_tail2 <- sum(df$or.1 > 1.0 & df$p.1 < 0.05)
    
    type2 <- sum(df$p.2 >= 0.05)
    type3 <- sum(df$or.2 > 1.0 & df$p.2 < 0.05)
    
    a <<- df
    
    resultOut <- function(text, count){
        return(paste0("\n\t", text, ": ", count, " (", 100 * count / monteCarloSims, "%)"))
    }
    cat("Simulated trials:", prettyNum(monteCarloSims, big.mark = " "),
        "\nAllocations per Trial Arm:", participantsPerArm,
        "\n-------------------------",
        "\nEqual Risks (", baselineRisk, "):", 
        resultOut("type 1 errors made: tail 1", type1_tail1),
        resultOut("type 1 errors made: tail 2", type1_tail2),
        resultOut("total type 1 errors made (FP)", type1_tail1 + type1_tail1),
        resultOut("includes type 1 with p < 0.02", sum(df$p.1 < 0.02)),
        "\nRisk Reduction (", sub, "):",
        resultOut("type 2 errors made", type2),
        resultOut("type 3 errors made", type3),
        resultOut("type 2 or 3 errors", type2 + type3),
        resultOut("no error (TP)", monteCarloSims - type2 - type3),"\n")
    
    cat("\n\n\nPoint effect estimate when type 1 error encountered (Odds Ratio) truth = 1.0:\n")
    print(summary(df$or.1[which(df$p.1 < 0.05 & df$or.1 > 1.0)]))
    print(summary(df$or.1[which(df$p.1 < 0.05 & df$or.1 < 1.0)]))
    cat("\nPoint effect estimate when no error (Odds Ratio) truth = ",
        rxRisk * (1 - baselineRisk) / (baselineRisk * (1 - rxRisk)), 
        ":\n")
    print(summary(df$or.2[which(df$p.2 < 0.05 & df$or.2 < 1.0)]))
    cat("confidence intervals of the odds ratios (lower then upper bounds)\n")
    print(summary(df$ci_lb.2[which(df$p.2 < 0.05 & df$or.2 < 1.0)]))
    print(summary(df$ci_ub.2[which(df$p.2 < 0.05 & df$or.2 < 1.0)]))
}

simTrials(
    monteCarloSims = 1000000L,
    participantsPerArm = 200L,
    baselineRisk = 0.5,
    absoluteRR = 0.1)
