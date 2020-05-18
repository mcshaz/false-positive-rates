rm(list = ls()) # clear all vars from the current workspace

simTrials <- function (
    monteCarloSims = 10000L,
    participantsPerArm = 100L,
    baselineRisk = 0.5,
    absoluteRR = 0.1) {
    
    ctrlOutcomes <- rbinom(monteCarloSims, participantsPerArm, baselineRisk)
    
    df <- monteCarloFisher(n = participantsPerArm,
                           outcomes = cbind(
                               rbinom(monteCarloSims, participantsPerArm, baselineRisk), # no dif interv
                               ctrlOutcomes, # control group
                               rbinom(monteCarloSims, participantsPerArm, baselineRisk - absoluteRR), # working Rx
                               ctrlOutcomes))
                                            
    
    sub <- "no difference"
    hist(df$p.1, sub=sub)
    hist(df$p.1[df$p.1 < 0.1], sub=sub)
    
    sub <- paste0("control:", baselineRisk, "; intervention:", baselineRisk - absoluteRR)
    hist(df$p.2, sub=sub)
    hist(df$p.2[df$p.2 < 0.1], sub=sub)
    
    false_rx_works <- sum(df$or.1 < 1 & df$p.1 < 0.05)
    false_rx_harmfull <- sum(df$or.1 > 1 & df$p.1 < 0.05)
    highly_false_pos <- sum(df$p.1 < 0.02)
    
    false_negative <- sum(df$or.2 > 1 | df$p.2 >= 0.05)
    false_rx_type3 <- sum(df$or.2 > 1 & df$p.2 < 0.05)
    
    resultOut <- function(text, count){
        return(paste0("\n\t", text, ": ", count, " (", 100 * count / monteCarloSims, "%)"))
    }
    cat("NO DIFFERENCE:", 
        resultOut("treatment shown to work", false_rx_works),
        resultOut("treatment shown to be harmful", false_rx_harmfull),
        resultOut("false positive rate (p < 0.05)", false_rx_works + false_rx_harmfull),
        resultOut("includes highly false positive (p < 0.02)", highly_false_pos),
        "\nDIFFERENCE (", sub, "):",
        resultOut("false negative rate ", false_negative),
        resultOut("including opposite finding of harm (p < 0.05)", false_rx_type3))
}

simTrials(
    monteCarloSims = 1000000L,
    participantsPerArm = 100L,
    baselineRisk = 0.5,
    absoluteRR = 0.1)
