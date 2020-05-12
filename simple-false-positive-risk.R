rm(list = ls()) # clear all vars from the current workspace

do_test <- function(placebo=0.5, treatment=0.5, test=chisq.test){
    test_1 <- rbinom(1,100,placebo)
    test_1_minus100 <- 100 - test_1
    dead_alive_1 <- c(test_1, test_1_minus100)
    test_2 <- rbinom(1,100,treatment)
    test_2_minus100 <- 100 - test_2
    dead_alive_2 <- c(test_2, test_2_minus100)
    matrix_dead_alive <- rbind(dead_alive_1,dead_alive_2)
    cs <- test(matrix_dead_alive)
    return(cs$p.value)
}

is_significant <- function(p_value){
    return(0.03 <= p_value && p_value <= 0.05)
}

placebo_prob = 0.5
treatment_prob = 0.6
test = fisher.test ## chisq.test

placebo_vs_placebo <- replicate(100000, do_test(placebo_prob,placebo_prob, test))
treatment_vs_placebo <- replicate(100000, do_test(placebo_prob,treatment_prob, test))

false_positives = length(Filter(is_significant, placebo_vs_placebo))
true_positives = length(Filter(is_significant, treatment_vs_placebo))

cat("False Positives: ", false_positives, "\n",
    "True Positives: ", true_positives, "\n",
    "False Positive Risk: ", false_positives / (true_positives + false_positives), "\n")

subtitle <- ifelse(identical(test, fisher.test), "Fisher's Exact", "Chi-squared")

hist(placebo_vs_placebo, sub=subtitle)
hist(treatment_vs_placebo, sub=subtitle)
