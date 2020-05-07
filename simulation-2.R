num_iterations = 1000
placebo_prob = 0.5
treatment_prob = 0.6

binom_matrix <- function(placebo=0.5, treatment=0.5){
    test_1 <- rbinom(1,100,placebo)
    test_1_minus100 <- 100 - test_1
    dead_alive_1 <- c(test_1, test_1_minus100)
    test_2 <- rbinom(1,100,treatment)
    test_2_minus100 <- 100 - test_2
    dead_alive_2 <- c(test_2, test_2_minus100)
    matrix_dead_alive <- rbind(dead_alive_1,dead_alive_2)
    return(matrix_dead_alive)
}

to_p_value <- function(m){
    return(chisq.test(m)$p.value)
}

is_significant <- function(p_value){
    return(p_value >= 0.03 && p_value <= 0.05)
}

is_significant_matrix <- function(m){
    return(is_significant(to_p_value(m)))
}

placebo_vs_placebo <- replicate(num_iterations, binom_matrix(placebo_prob,placebo_prob), simplify = FALSE)
treatment_vs_placebo <- replicate(num_iterations, binom_matrix(placebo_prob,treatment_prob), simplify = FALSE)

false_positives = Filter(is_significant_matrix, placebo_vs_placebo)
true_positives = Filter(is_significant_matrix, treatment_vs_placebo)

num_false_positives = length(false_positives)
num_true_positives = length(true_positives)
false_positive_risk = num_false_positives / (num_true_positives + num_false_positives)

cat("False Positives: ", num_false_positives, "\n",
    "True Positives: ", num_true_positives, "\n",
    "False Positive Risk: ", false_positive_risk)

