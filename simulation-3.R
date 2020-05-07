num_iterations = 1000
placebo_prob = 0.5
treatment_prob = 0.6
test_size = 100

#install.packages("functional")
#library(functional)

run_sim <- function(p_placebo=0.5, p_treatment=0.5, sim_size=100){
    test_1 <- rbinom(1,sim_size,p_placebo)
    test_1_minus100 <- sim_size - test_1
    dead_alive_1 <- c(test_1, test_1_minus100)
    test_2 <- rbinom(1,sim_size,p_treatment)
    test_2_minus100 <- sim_size - test_2
    dead_alive_2 <- c(test_2, test_2_minus100)
    matrix_dead_alive <- rbind(dead_alive_1,dead_alive_2)
    return(matrix_dead_alive)
}

to_chisq <- function(m){
    return(chisq.test(m, correct = TRUE))
}

is_significant <- function(chi){
    return(chi$p.value >= 0.03 && chi$p.value <= 0.05)
}

is_significant_matrix <- function(m){
    return(is_significant(to_chisq(m)))
}

to_phi <- function(chi){
    chi$statistic
}

generate_table <- function(num_sims=3, 
                           sim_size=100, 
                           p_placebo=0.5, 
                           p_treatment=0.8){
    matrices = replicate(num_sims, 
                         run_sim(sim_size=sim_size, 
                                 p_placebo=p_placebo, 
                                 p_treatment=p_treatment), 
                         simplify = FALSE)
    sims <- cbind(simulation=matrices)
    sims <- cbind(sims, chisq=Map(to_chisq, sims[,"simulation"]))
    sims <- cbind(sims, phi=Map(to_phi, sims[,"chisq"]))
    sims <- cbind(sims, significant=Map(is_significant, sims[,"chisq"]))
    return(sims)
}

placebo_vs_placebo = generate_table(p_treatment=placebo_prob)
treatment_vs_placebo = generate_table(p_treatment=treatment_prob)



# sims[,"simulation"]


placebo_vs_placebo <- replicate(num_iterations, binom_matrix(placebo_prob,placebo_prob), simplify = FALSE)
treatment_vs_placebo <- replicate(num_iterations, binom_matrix(placebo_prob,treatment_prob), simplify = FALSE)

placebo_vs_placebo_chi <- Map(to_chisq, placebo_vs_placebo)
treatment_vs_placebo_chi <- Map(to_chisq, treatment_vs_placebo)

false_positives = Filter(is_significant, placebo_vs_placebo_chi)
true_positives = Filter(is_significant, treatment_vs_placebo_chi)

num_false_positives = length(false_positives)
num_true_positives = length(true_positives)
false_positive_risk = num_false_positives / (num_true_positives + num_false_positives)

cat("False Positives: ", num_false_positives, "\n",
    "True Positives: ", num_true_positives, "\n",
    "False Positive Risk: ", false_positive_risk)




as.vector(replicate(3, binom_matrix(placebo_prob,placebo_prob), simplify = FALSE))


