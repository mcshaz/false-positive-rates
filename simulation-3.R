rm(list = ls()) # clear all vars from the current workspace

do_trial <- function(p_control=0.5, p_treatment=0.8, num_subjects=100){
    binom_control = rbinom(1,num_subjects,p_control)
    binom_treatment = rbinom(1,num_subjects,p_treatment)
    rbind(control=c(live=binom_control, 
                    die=(num_subjects - binom_control)),
        treatment=c(live=binom_treatment, 
                    die=(num_subjects - binom_treatment)))
}


# do_trial <- function(p_control=0.5, p_treatment=0.8, num_subjects=100){
#     binom_control = rbinom(1,num_subjects,p_control)
#     binom_treatment = rbinom(1,num_subjects,p_treatment)
#     rbind(live=c(control=binom_control, 
#                 treatment=binom_treatment),
#           die=c(control=(num_subjects - binom_control), 
#                 treatment=(num_subjects - binom_treatment)))
# }

chisq_parts <- function(trial){
    chi <- chisq.test(trial, correct = TRUE)
    squared <- chi$statistic[["X-squared"]]
    p <- chi$p.value
    result <- c(p_value=p, x_squared=squared)
    return(result)
}

matrix_to_p_value <- function(m){
    chi = chisq.test(m)
    return(chi$p.value)
}

matrix_to_x_squared <- function(m){
    chi = chisq.test(m)
    return(chi$statistic[["X-squared"]])
}

is_significant <- function(p_value){
    return(p_value >= 0.03 && p_value <= 0.05)
}

calculate_phi <- function(x_squared, num_subjects){
    sqrt((x_squared * x_squared)/num_subjects)
}

analyze_trials <- function(trial_list){
    chi <- lapply(trial_list, chisq_parts)
    df <- data.frame(do.call(rbind, chi))
    calculate_phi <- function(x_squared){
        sqrt((x_squared * x_squared)/num_subjects)
    }
    df$phi <- sapply(df$x_squared, calculate_phi)
    df$significant <- sapply(df$p_value, is_significant)
    return(df)
}

display_comparison <- function(control_statistics, treatment_statistics){
    control_significant = subset(control_statistics, significant==TRUE)
    treatment_significant = subset(treatment_statistics, significant==TRUE)
    num_false_positives = nrow(control_significant)
    num_true_positives = nrow(treatment_significant)
    false_positive_risk = num_false_positives / (num_true_positives + num_false_positives)
    
    cat( "\n",
        "False Positives: ", num_false_positives, "\n",
        "True Positives: ", num_true_positives, "\n",
        "False Positive Risk: ", false_positive_risk, "\n")
    #TODO: graph the differences
}

run_trials <- function(num_trials,
                       num_subjects,
                       p_control, 
                       p_treatment){
    replicate(num_trials, 
              do_trial(p_control=p_control, 
                       p_treatment=p_treatment,
                       num_subjects=num_subjects), 
              simplify = FALSE)
}

num_trials=100000
num_subjects=100
p_control=0.5
p_treatment=0.8

control_trials = run_trials(num_trials, num_subjects, p_control, p_control)
treatment_trials = run_trials(num_trials, num_subjects, p_control, p_treatment)

control_statistics = analyze_trials(control_trials)
treatment_statistics = analyze_trials(treatment_trials)

display_comparison(control_statistics, treatment_statistics)

hist(control_statistics$p_value[control_statistics$significant == TRUE])
hist(treatment_statistics$p_value[treatment_statistics$significant == TRUE])
