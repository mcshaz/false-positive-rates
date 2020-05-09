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
    return(p_value >= 0.03 & p_value <= 0.05)
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

    x_squared_mean <- mean(control_statistics$x_squared[control_statistics$p_value <= 0.05])
    x_squared_sd <- sd(control_statistics$x_squared[control_statistics$p_value <= 0.05])

    phi_mean <- mean(control_statistics$phi[control_statistics$p_value <= 0.05])
    phi_sd <- sd(control_statistics$phi[control_statistics$p_value <= 0.05])

    x_squared_mean <- mean(treatment_statistics$x_squared[treatment_statistics$p_value <= 0.05])
    x_squared_sd <- sd(treatment_statistics$x_squared[treatment_statistics$p_value <= 0.05])

    phi_mean <- mean(treatment_statistics$phi[treatment_statistics$p_value <= 0.05])
    phi_sd <- sd(treatment_statistics$phi[treatment_statistics$p_value <= 0.05])

    cat( "\nControl vs Control\n",
    "    x_squared mean,", x_squared_mean, "\n",
    "    x_squared sd:", x_squared_sd, "\n",
    "    phi mean,", phi_mean, "\n",
    "    phi sd:", phi_sd, "\n",
    "\nTreatment vs Control\n",
    "    x_squared mean,", x_squared_mean, "\n",
    "    x_squared sd:", x_squared_sd, "\n",
    "    phi mean,", phi_mean, "\n",
    "    phi sd:", phi_sd, "\n")
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

# Display p_values for significant control vs control trials:
hist(control_statistics$p_value[control_statistics$significant == TRUE])

# Display p_values for significant treatment vs control trials:
hist(treatment_statistics$p_value[treatment_statistics$significant == TRUE])

# Example histogram of phi specifying a p_value range
# hist(treatment_statistics$phi[treatment_statistics$p_value >= 0.03 & treatment_statistics$p_value <= 0.05])

