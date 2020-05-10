rm(list = ls()) # clear all vars from the current workspace

# Generates a 2x2 matrix binomial distribution live vs die for treatment and control
do_trial <- function(p_control=0.5, p_treatment=0.8, num_subjects=100){
    binom_control = rbinom(1,num_subjects,p_control)
    binom_treatment = rbinom(1,num_subjects,p_treatment)
    rbind(control=c(live=binom_control, 
                    die=(num_subjects - binom_control)),
        treatment=c(live=binom_treatment, 
                    die=(num_subjects - binom_treatment)))
}

# Generates a list of trial matrices
run_trials <- function(num_trials=1000,
                       num_subjects=100,
                       p_control=0.5, 
                       p_treatment=0.8){
    replicate(num_trials, 
              do_trial(p_control=p_control, 
                       p_treatment=p_treatment,
                       num_subjects=num_subjects), 
              simplify = FALSE)
}

# Calculates the chi square of a trial and returns the p and x^2
matrix_to_chisq_parts <- function(trial){
    chi <- chisq.test(trial, correct = TRUE)
    squared <- chi$statistic[["X-squared"]]
    p <- chi$p.value
    result <- c(p_value=p, x_squared=squared)
    return(result)
}

# Calculates the chi square of a trial and returns the p
matrix_to_p_value <- function(m){
    chi = chisq.test(m)
    return(chi$p.value)
}

# Calculates the chi square of a trial and returns the x^2
matrix_to_x_squared <- function(m){
    tryCatch(
        {
            chi = chisq.test(m)
            return(chi$statistic[["X-squared"]])
        },
        error=function(cond) {
            message("Unable to calculate chi square for " + m)
            message("Here's the original error message:")
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            message("Matrix caused a warning during chi square: " + m)
            message("Here's the original warning message:")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={}
    )    
}

# Returns TRUE if a p value is significant
p_value_is_significant <- function(p_value){
    return(p_value <= 0.05)
    #return(p_value >= 0.03 & p_value <= 0.05)
}

# Converts a list of trial matrices into a dataframe with columns for
# p_value, x_squared, phi, significant
analyze_trials <- function(trial_list){
    chi <- lapply(trial_list, matrix_to_chisq_parts)
    df <- data.frame(do.call(rbind, chi))
    calculate_phi <- function(x_squared){
        sqrt((x_squared * x_squared)/(2*num_subjects))
    }
    df$phi <- sapply(df$x_squared, calculate_phi)
    df$significant <- sapply(df$p_value, p_value_is_significant)
    return(df)
}

# Calculate and display the false positive risk for a control vs control and
# a control vs treatment dataframe
display_false_positive_risk <- function(control_statistics, treatment_statistics){
    control_significant = subset(control_statistics, significant==TRUE)
    treatment_significant = subset(treatment_statistics, significant==TRUE)
    num_false_positives = nrow(control_significant)
    num_true_positives = nrow(treatment_significant)
    false_positive_risk = num_false_positives / (num_true_positives + num_false_positives)
    
    cat( "\n",
        "False Positives: ", num_false_positives, "\n",
        "True Positives: ", num_true_positives, "\n",
        "False Positive Risk: ", false_positive_risk, "\n\n")
}

# Calculate and display phi and chi averages for the control and treatment dataframes
display_phi_and_chi_averages = function(control_statistics, treatment_statistics){
    control_x_squared_mean <- mean(control_statistics$x_squared[control_statistics$p_value <= 0.05])
    control_x_squared_sd <- sd(control_statistics$x_squared[control_statistics$p_value <= 0.05])
    control_phi_mean <- mean(control_statistics$phi[control_statistics$p_value <= 0.05])
    control_phi_sd <- sd(control_statistics$phi[control_statistics$p_value <= 0.05])

    treatment_x_squared_mean <- mean(treatment_statistics$x_squared[treatment_statistics$p_value <= 0.05])
    treatment_x_squared_sd <- sd(treatment_statistics$x_squared[treatment_statistics$p_value <= 0.05])
    treatment_phi_mean <- mean(treatment_statistics$phi[treatment_statistics$p_value <= 0.05])
    treatment_phi_sd <- sd(treatment_statistics$phi[treatment_statistics$p_value <= 0.05])

    cat( "\nControl vs Control\n",
    "    x_squared mean,", control_x_squared_mean, "\n",
    "    x_squared sd:", control_x_squared_sd, "\n",
    "    phi mean,", control_phi_mean, "\n",
    "    phi sd:", control_phi_sd, "\n",
    "\nTreatment vs Control\n",
    "    x_squared mean,", treatment_x_squared_mean, "\n",
    "    x_squared sd:", treatment_x_squared_sd, "\n",
    "    phi mean,", treatment_phi_mean, "\n",
    "    phi sd:", treatment_phi_sd, "\n\n")
}

bonus_code <- function(control_statistics, treatment_statistics){

    # Display p_values for significant control vs control trials:
    hist(control_statistics$p_value[control_statistics$significant == TRUE])

    # Display p_values for significant treatment vs control trials:
    hist(treatment_statistics$p_value[treatment_statistics$significant == TRUE])

    # num control vs control trials where p <= 0.05
    length(control_statistics$p_value[control_statistics$p_value <= 0.05])

    # num treatment vs control trials where p <= 0.05
    length(treatment_statistics$p_value[treatment_statistics$p_value <= 0.05])

    # Example histogram of phi specifying a p_value range
    hist(treatment_statistics$phi[treatment_statistics$p_value >= 0.03 & treatment_statistics$p_value <= 0.05])
}


num_trials=100000
num_subjects=100
p_control=0.5
p_treatment=0.8

control_trials = run_trials(num_trials, num_subjects, p_control, p_control)
treatment_trials = run_trials(num_trials, num_subjects, p_control, p_treatment)

control_statistics = analyze_trials(control_trials)
treatment_statistics = analyze_trials(treatment_trials)

display_false_positive_risk(control_statistics, treatment_statistics)
display_phi_and_chi_averages(control_statistics, treatment_statistics)

