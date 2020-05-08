rm(list = ls()) # clear all vars from the current workspace

do_trial <- function(p_control=0.5, p_treatment=0.8, num_subjects=100){
    binom_control = rbinom(1,num_subjects,p_control)
    binom_treatment = rbinom(1,num_subjects,p_treatment)
    rbind(control=c(live=binom_control, 
                    die=(num_subjects - binom_control)),
        treatment=c(live=binom_treatment, 
                    die=(num_subjects - binom_treatment)))
}

chisq_parts <- function(trial){
    chi <- chisq.test(trial, correct = TRUE)
    squared <- chi$statistic[["X-squared"]]
    p <- chi$p.value
    result <- c(p_value=p, x_squared=squared)
    return(result)
}

is_significant <- function(p_value){
    return(p_value >= 0.03 && p_value <= 0.05)
}

calculate_phi <- function(x_squared, num_subjects){
    sqrt((x_squared * x_squared)/num_subjects)
}

analyze_trials <- function(trial_list){
    chi = Map(chisq_parts, trial_list)
    parts = c('p_value', 'x_squared')
    df <- data.frame(do.call(rbind, lapply(chi, `[`, parts)))

    calculate_phi <- function(x_squared){
        sqrt((x_squared * x_squared)/num_subjects)
    }

    df$phi <- do.call(calculate_phi, df['x_squared'])
    df$significant <- do.call(is_significant, df['p_value'])
    return(df)
}

display_comparison <- function(control_statistics, treatment_statistics){
    control_significant = subset(control_statistics, significant==TRUE)
    treatment_significant = subset(treatment_statistics, significant==TRUE)
    num_false_positives = nrow(control_significant)
    num_true_positives = nrow(num_true_positives)
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

run_simulation <- function(num_trials=3,
                           num_subjects=100,
                           p_control=0.5,
                           p_treatment=0.8){
    control_trials = run_trials(num_trials, num_subjects, p_control, p_control)
    treatment_trials = run_trials(num_trials, num_subjects, p_control, p_treatment)

    control_statistics = analyze_trials(control_trials)
    treatment_statistics = analyze_trials(treatment_trials)

    display_comparison(control_statistics, treatment_statistics)
}

num_trials=10000
num_subjects=100
p_control=0.5
p_treatment=0.8

run_simulation(num_trials,num_subjects,p_control,p_treatment)




