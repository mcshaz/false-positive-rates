# Clear the environment
# rm(list = ls())

num_trials=3
num_subjects=100
p_control=0.5
p_treatment=0.8


do_trial <- function(p_control=0.5, p_treatment=0.5, num_subjects=100){
    binom_control = rbinom(1,num_subjects,p_control)
    binom_treatment = rbinom(1,num_subjects,p_treatment)
    rbind(control=c(live=binom_control, 
                    die=(num_subjects - binom_control)),
        treatment=c(live=binom_treatment, 
                    die=(num_subjects - binom_treatment)))
}


# do_trial_2 <- function(p_control=0.5, p_treatment=0.5, num_subjects=100){
#     binom_control = rbinom(1,num_subjects,p_control)
#     binom_treatment = rbinom(1,num_subjects,p_treatment)
#     rbind(live=c(control=binom_control, 
#                  treatment=binom_treatment),
#           die=c(control=(num_subjects - binom_control), 
#                 treatment=(num_subjects - binom_treatment)))
# }

matrices = replicate(num_trials, 
                     do_trial(p_control=p_control, 
                             p_treatment=p_treatment,
                             num_subjects=num_subjects), 
                     simplify = FALSE)

# matrices = replicate(num_trials, 
#                      do_trial_2(p_control=p_control, 
#                              p_treatment=p_treatment,
#                              num_subjects=num_subjects), 
#                      simplify = FALSE)

data.frame(simulation=matrices)

t <- data.frame(id=seq(1,3,1))
t <- cbind(t, simulation=matrices)


# sims <- do.call(rbind.data.frame, matrices)
# sims <- cbind(simulation=matrices)
# sims <- cbind(sims, chisq=Map(to_chisq, sims[,"simulation"]))
# sims <- cbind(sims, phi=Map(to_phi, sims[,"chisq"]))
# sims <- cbind(sims, significant=Map(is_significant, sims[,"chisq"]))