# Clear the environment
# rm(list = ls())

num_sims=3
sim_size=100
p_placebo=0.5
p_treatment=0.8

matrices = replicate(num_sims, 
run_sim(sim_size=sim_size, 
 p_placebo=p_placebo, 
 p_treatment=p_treatment), 
simplify = FALSE)


data.frame(simulation=matrices)

t <- data.frame(id=seq(1,3,1))
t <- cbind(t, simulation=matrices)


# sims <- do.call(rbind.data.frame, matrices)
# sims <- cbind(simulation=matrices)
# sims <- cbind(sims, chisq=Map(to_chisq, sims[,"simulation"]))
# sims <- cbind(sims, phi=Map(to_phi, sims[,"chisq"]))
# sims <- cbind(sims, significant=Map(is_significant, sims[,"chisq"]))