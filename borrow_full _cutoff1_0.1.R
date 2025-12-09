# Clean R workspace
rm(list=ls())

# Load required packages
library(tidyr)
library(dplyr)
library(startupmsg)
library(distr)
library(SAMprior)
library(rmutil)
library(ggplot2)
library(RBesT)

# Bayesian decision function setup
decision <- decision2S(0.9, 0, lower.tail=FALSE) # Pr(dORR > 0) > 90%
decision3 <- decision2S(0.9, 0.05, lower.tail=FALSE) # Pr(dORR > 0.05) > 90%

#----------------------------------------------
# Posterior calculation with dose selection first
#----------------------------------------------

# sample_b: Binary outcome vector for a dropped or borrowed treatment arm (e.g., dose level 1); used to estimate historical or borrowed response rate incorporated as prior information.
# sample_c: Binary outcome vector for the control arm; represents observed responses in the control group for analysis and model fitting.
# sample_t: Binary outcome vector for the treatment arm under current evaluation (e.g., dose level 2); represents observed responses for the dose currently compared.
# delta: Clinically significant effect size or minimum clinically meaningful difference in response rates; used as a key parameter in Bayesian prior construction and dose selection criteria.
p.f <- function(sample_b, sample_c, sample_t, delta) {
  # Calculate ORR and sample sizes
  orr_low <- mean(sample_b)
  n_b <- length(sample_b)
  orr_c <- mean(sample_c)
  ss_trt <- length(sample_t)
  ss_control <- length(sample_c)
  
  # Set informative priors
  nf.prior <- mixbeta(nf.prior = c(1, 0.5, 0.5))
  post_c   <- postmix(nf.prior, n = ss_control, r = sum(sample_c))
  prior.hist <- mixbeta(c(1, n_b * orr_low, n_b * (1 - orr_low)))
  
  post_t <- postmix(prior.hist, n = ss_trt, r = sum(sample_t))
  
  # Robust Mixture Prior calculations
  w <- 0.3
  orr.trt <- sum(sample_t) / ss_trt
  Lf1 <- dbeta(orr.trt, prior.hist[2, 1], prior.hist[3, 1])
  Lfv <- dbeta(orr.trt, 0.5, 0.5)
  w.new <- w * Lf1 / (w * Lf1 + (1 - w) * Lfv)
  
  RMP.prior <- SAM_prior(
    if.prior = prior.hist,
    nf.prior = nf.prior,
    weight = w.new,
    delta = delta
  )
  RMP.post <- postmix(RMP.prior, n = sum(ss_trt), r = sum(sample_t))
  p.DRMP   <- decision(RMP.post, post_c)
  drmp_w_cum <- ess(RMP.post, method = "moment")
  
  p.t <- decision(post_t, post_c)
  
  # Return selected outputs
  c(p.DRMP, p.t, drmp_w_cum, w.new)
}

#----------------------------------------------
# Simulation function setup for study designs
#----------------------------------------------

# ss_control_IA: Planned sample size for the control arm at interim analysis (IA); number of control subjects evaluated in the interim stage.
# ss_trt_IA: Planned sample size for the treatment arm at interim analysis (IA); number of treatment subjects evaluated in the interim stage.
# ss_control_FA: Planned sample size for the control arm at final analysis (FA); total control subjects for the full study duration.
# ss_trt_FA: Planned sample size for the treatment arm at final analysis (FA); total treatment subjects for the full study duration.
# orr_trt1: True or expected response rate (ORR) for treatment dose level 1; probability of treatment response in simulation for dose 1.
# orr_trt2: True or expected response rate (ORR) for treatment dose level 2; probability of treatment response in simulation for dose 2.
# orr_control: True or expected response rate (ORR) for the control arm; probability of response in the control group used for simulation.
simFunc <- function(
    ss_control_IA, ss_trt_IA, ss_control_FA, ss_trt_FA,
    orr_trt1, orr_trt2, orr_control
) {
  # Parameters & containers
  scan2_t <- 3         # Months to confirm response
  acc_rate <- 5        # Accrual rate (per month)
  delta <- 0.15        # Clinically significant difference
  nk <- 10000          # Simulation iterations
  cutoff1 <- 0.1       # Dose-selection cutoff
  
  # Results containers
  sel.lower1 <- sel.higher1 <- sel.lower2 <- sel.higher2 <- nogo <- 0
  rmp_w_cum <- rmp_w <- test.freq <- res_Non_tmp <- res_RMP_tmp <- numeric(0)
  res_Non_high <- res_Non_low <- res_RMP_high <- res_RMP_low <- numeric(0)
  
  # Start simulations
  for (i in 1:nk) {
    print(i)
    tryCatch({
      set.seed(i + 2024)
      
      # Accrual times for IA samples
      acc_t0    <- runif(ss_control_IA, 0, ss_control_IA / acc_rate)
      acc_t1.dl1 <- runif(ss_trt_IA, 0, ss_trt_IA / acc_rate)
      acc_t1.dl2 <- runif(ss_trt_IA, 0, ss_trt_IA / acc_rate)
      acc.max <- max(acc_t1.dl1, acc_t1.dl2, acc_t0)
      
      # Observed patients for analysis after scan period
      ss_control <- acc_t0       <= (acc.max - scan2_t)
      ss_trt1   <- acc_t1.dl1    <= (acc.max - scan2_t)
      ss_trt2   <- acc_t1.dl2    <= (acc.max - scan2_t)
      
      # Generate efficacy samples
      sample_c_FA <- rbinom(ss_control_FA, 1, orr_control)
      sample_c_IA <- sample_c_FA[1:ss_control_IA]
      sample_c    <- sample_c_IA[ss_control]
      
      sample_t1_FA <- rbinom(ss_trt_FA, 1, orr_trt1)
      sample_t1_IA <- sample_t1_FA[1:ss_trt_IA]
      sample_t1    <- sample_t1_IA[ss_trt1]
      
      sample_t2_FA <- rbinom(ss_trt_FA, 1, orr_trt2)
      sample_t2_IA <- sample_t2_FA[1:ss_trt_IA]
      sample_t2    <- sample_t2_IA[ss_trt2]
      
      # Calculate differences & cutoffs
      test.cutoff1 <- mean(sample_t1)    - mean(sample_c)    # prior to IA: calculate Delta_1 = ORR1 - ORRc
      test.cutoff2 <- mean(sample_t2)    - mean(sample_c)    # prior to IA: calculate Delta_2 = ORR2 - ORRc
      test.cutoff3 <- mean(sample_t1_IA) - mean(sample_c_IA) # Interim Analysis: ORR1 - ORRc
      test.cutoff4 <- mean(sample_t2_IA) - mean(sample_c_IA) # Interim Analysis: ORR2 - ORRc
      
      # Define non-informative & informative prior values
      nf.prior <- mixbeta(nf.prior = c(1, 0.5, 0.5))
      post_t1 <- postmix(nf.prior, n = sum(ss_trt1), r = sum(sample_t1))
      post_t2 <- postmix(nf.prior, n = sum(ss_trt2), r = sum(sample_t2))
      
      p.diff  <- decision3(post_t2, post_t1)
      p.diff2 <- decision3(post_t1, post_t2)
      post_t3 <- postmix(nf.prior, n = length(sample_t1_IA), r = sum(sample_t1_IA))
      post_t4 <- postmix(nf.prior, n = length(sample_t2_IA), r = sum(sample_t2_IA))
      p.diff3 <- decision3(post_t4, post_t3)
      p.diff4 <- decision3(post_t3, post_t4)
      
      # Dose selection decision tree (IA1 and IA2)
      # See logic comment blocks
      p.res <- c(0, 0, 0, 0)
      ind <- 0
      
      # IA1 decision logic
      if (min(test.cutoff1, test.cutoff2) > cutoff1 & p.diff == 1) {
        ind <- 1
        p.res <- p.f(sample_t1_IA, sample_c_FA, sample_t2_FA, delta=delta)
        sel.higher1 <- sel.higher1 + 1
        res_Non_high <- c(res_Non_high, p.res[2])
        res_RMP_high <- c(res_RMP_high, p.res[1])
      } else if (min(test.cutoff1, test.cutoff2) > cutoff1 & p.diff2 == 1) {
        ind <- 2
        p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
        sel.lower1 <- sel.lower1 + 1
        res_Non_low <- c(res_Non_low, p.res[2])
        res_RMP_low <- c(res_RMP_low, p.res[1])
      } else if (min(test.cutoff1, test.cutoff2) > cutoff1) {
        ind <- 2
        p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
        sel.lower1 <- sel.lower1 + 1
        res_Non_low <- c(res_Non_low, p.res[2])
        res_RMP_low <- c(res_RMP_low, p.res[1])
      } else if (test.cutoff1 > cutoff1) {
        ind <- 2
        p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
        sel.lower1 <- sel.lower1 + 1
        res_Non_low <- c(res_Non_low, p.res[2])
        res_RMP_low <- c(res_RMP_low, p.res[1])
      } else if (test.cutoff2 > cutoff1) {
        ind <- 1
        p.res <- p.f(sample_t1_IA, sample_c_FA, sample_t2_FA, delta=delta)
        sel.higher1 <- sel.higher1 + 1
        res_Non_high <- c(res_Non_high, p.res[2])
        res_RMP_high <- c(res_RMP_high, p.res[1])
      } else {
        # IA2 decision logic
        if (min(test.cutoff3, test.cutoff4) > cutoff1 & p.diff3 == 1) {
          ind <- 1
          p.res <- p.f(sample_t1_IA, sample_c_FA, sample_t2_FA, delta=delta)
          sel.higher2 <- sel.higher2 + 1
          res_Non_high <- c(res_Non_high, p.res[2])
          res_RMP_high <- c(res_RMP_high, p.res[1])
        } else if (min(test.cutoff3, test.cutoff4) > cutoff1 & p.diff4 == 1) {
          ind <- 2
          p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
          sel.lower2 <- sel.lower2 + 1
          res_Non_low <- c(res_Non_low, p.res[2])
          res_RMP_low <- c(res_RMP_low, p.res[1])
        } else if (min(test.cutoff3, test.cutoff4) > cutoff1) {
          ind <- 2
          p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
          sel.lower2 <- sel.lower2 + 1
          res_Non_low <- c(res_Non_low, p.res[2])
          res_RMP_low <- c(res_RMP_low, p.res[1])
        } else if (test.cutoff3 > cutoff1) {
          ind <- 2
          p.res <- p.f(sample_t2_IA, sample_c_FA, sample_t1_FA, delta=delta)
          sel.lower2 <- sel.lower2 + 1
          res_Non_low <- c(res_Non_low, p.res[2])
          res_RMP_low <- c(res_RMP_low, p.res[1])
        } else if (test.cutoff4 > cutoff1) {
          ind <- 1
          p.res <- p.f(sample_t1_IA, sample_c_FA, sample_t2_FA, delta=delta)
          sel.higher2 <- sel.higher2 + 1
          res_Non_high <- c(res_Non_high, p.res[2])
          res_RMP_high <- c(res_RMP_high, p.res[1])
        } else {
          ind <- 0
          nogo <- nogo + 1
        }
      }
      
      # Frequentist test: t-test (unpaired, one-sided)
      if (ind == 1) {
        test.freq[i] <- t.test(sample_t2_FA, sample_c_FA,
                               paired=FALSE, alternative="greater")$p.value
      } else if (ind == 2) {
        test.freq[i] <- t.test(sample_t1_FA, sample_c_FA,
                               paired=FALSE, alternative="greater")$p.value
      } else {
        test.freq[i] <- 1
      }
      
      # Collect simulation outputs
      rmp_w_cum <- c(rmp_w_cum, p.res[3])
      rmp_w <- c(rmp_w, p.res[4])
      res_RMP_tmp <- c(res_RMP_tmp, p.res[1])
      res_Non_tmp <- c(res_Non_tmp, p.res[2])
      
    }, error = function(e) {
      cat("Simulation", i, "encountered an error:", conditionMessage(e), "\n")
      NULL # Skip on error
    })
  }
  
  # Update nk (if errors skipped)
  nk <- length(res_Non_tmp)
  out <- c(
    prop.table(table(res_Non_tmp))[2],
    table(res_Non_low)[2] / nk,
    table(res_Non_high)[2] / nk,
    prop.table(table(res_Non_tmp))[1]
  ) * 100
  out <- t(data.frame(out))
  
  # Dose selection percentages: IA1, IA2, Total
  dosepick1 <- 100 * c(sel.lower1/nk, sel.higher1/nk, 1-(sel.lower1+sel.higher1)/nk)
  dosepick2 <- 100 * c(sel.lower2/nk, sel.higher2/nk,
                       (1-(sel.lower1+sel.higher1)/nk)-(sel.lower2+sel.higher2)/nk)
  dosepicktot <- 100 * c(
    sel.lower1/nk+sel.lower2/nk, sel.higher1/nk+sel.higher2/nk,
    1-(sel.lower1/nk+sel.lower2/nk + sel.higher1/nk+sel.higher2/nk)
  )
  
  dosepick <- t(data.frame(c(dosepick1, dosepick2, dosepicktot)))
  
  # Final output
  return(list(out = out, dosepick = dosepick))
}


#----------------------------------------------
# Example scenario using 
#  (1) IA n = 35 and FA n = 60
#  (2) S2: (theta_L, theta_H, theta_C) = (0.55, 0.75, 0.55)
#----------------------------------------------

t1 <- Sys.time()
S2 = simFunc(
  ss_control_IA=35, ss_trt_IA=35, ss_control_FA=60, ss_trt_FA=60, 
  orr_trt1=0.55, orr_trt2=0.75, orr_control=0.55
  )
t2 <- Sys.time()
t2 - t1
# Time difference of 2.914868 mins

rownames_trueORR = "(0.55, 0.75, 0.55)"
colnames_poc = c("Go", "DL1 Go", "DL2 Go", "No Go")
colnames_do = c("IA1:DL1", "IA1:DL2", "IA1:NoGo", "IA2:DL1", "IA2:DL2", "IA2:NoGo", "Total:DL1", "Total:DL2", "Total:NoGo")

poc_full = data.frame(S2$out)
rownames(poc_full) = rownames_trueORR
colnames(poc_full) = colnames_poc

dosepick_full = data.frame(S2$dosepick)
rownames(dosepick_full) = rownames_trueORR
colnames(dosepick_full) = colnames_do

# save(poc_full, file = "poc_full.RData")
# save(dosepick_full, file = "dosepick_full.RData")









