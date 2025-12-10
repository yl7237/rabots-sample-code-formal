##Traditional fixed three arm design
##Date: 8Dec2025
##Use information borrowing at the FA and select dose at the IA in dose optimization

rm(list=ls())

library(tidyr)
library(dplyr)
library(startupmsg)
library(distr)
library(rmutil)
library(ggplot2)
library(RBesT)
library(SAMprior)
library(tidyverse)
library(stringr)
library(patchwork)

#This is when to calculate dual criteria

#decision2S: Bayesian probability decision function

#overall compared with control
decision = decision2S(0.9, 0, lower.tail=FALSE) #Pr(dORR>0)>90%
# decision.1 = decision2S(0.5, 0.1, lower.tail=FALSE) #Pr(dORR>0.1)>50%
decision3 = decision2S(0.9, 0.05, lower.tail=FALSE) #Pr(dORR>0.05)>90% decision tree when two ORR are close

p.f <- function(sample_c, sample_t, delta) {
  
  orr_c <- mean(sample_c) # ORR for control
  nf.prior <- mixbeta(nf.prior = c(1,0.5,0.5)) # full-informative prior
  
  ss_trt <- length(sample_t)
  ss_control <- length(sample_c)
  ## Assume beta(0.5,0.5) as the full-informative prior used for mixture
  post_c <- postmix(nf.prior, n = ss_control, r = sum(sample_c)) #control: full-informative
  post_t <- postmix(nf.prior, n = ss_trt, r = sum(sample_t)) #trt DL2: full-informative
  
  #comment out when only consider one decision rule
  p.t <- decision(post_t, post_c)#+decision.1(post_t, post_c) #full-informative post-t
  
  return(p.t)
}
simFunc = function(ss_control_FA, ss_trt_FA, orr_trt1, orr_trt2, orr_control){
  
  nk <- 10000 # number of iterations

  sel.lower <- 0
  sel.higher <- 0

  res_3arms_tmp <- c()
  res_3arms_high <- res_3arms_low <- c()
  
  ## save ind as a list
  
  ## ind: 
  ## 0: no-go
  ## 1: indeterminant
  ## 2: go

  for (i in 1:nk) {  
    print(i)
    tryCatch({
      
      set.seed(i+2024)
      #note ss_xxx_IA is the planned sample size enrolled at IA, and ss_xxx are patients who finished the second scan
      sample_c_FA <- rbinom(ss_control_FA,1,orr_control) # efficacy sample generation
      
      sample_t1_FA <- rbinom(ss_trt_FA,1,orr_trt1)
      sample_t2_FA <- rbinom(ss_trt_FA,1,orr_trt2)

      nf.prior <- mixbeta(nf.prior = c(1,0.5,0.5))
      post_t1 <- postmix(nf.prior, n = sum(ss_trt_FA), r = sum(sample_t1_FA))
      post_t2 <- postmix(nf.prior, n = sum(ss_trt_FA), r = sum(sample_t2_FA))
     
      p.diff = decision3(post_t2, post_t1)

     if (p.diff == 1) {
        ind <- 1 
        sel.higher <- sel.higher + 1
        p.res <- p.f(sample_c_FA, sample_t2_FA, delta=delta)
        res_3arms_high <- c(res_3arms_high, p.res)

      }else {
        ind <- 2
        sel.lower <- sel.lower + 1
        p.res <- p.f(sample_c_FA, sample_t1_FA, delta=delta)
        res_3arms_low <- c(res_3arms_low, p.res)
      }
      res_3arms_tmp <- c(res_3arms_tmp, p.res)
     }, error = function(e) {
      # Handle errors (skip this simulation)
      cat("Simulation", i, "encountered an error:", conditionMessage(e), "\n")
      NULL  # Return NULL as a placeholder for failed run
    })
  }
  

  #update nk (simulations with error skipped)
  nk = length(res_3arms_tmp)
  # output
  out = c(table(res_3arms_tmp)[2]/nk, 
          table(res_3arms_low)[2]/nk, 
          table(res_3arms_high)[2]/nk,
          table(res_3arms_tmp)[1]/nk) * 100
  out = t(data.frame(out))

  
  return(list(out = out))
  
}


rownames_trueORR = c("(0.55, 0.55, 0.55)", "(0.55, 0.75, 0.55)", "(0.60, 0.75, 0.55)", "(0.65, 0.75, 0.55)", 
                     "(0.73, 0.75, 0.55)", "(0.75, 0.75, 0.55)")
colnames_poc = c("Go", "DL1 Go", "DL2 Go", "No Go")

##############################################    cont_IA = 43, trt_IA = 31, FA = 60    ###########################################

t1 = Sys.time()

S1 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.55, orr_trt2=0.55, orr_control=0.55)
S2 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.55, orr_trt2=0.75, orr_control=0.55)
S3 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.60, orr_trt2=0.75, orr_control=0.55)
S4 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.65, orr_trt2=0.75, orr_control=0.55)
S5 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.73, orr_trt2=0.75, orr_control=0.55)
S6 = simFunc(ss_control_FA=60, ss_trt_FA=60, orr_trt1=0.75, orr_trt2=0.75, orr_control=0.55)

poc_full = data.frame(rbind(S1$out, S2$out, S3$out, S4$out, S5$out, S6$out))
rownames(poc_full) = rownames_trueORR
colnames(poc_full) = colnames_poc


## save outputs
poc_3arms = poc_full
save(poc_3arms, file = "./results/poc_3arms.RData")
## remove temporary data
rm(list = c(ls(pattern = "^S[1-8]$"), "poc_full", "dosepick_full"))

t2 = Sys.time()
t2 - t1

######################Analysis#############################
# load(file = "./results/poc_3arms.RData")

load(file = "./results/plot_005/poc_sam_list_cutoff1_0.05.RData")
poc_sam_005_35 <- poc_sam_list$poc_sam_ia35
poc_sam_005_45 <- poc_sam_list$poc_sam_ia45

load(file = "./results/plot_010/poc_sam_list_cutoff1_0.10.RData")
poc_sam_010_35 <- poc_sam_list$poc_sam_ia35
poc_sam_010_45 <- poc_sam_list$poc_sam_ia45

load(file = "./results/plot_015/poc_sam_list_cutoff1_0.15.RData")
poc_sam_015_35 <- poc_sam_list$poc_sam_ia35
poc_sam_015_45 <- poc_sam_list$poc_sam_ia45

poc_sam_005_35$Methods <- "N35phi005"
poc_sam_005_45$Methods <- "N45phi005"
poc_sam_010_35$Methods <- "N35phi010"
poc_sam_010_45$Methods <- "N45phi010"
poc_sam_015_35$Methods <- "N35phi015"
poc_sam_015_45$Methods <- "N45phi015"
poc_3arms$Methods <- "3arms"

df <- rbind(poc_sam_005_35,poc_sam_005_45,poc_sam_010_35,poc_sam_010_45,poc_sam_015_35,poc_sam_015_45,poc_3arms[c(1,2,6),])

legend_label_text <- c(
  N35phi005 = "paste('N'[1], ' = 105, ', 'c'[theta], ' = 0.05')",
  N45phi005 = "paste('N'[1], ' = 135, ', 'c'[theta], ' = 0.05')",
  N35phi010 = "paste('N'[1], ' = 105, ', 'c'[theta], ' = 0.10')",
  N45phi010 = "paste('N'[1], ' = 135, ', 'c'[theta], ' = 0.10')",
  N35phi015 = "paste('N'[1], ' = 105, ', 'c'[theta], ' = 0.15')",
  N45phi015 = "paste('N'[1], ' = 135, ', 'c'[theta], ' = 0.15')",
  `3arms`   = "paste('3-arm design')"
)

legend_lab_fun <- function(x) parse(text = legend_label_text[x])

df2 <- df %>%
  mutate(
    scenario = case_when(
      grepl("\\(0.55, 0.55, 0.55\\)", rownames(df)) ~ "Null",
      grepl("\\(0.55, 0.75, 0.55\\)", rownames(df)) ~ "High only effective",
      grepl("\\(0.75, 0.75, 0.55\\)", rownames(df)) ~ "Both effective",
      TRUE ~ NA_character_
    )
  )


## 2. Make 'Methods' an ordered factor so bars line up nicely
df2 <- df2 %>%
  mutate(
    Methods = factor(
      Methods,
      levels = c("N35phi005","N45phi005",
                 "N35phi010","N45phi010",
                 "N35phi015","N45phi015",
                 "3arms")
    )
  )

## 3. Split into three datasets
df_type1      <- df2 %>% filter(scenario == "Null")            %>% select(Methods, Go)
df_power_high <- df2 %>% filter(scenario == "High only effective") %>% select(Methods, `DL2 Go`)
df_power_both <- df2 %>% filter(scenario == "Both effective")  %>% select(Methods, Go)

## Helper: common theme and scales so bars are close & clean
## 4. Individual plots

# (a) Bayesian Type I error (null scenario)
p_type1 <- ggplot(df_type1, aes(x = Methods, y = Go, fill = Methods)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = NULL, y = "Bayesian Type I Error Rate") +
  theme_bw() +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 16),
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  ) +
  scale_fill_discrete(labels = legend_lab_fun)

# (b) Bayesian power when only high dose is truly effective
p_highonly <- ggplot(df_power_high, aes(x = Methods, y = `DL2 Go`, fill = Methods)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = NULL, y = "Bayesian Power \n(high dose effective)") +
  theme_bw() +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 16),
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  ) +
  scale_fill_discrete(labels = legend_lab_fun)


# (c) Bayesian power when both doses are truly effective
p_power <- ggplot(df_power_both, aes(x = Methods, y = Go, fill = Methods)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = NULL, y = "Bayesian Power \n(both doses effective)") +
  theme_bw() +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 16),
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  ) +
  scale_fill_discrete(labels = legend_lab_fun)


## 5. Combine with patchwork (3 panels side by side)
combined_plot <- (p_type1 | p_highonly| p_power) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  )
# Example: save to PDF
print(combined_plot) #10*5
