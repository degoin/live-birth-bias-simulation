rm(list=ls())
library(survival)
library(tidyverse)
#library(lattice)
#library(latticeExtra)
# simulate 20% of full population who had either induced or spontaneous abortion before week 22 

# generate time to event (before 22 weeks) based off exposure status and relationship with exposure 

# vary strength of association with homicide exposure 

# calculate association using time to event models that include left-truncated simulated data 

# compare with estimate without left-truncated data to get a range of possible bias scenarios


# define simulation function ------------------------------------------------------------------

surv_sim <- function(iteration, s_e, lb_e, N) {
  print(paste0("s_e =",s_e," lb_e = ",lb_e)) 
  print(iteration)
  
  # generate exposure and unmeasured confounder
  dat <- data.frame(A =rbinom(n=N, size=1, prob=0.5))
  
  # generate underlying risk of SA 
  dat$sa_risk <- rnorm(n=N, mean=0.3, sd = 0.05)
  dat$sa_risk <- ifelse(dat$sa_risk<0, 0.01, ifelse(dat$sa_risk>1, 0.99, dat$sa_risk))
  
  # generate spont time 
  dat$spont_time <- 5 + rweibull(N,shape=1.5,scale=0.5)*10

  # identify spontaneous abortions 
  dat$spont <- ifelse(dat$A==1, rbinom(n=N, size=1, prob = dat$sa_risk + s_e), 
                      rbinom(n=N, size=1, prob = dat$sa_risk))
  
  
  # generate underlying live birth time 
  dat$lb_time_u <- 43 - rweibull(N,shape=1.75,scale=1.2)*3
  
  
  # generate underlying live birth time 
  dat$lb_time <- ifelse(dat$A==1, dat$lb_time_u - lb_e, dat$lb_time_u)
                       
  # if live birth time is less than 21 weeks, recode to 21 weeks
  dat$lb_time <- ifelse(dat$lb_time<21, 21, dat$lb_time)
  
  
  dat$time <- ifelse(dat$spont==1, dat$spont_time, dat$lb_time)
  
  dat$event_type <- ifelse(dat$spont==1, "spontaneous abortion","live birth")
  
  dat$event <- ifelse((dat$time<37 & dat$event_type=="live birth") | dat$spont==1, 1,0)
  
  
  # identify preterm births  
  dat$ptb <- ifelse(dat$event==1 & dat$event_type=="live birth", 1, 
                    ifelse(dat$event_type=="spontaneous abortion", NA, 0))
  
  
  # identify those who would have delivered preterm regarless of whether they experienced a spontaneous abortion 
  dat$cf_ptb <- ifelse(dat$lb_time<37,1,0)
  
  # additive effect of exposure on spontaneous abortion or preterm birth 
  rd_ffit <- glm(event ~ A, data=dat)

  # additive effect of exposure on preterm birth among births that survive
  rd_tfit <- glm(ptb ~ A, data = dat[dat$event_type=="live birth",])

  # additive effect of exposure on preterm birth if spontaneous abortions had not occurred 
  rd_cffit <- glm(cf_ptb ~ A, data = dat)

  
  results <- data.frame(cbind(full_rd = coef(rd_ffit)["A"], full_rd_SE = coef(summary(rd_ffit))["A","Std. Error"], 
                              cf_rd = coef(rd_cffit)["A"], cf_rd_SE = coef(summary(rd_cffit))["A","Std. Error"],
                              trunc_rd = coef(rd_tfit)["A"], trunc_rd_SE = coef(summary(rd_cffit))["A", "Std. Error"], 
                              s_e = s_e, lb_e=lb_e))
  
  
  return(results) 
}

# effect on probability of induced abortion range from 0 to 0.005 by 0.001
# effect on probability of spontaneous abortion range from 0 to 0.01 by 0.0025
# effect on gestational age on scale of live birth range from -0.05 to 0.05  by 0.05

# run the simulation n times and take the mean of the parameters 
mc <- function(n, s_e, lb_e, N) {
ls <- lapply(1:n, function(x) surv_sim(x, s_e = s_e, lb_e = lb_e, N = N))
df <- do.call(rbind, ls)
r <- colMeans(df)
return(r)
}


# run simulation across values for spontaneous and live birth effects ------------------------------------------------------------------


surv_ls <- lapply(seq(0, 0.2, by=0.01), function(y) 
                    lapply(seq(-0.05, 0.05, by=0.05), function(z) mc(n = 1000, s_e = y, lb_e = z, N=625000)))

surv_ls2 <- lapply(1:length(seq(0, 0.2, by=0.01)), function(x) do.call(rbind, surv_ls[[x]]))

surv_m  <- data.frame(do.call(rbind, surv_ls2))

#surv_m$s_e <- c(rep(0.00, 3), rep(0.01, 3), rep(0.02, 3), rep(0.03, 3), rep(0.04, 3), rep(0.05, 3), rep(0.06, 3), rep(0.07, 3), rep(0.08, 3), rep(0.09,3), rep(0.10, 3), 
 #               rep(0.11, 3), rep(0.12, 3), rep(0.13, 3), rep(0.14, 3), rep(0.15, 3), rep(0.16, 3), rep(0.17, 3), rep(0.18, 3), rep(0.19, 3), rep(0.20, 3))
#surv_m$lb_e <- rep(c(-0.05,0,0.05), length(seq(0,0.2, by=0.01)))


write.csv(surv_m, file="/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/sim_results_nobias.csv", row.names = F)



surv_m <- read.csv("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/sim_results_nobias.csv")



#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------



p <- ggplot(surv_m) +  
  geom_line(aes(x=s_e, y=full_rd, linetype="Full"), size=1.5) +
  geom_line(aes(x=s_e, y=trunc_rd, linetype="Truncated"), size=1.5) + 
  theme_bw() + 
  facet_wrap(~lb_e, labeller = labeller(lb_e = c("-0.05" ="Lengthen gestation time", "0" ="No effect on gestation time", "0.05"="Shorten gestation time"))) + 
  labs(title="", 
       x="Effect of exposure on \n probability of spontaneous abortion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12)) + 
  scale_linetype_manual("", breaks = c("Full", "Truncated"), values=c(4,5)) + 
  geom_ribbon(aes(x = s_e, ymin=full_rd - 1.96*full_rd_SE, ymax=full_rd + 1.96*full_rd_SE), alpha=0.4)  + 
  geom_ribbon(aes(x = s_e, ymin=trunc_rd - 1.96*trunc_rd_SE, ymax=trunc_rd + 1.96*trunc_rd_SE), alpha=0.4)  + 
  geom_hline(yintercept=0, linetype=3) + scale_y_continuous("risk difference per 100 women", labels=function(x) sprintf("%.1f", x*100))


ggsave(p, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_nobias.pdf"), width=15)



p2 <- ggplot(surv_m %>% filter(lb_e==0.05)) +  
  geom_line(aes(x=s_e, y=full_rd, linetype="All conceptions"), size=1.5) +
  geom_line(aes(x=s_e, y=trunc_rd, linetype="Live births"), size=1.5) + 
  theme_bw() + 
  labs(title="", 
       x="Effect of exposure on \n probability of spontaneous abortion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=14), 
        axis.title.x=element_text(size=14),
        axis.text.y = element_text(size=14), axis.title.y= element_text(size=14), 
        strip.text.x = element_text(size = 14), 
        legend.text = element_text(size=14)) + 
  scale_linetype_manual("", breaks = c("All conceptions", "Live births"), values=c(4,5))  + 
  geom_hline(yintercept=0, linetype=3) + scale_y_continuous("risk difference per 100 women",  labels=function(x) sprintf("%.1f", x*100))


ggsave(p2, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_nobias_lb05.pdf"), width=10)



