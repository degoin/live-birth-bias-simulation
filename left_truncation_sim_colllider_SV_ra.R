rm(list=ls())
library(survival)
library(tidyverse)

# simulation of live birth bias due to collider stratification 
# with social vulnerability factor that 
#   1. increases probability of exposue and 
#   2. has higher underlying risk of ptb 

# define simulation function ------------------------------------------------------------------

surv_sim <- function(iteration, s_e, lb_e, u_s_e, u_lb_e, N, superN) {
  print(paste0("s_e =",s_e," lb_e = ",lb_e)) 
  print(iteration)
  
  # generate exposure, SV factor, and unmeasured confounder
  fulldat <- data.frame(SV = rbinom(n=superN, size=1, prob=0.5))
  fulldat$U <- rbinom(n=superN, size=1, prob=0.2)            
  fulldat$A <- rbinom(n=superN, size=1, prob= 0.1 + 0.25*fulldat$SV)
  
  # generate underlying risk of SA 
  #dat$sa_risk <- ifelse(dat$SV==1, rnorm(n=N, mean=0.4, sd=0.05), rnorm(n=N, mean=0.3, sd = 0.05))
  fulldat$sa_risk <- rnorm(n=superN, mean=0.2, sd = 0.05)
  fulldat$sa_risk <- ifelse(fulldat$sa_risk<0, 0.01, ifelse(fulldat$sa_risk>1, 0.99, fulldat$sa_risk))
  
  # generate spont time 
  fulldat$spont_time <- 5 + rweibull(superN,shape=1.5,scale=0.5)*10
  
  
  # identify spontaneous abortions 
  fulldat$spont <- ifelse(fulldat$A==1, rbinom(n=superN, size=1, prob = fulldat$sa_risk + s_e + u_s_e*fulldat$U),
                      rbinom(n=superN, size=1, prob = fulldat$sa_risk + u_s_e*fulldat$U))
  
  
  # generate underlying live birth time 
  fulldat$lb_time_u <- ifelse(fulldat$SV==1, 42.5 - rweibull(superN,shape=1.75,scale=1.2)*3, 43 - rweibull(superN,shape=1.75,scale=1.2)*3)
  
  
  # generate live birth time account for exposure and social vulnerability factor 
  #dat$lb_time <- ifelse(dat$A==1 & dat$SV==1, dat$lb_time_u - lb_e - u_lb_e*dat$U - 0.05, 
  #                      ifelse(dat$A==1 & dat$SV==0, dat$lb_time_u - lb_e - u_lb_e*dat$U, 
  #                             dat$lb_time_u - u_lb_e*dat$U))
  
  fulldat$lb_time <- ifelse(fulldat$A==1, fulldat$lb_time_u - lb_e - u_lb_e*fulldat$U, 
                               fulldat$lb_time_u - u_lb_e*fulldat$U)
                       
  # if live birth time is less than 21 weeks, recode to 21 weeks
  fulldat$lb_time <- ifelse(fulldat$lb_time<21, 21, fulldat$lb_time)
  
  
  fulldat$time <- ifelse(fulldat$spont==1, fulldat$spont_time, fulldat$lb_time)
  
  fulldat$event_type <- ifelse(fulldat$spont==1, "spontaneous abortion","live birth")
  
  fulldat$event <- ifelse((fulldat$time<37 & fulldat$event_type=="live birth") | fulldat$spont==1, 1,0)
  
  
  # identify preterm births  
  fulldat$ptb <- ifelse(fulldat$time<37 & fulldat$event_type=="live birth", 1, 
                    ifelse(fulldat$event_type=="spontaneous abortion", NA, 0))
  
  
  # identify those who would have delivered preterm regarless of whether they experienced a spontaneous abortion 
  fulldat$cf_ptb <- ifelse(fulldat$lb_time<37,1,0)
  

  # additive effect of exposure on preterm birth if spontaneous abortions had not occurred 
  c_truth <- glm(cf_ptb ~ A + SV , data = fulldat)

  # stratify by social vulnerability factor 
  
  # additive effect of exposure on preterm birth if spontaneous abortions had not occurred 
  c_truth_SV0 <- glm(cf_ptb ~ A, data = fulldat[fulldat$SV==0,])
  c_truth_SV1 <- glm(cf_ptb ~ A, data = fulldat[fulldat$SV==1,])
  
  
  # now take sample 
  
  dat <- fulldat[sample(nrow(fulldat), size=N, replace = F),]
  
  # and calculate effects 
  
  # additive effect of exposure on spontaneous abortion or preterm birth 
  rd_ffit <- glm(event ~ A + SV, data=dat)
  
  # additive effect of exposure on preterm birth among births that survive
  rd_tfit <- glm(ptb ~ A + SV, data = dat[dat$event_type=="live birth",])
  
  # additive effect of exposure on preterm birth if spontaneous abortions had not occurred 
  rd_cffit <- glm(cf_ptb ~ A + SV, data = dat)
  
  # stratify by social vulnerability factor 
  
  # additive effect of exposure on spontaneous abortion or preterm birth 
  rd_ffit_SV0 <- glm(event ~ A, data=dat[dat$SV==0,])
  rd_ffit_SV1 <- glm(event ~ A, data=dat[dat$SV==1,])
  
  # additive effect of exposure on preterm birth among births that survive
  rd_tfit_SV0 <- glm(ptb ~ A, data = dat[dat$event_type=="live birth" & dat$SV==0,])
  rd_tfit_SV1 <- glm(ptb ~ A, data = dat[dat$event_type=="live birth" & dat$SV==1,])
  
  # additive effect of exposure on preterm birth if spontaneous abortions had not occurred 
  rd_cffit_SV0 <- glm(cf_ptb ~ A, data = dat[dat$SV==0,])
  rd_cffit_SV1 <- glm(cf_ptb ~ A, data = dat[dat$SV==1,])
  
  results <- data.frame(cbind(full_rd = coef(rd_ffit)["A"], full_rd_SE = coef(summary(rd_ffit))["A","Std. Error"], 
                              cf_rd = coef(rd_cffit)["A"], cf_rd_SE = coef(summary(rd_cffit))["A","Std. Error"],
                              trunc_rd = coef(rd_tfit)["A"], trunc_rd_SE = coef(summary(rd_tfit))["A", "Std. Error"], 
                              cov95 = ifelse(coef(rd_tfit)["A"] - 1.96*coef(summary(rd_tfit))["A", "Std. Error"] <= coef(c_truth)["A"] & 
                                             coef(rd_tfit)["A"] + 1.96*coef(summary(rd_tfit))["A", "Std. Error"] >=coef(c_truth)["A"], 1, 0),
                              
                              full_rd_SV0 = coef(rd_ffit_SV0)["A"], full_rd_SE_SV0 = coef(summary(rd_ffit_SV0))["A","Std. Error"], 
                              full_rd_SV1 = coef(rd_ffit_SV1)["A"], full_rd_SE_SV1 = coef(summary(rd_ffit_SV1))["A","Std. Error"], 
                              
                              cf_rd_SV0 = coef(rd_cffit_SV0)["A"], cf_rd_SE_SV0 = coef(summary(rd_cffit_SV0))["A","Std. Error"],
                              cf_rd_SV1 = coef(rd_cffit_SV1)["A"], cf_rd_SE_SV1 = coef(summary(rd_cffit_SV1))["A","Std. Error"],
                              
                              trunc_rd_SV0 = coef(rd_tfit_SV0)["A"], trunc_rd_SE_SV0 = coef(summary(rd_tfit_SV0))["A", "Std. Error"],
                              trunc_rd_SV1 = coef(rd_tfit_SV1)["A"], trunc_rd_SE_SV1 = coef(summary(rd_tfit_SV1))["A", "Std. Error"], 
                              
                              cov95_SV0 = ifelse(coef(rd_tfit_SV0)["A"] - 1.96*coef(summary(rd_tfit_SV0))["A", "Std. Error"] <= coef(c_truth_SV0)["A"] & 
                                               coef(rd_tfit_SV0)["A"] + 1.96*coef(summary(rd_tfit_SV0))["A", "Std. Error"] >=coef(c_truth_SV0)["A"], 1, 0),
                               
                              cov95_SV1 = ifelse(coef(rd_tfit_SV1)["A"] - 1.96*coef(summary(rd_tfit_SV1))["A", "Std. Error"] <= coef(c_truth_SV1)["A"] & 
                                                   coef(rd_tfit_SV1)["A"] + 1.96*coef(summary(rd_tfit_SV1))["A", "Std. Error"] >=coef(c_truth_SV1)["A"], 1, 0),
                              
                              s_e = s_e, lb_e=lb_e, u_s_e=u_s_e, u_lb_e = u_lb_e, mean_ptb = mean(dat$ptb[dat$event_type=="live birth"]), 
                              mean_sa = mean(dat$spont), mean_A_all = mean(dat$A), mean_A_lb = mean(dat$A[dat$event_type=="live birth"]), 
                              mean_SV_all = mean(dat$SV), mean_SV_lb = mean(dat$SV[dat$event_type=="live birth"])))
  
  
  return(results) 
}

# effect on probability of induced abortion range from 0 to 0.005 by 0.001
# effect on probability of spontaneous abortion range from 0 to 0.01 by 0.0025
# effect on gestational age on scale of live birth range from -0.05 to 0.05  by 0.05

# run the simulation n times and take the mean of the parameters 
mc <- function(n, s_e, lb_e, u_s_e, u_lb_e, N, superN=superN) {
ls <- lapply(1:n, function(x) surv_sim(x, s_e = s_e, lb_e = lb_e, u_s_e=u_s_e, u_lb_e = u_lb_e, N = N, superN=superN))
df <- do.call(rbind, ls)
r <- colMeans(df)
return(r)
}


# run simulation across values for spontaneous and live birth effects ------------------------------------------------------------------


surv_ls <- lapply(seq(0, 0.2, by=0.01), function(y) 
                    lapply(seq(-0.05, 0.05, by=0.05), function(z) mc(n = 20, s_e = y, lb_e = z, u_s_e=0.2, u_lb_e=0.5, N=625000, superN=5000000)))

surv_ls2 <- lapply(1:length(seq(0, 0.2, by=0.01)), function(x) do.call(rbind, surv_ls[[x]]))

surv_m  <- data.frame(do.call(rbind, surv_ls2))


write.csv(surv_m, file="/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/sim_results_collider_SV_ra_1000_v2.csv", row.names = F)

surv_m <- read.csv("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/sim_results_collider_SV_ra_1000_v2.csv")


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# plots 


library(gridExtra)


p1 <-   ggplot(surv_m) + 
  theme_bw() + 
  facet_wrap(~lb_e, labeller = labeller(lb_e = c("-0.05" ="Lengthen gestation time", "0" ="No effect on gestation time", "0.05"="Shorten gestation time"))) + 
  labs(title="Effect of exposure on time to live birth", 
       x="Effect of exposure on \n probability of spontaneous abortion") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        plot.margin = unit(c(1,1.5,1,1.5), "cm"), 
        axis.text.x = element_text(size=14), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=14), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=14)) + 
  scale_linetype_manual("", labels = c("Counterfactual", "Truncated"), values=c(6,5)) + 
  scale_fill_manual("", labels = c("Counterfactual", "Truncated"), values=c("#1a9850","#4575b4")) +
  geom_ribbon(aes(x = s_e, ymin=cf_rd - 1.96*cf_rd_SE, ymax=cf_rd + 1.96*cf_rd_SE, fill="Counterfactual"), alpha=0.6) + 
  geom_ribbon(aes(x = s_e, ymin=trunc_rd - 1.96*trunc_rd_SE, ymax=trunc_rd + 1.96*trunc_rd_SE, fill="Truncated"), alpha=0.5) + 
  geom_hline(yintercept=0, linetype=3) + 
  geom_line(aes(x=s_e, y=cf_rd, linetype="Counterfactual"), size=1.5) + 
  geom_line(aes(x=s_e, y=trunc_rd, linetype="Truncated"), size=1.5)  + 
  scale_y_continuous("Preterm birth risk difference per 100 women",  labels=function(x) sprintf("%.1f", x*100), limits=c(-0.012, 0.012))
p1
ggsave(p1, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_SV_ra_1000.pdf"), width=15)


p2a <- ggplot(surv_m %>% filter(lb_e==0.05)) + 
  theme_bw() + 
  labs(title="Without social vulnerability", x="Effect of exposure on \n probability of spontaneous abortion") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        plot.margin = unit(c(1,1.5,1,1.5), "cm"), 
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16)) + 
  scale_linetype_manual("", breaks = c("Counterfactual", "Truncated"), values=c(6,5)) + 
  scale_fill_manual("", labels = c("Counterfactual", "Truncated"), values=c("#636363","#bdbdbd")) +
  geom_ribbon(aes(x = s_e, ymin=cf_rd_SV0 - 1.96*cf_rd_SE_SV0, ymax=cf_rd_SV0 + 1.96*cf_rd_SE_SV0, fill="Counterfactual"), alpha=0.6) + 
  geom_ribbon(aes(x = s_e, ymin=trunc_rd_SV0 - 1.96*trunc_rd_SE_SV0, ymax=trunc_rd_SV0 + 1.96*trunc_rd_SE_SV0, fill="Truncated"), alpha=0.4) +
 guides(linetype=FALSE, fill=FALSE) + 
  geom_line(aes(x=s_e, y=cf_rd_SV0, linetype="Counterfactual"), size=1.5) + 
  geom_line(aes(x=s_e, y=trunc_rd_SV0, linetype="Truncated"), size=1.5) + 
  scale_y_continuous("Preterm birth risk difference per 100 women",  labels=function(x) sprintf("%.1f", x*100), limits=c(-0.005, 0.01)) 

p2a

ggsave(p2a, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_SV0_ra_1000_lb05.pdf"), width=10)


p3a <- ggplot(surv_m %>% filter(lb_e == 0.05)) + 
  theme_bw() + 
  labs(title="With social vulnerability", x="Effect of exposure on \n probability of spontaneous abortion",  
       y="PTB risk difference") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16), 
        plot.margin = unit(c(1,0,1,0), "cm")) + 
  scale_linetype_manual("", labels = c("Counterfactual", "Truncated"), values=c(6,5)) + 
  scale_fill_manual("", labels = c("Counterfactual", "Truncated"), values=c("#636363","#bdbdbd")) +
  geom_ribbon(aes(x = s_e, ymin=cf_rd_SV1 - 1.96*cf_rd_SE_SV1, ymax=cf_rd_SV1 + 1.96*cf_rd_SE_SV1, fill="Counterfactual"), alpha=0.6) + 
  geom_ribbon(aes(x = s_e, ymin=trunc_rd_SV1 - 1.96*trunc_rd_SE_SV1, ymax=trunc_rd_SV1 + 1.96*trunc_rd_SE_SV1, fill="Truncated"), alpha=0.4) + 
  geom_line(aes(x=s_e, y=cf_rd_SV1, linetype="Counterfactual"), size=1.5) + 
  geom_line(aes(x=s_e, y=trunc_rd_SV1, linetype="Truncated"), size=1.5) + 
  scale_y_continuous("Preterm birth risk difference per 100 women",  labels=function(x) sprintf("%.1f", x*100), limits=c(-0.005, 0.01)) 


p3a
ggsave(p3a, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_SV1_ra_1000_lb05.pdf"), width=10)


p_ab <- grid.arrange(p2a, p3a, nrow=1)

ggsave(p_ab, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_SV_ra_1000_lb05.pdf"), width=15)


p4 <- ggplot(surv_m) + 
  geom_line(aes(x=s_e, y=cf_rd_SV0 - trunc_rd_SV0, linetype="Without social vulnerability"), size=1.5) + 
  geom_line(aes(x=s_e, y=cf_rd_SV1 - trunc_rd_SV1, linetype="With social vulnerability"), size=1.5) + 
  theme_bw() + 
  facet_wrap(~lb_e, labeller = labeller(lb_e = c("-0.05" ="Lengthen gestation time", "0" ="No effect on gestation time", "0.05"="Shorten gestation time"))) + 
  labs(title="Effect of exposure on time to live birth", 
       x="Effect of exposure on \n probability of spontaneous abortion",  
       y="Absolute bias") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16)) +
  scale_linetype_manual("", breaks = c("Without social vulnerability", "With social vulnerability"), labels=c("Without social vulnerability","With social vulnerability"), values=c(1,4)) + 
  geom_hline(yintercept=0, linetype=3) + ylim(-0.0025, 0.005)

p4
ggsave(p4, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_ra_bias_1000.pdf"), width=15)


p4a <- ggplot(surv_m %>% filter(lb_e==0.05)) + 
  geom_line(aes(x=s_e, y=cf_rd_SV0 - trunc_rd_SV0, linetype="Without social vulnerability"), size=1.5) + 
  geom_line(aes(x=s_e, y=cf_rd_SV1 - trunc_rd_SV1, linetype="With social vulnerability"), size=1.5) + 
  theme_bw() + 
  labs(x="Effect of exposure on \n probability of spontaneous abortion",  
       y="Absolute bias") +
  theme(plot.title = element_text(hjust = 0.5),  
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16)) + 
  scale_linetype_manual("", breaks =c("Without social vulnerability","With social vulnerability"), values=c(1,4)) + 
  geom_hline(yintercept=0, linetype=3) + ylim(-0.001, 0.002)

p4a
ggsave(p4a, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_ra_bias_1000_lb05.pdf"), width=10)

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# coverage plots


p1cov <-   ggplot(surv_m) + 
  geom_line(aes(x=s_e, y=cov95, linetype="coverage"), size=1.5) + 
  theme_bw() + 
  facet_wrap(~lb_e, labeller = labeller(lb_e = c("-0.05" ="Lengthen gestation time", "0" ="No effect on gestation time", "0.05"="Shorten gestation time"))) + 
  labs(title="Effect of exposure on time to live birth", 
       x="Effect of exposure on \n probability of spontaneous abortion",  
       y="95% CI coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        plot.margin = unit(c(1,1.5,1,1.5), "cm"), 
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16)) + 
  scale_linetype_manual("", breaks = c("coverage"), values=c(1)) + 
  geom_hline(yintercept=0.95, linetype=3) + ylim(0.8, 1)

p1cov
ggsave(p1cov, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_ra_cov95_1000.pdf"), width=15)


p2cov <- ggplot(surv_m) + 
  geom_line(aes(x=s_e, y=cov95_SV0, linetype="Without social vulnerability"), size=1.5) + 
  geom_line(aes(x=s_e, y=cov95_SV1, linetype="With social vulnerability"), size=1.5) + 
  theme_bw() + 
  facet_wrap(~lb_e, labeller = labeller(lb_e = c("-0.05" ="Lengthen gestation time", "0" ="No effect on gestation time", "0.05"="Shorten gestation time"))) + 
  labs(title="Effect of exposure on time to live birth", 
       x="Effect of exposure on \n probability of spontaneous abortion",  
       y="95% CI coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        plot.margin = unit(c(1,1.5,1,1.5), "cm"), 
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16)) + 
  scale_linetype_manual("", breaks = c("Without social vulnerability", "With social vulnerability"), values=c(1,4)) + 
  geom_hline(yintercept=0.95, linetype=3) + ylim(0.8, 1)

p2cov
ggsave(p2cov, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_ra_cov95_SV_1000.pdf"), width=15)


p2cov_a <- ggplot(surv_m %>% filter(lb_e==0.05)) + 
  geom_line(aes(x=s_e, y=cov95_SV0, linetype="Without social vulnerability"), size=1.5) + 
  geom_line(aes(x=s_e, y=cov95_SV1, linetype="With social vulnerability"), size=1.5) + 
  theme_bw() + 
  labs(x="Effect of exposure on \n probability of spontaneous abortion",  
       y="95% confidence interval coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),  
        axis.text.x = element_text(size=16), 
        axis.title.x=element_text(size=18),
        axis.text.y = element_text(size=16), axis.title.y= element_text(size=18), 
        strip.text.x = element_text(size = 18), 
        legend.text = element_text(size=16))  + 
  scale_linetype_manual("", breaks = c("Without social vulnerability", "With social vulnerability"), values=c(1,4)) + 
  geom_hline(yintercept=0.95, linetype=3) + ylim(0.8, 1)

ggsave(p2cov_a, file=paste0("/Users/danagoin/Documents/Research projects/Left truncation of birth cohorts simulation study/results/plots/sim_plot_collider_ra_cov95_SV_1000_lb05.pdf"), width=10)


    