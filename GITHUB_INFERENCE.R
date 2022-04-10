library(tidyverse)
library(data.table)
library(psych)
library(EnvStats)
library(xtable)
library(knitr)


library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(knitr)
library(dplyr)
library(data.table)
library(tidyverse)
library(bit64)
library(dplyr)
library(readxl)
library(tidyverse)
library(data.table)

rm(list = ls())


dir_path <- getwd()

 likelihood2 <- function(x1,x2,p2,rho0){
    thetai <- p2 + rho0 *(x1 - p2)
    return(sum(x2*log(thetai)+(1-x2)*log(1-thetai),na.rm=T))
 }
 
 
 
 # Unifor hyper ------------
 
 set.seed(10032021)
 k    <- 1000
 p1   <- 0.50
 p2   <- 0.75
 rho1 <- 0.70
  i <- 1.00
  j <- 1.00
 p2 <- 0.25
 p1 <- 0.50
 base_time_full <- NULL
 base_time_posteriori_full <- NULL
 for(p1 in seq(0,1,by=.25)){
    for(p2 in seq(0,1,by=.25)){
       for(rho1 in seq(0,1,by=.2)){
          if(2>1){
          base_time_seq <- NULL
          for(i in 1:(k/100)){
             
             x1 <- rbinom(k,1,p1)
             x2 <- rbinom(k,1,p2 + rho1*(x1*(1-p2) - (1-x1)*p2))

             base_time <- data.frame(   rhhov=rho1,
                                         p1=p1,p2=p2,
                                         mx2obs=mean(x2),
                                         mx2theory=p2 + rho1*(p1-p2),
                                         vx2obs= var(x2),
                                         vx2theory= p2  - p2^2 + rho1*(p1 - p2 - 2*p2*(p1-p2)) - (rho1^2)*((p1-2*p1*p2+p2^2)-p1*(1-p1)),
                                         exytheory = p1*(p2 + rho1*( 1 - p2)),
                                         exyobs = mean(x1*x2),
                                         rhoestmomentos = (mean(x1*x2) - p1*p2)/(p1*(1-p2)),
                                         it = i)
             
             
             # sampling importance
             # distribution support
             rho0 <- runif(k)
             res2 <- matrix(NA,nrow = length(rho0),ncol = 1)
             for(j in 1:nrow(res2)){
                res2[j,1] <- likelihood2(x1,x2,p2,rho0[j])
             }
             a <- res2 %>% as.data.frame() 
             #  importance
             a$probboots <- exp(a$V1)/sum(exp(a$V1),na.rm = T)
             #  resampling 
             samplebootstrap <- sample(1:length(a$probboots),size=k,prob=a$probboots,replace = TRUE)
             rho0 <- rho0[samplebootstrap]
             base_time$simrho_med <- median(rho0)
             base_time$simrho_0250 <- quantile(rho0,p=.025)
             base_time$simrho_0975 <- quantile(rho0,p=.975)
             
             base_time_seq <- bind_rows(base_time_seq,base_time)
             
             if(i == 5){
                base_time_posteriori_full <- bind_rows(base_time_posteriori_full,
                                                       data.frame(p1=p1,
                                                                  p2=p2,
                                                                  rho1=rho1,
                                                                  simula = rho0))   
             }
          }
          base_time_full <- bind_rows(base_time_full,base_time_seq)
       }
       }
    }
    print(p1)
 }

 base_time_full_t <- base_time_full %>%
    select(it,rhhov,rhoestmomentos,simrho_med,p1,p2) %>%
    gather(key="methods",value="value",-rhhov,-it,-p1,-p2)
 
 
 
 vs <- c("Moments","Bayesian")
 names(vs) <- c("rhoestmomentos","simrho_med")
 
 Ss <- c(expression(delta[2],"=",0),expression(delta[2],"=",0),expression(delta[2],"=",0),expression(delta[2],"=",0),expression(delta[2],"=",0))
 Ss <- paste0("P2=",c(0,.25,.50,.75,1))
 names(Ss) <- c("0" ,"0.25","0.5", "0.75", "1")
 
 
 Ss2 <- paste0("P1=",c(0,.25,.50,.75,1))
 names(Ss2) <- c("0" ,"0.25","0.5", "0.75", "1")
 
 
   ggplot(base_time_full_t %>% filter(rhhov!=1) %>%
             mutate(p2 = as.character(p2)) %>%
             filter(!((p2==1 & p1==0) | (p2==1 & p1==1)| (p2==0 & p1==1)| (p2==0 & p1==0))),
          aes(x=(rhhov) %>% as.character(),y=value-rhhov,
              fill=methods %>% as.character())) +
      scale_fill_manual(breaks = names(vs),
                        labels = as.character(vs),
                        values = c("salmon","lightskyblue"))+
      geom_boxplot(outlier.size = .1)+
      geom_hline(yintercept = 0,linetype="dashed",col="red")+
      xlab(expression(delta[1])) + ylab(expression(hat(delta)-delta[1])) +
      facet_grid(p2~p1,scales = "free",
                 labeller = labeller(p1=Ss2,
                                     p2=Ss))+
      labs(fill=NULL)+ theme(legend.position="bottom")
   
   