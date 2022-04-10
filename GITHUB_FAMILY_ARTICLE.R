rm(list = ls());gc()
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



dir_path <- getwd()

sourceCpp("correl_Bernoulli.cpp")

data_death_rp_2000 <- 1



github_mortality_woman_rp_2000 <- readRDS(file.path(dir_path,"bases","github_mortality_woman_rp_2000.rds"))
github_mortality_man_rp_2000 <- readRDS(file.path(dir_path,"bases","github_mortality_man_rp_2000.rds"))

gender_loop <- 1 #male
datatemp_full <- NULL
age <- 40
age2 <- age-3
rho1 <- 0
rho2 <- 0
it <- 1
set.seed(5032021)
for(rho1 in c(0,0.01,0.1,.5,1)){
  print(paste0(rho1," e iteration:",it))
  t0 <- Sys.time()
    for(rho2 in c(0,0.01,0.1)){
      age_part <- age - 1
      age_max_conjuge <- 120
      age_max_daughter <- 30
      gender_conj_loop <- 3-gender_loop 
      age_conj_orig <-ifelse(gender_loop==1,age-3,age+3)
      age_conj_orig <- ifelse(age_conj_orig>age_max_conjuge,-1,age_conj_orig-1)
      
      
      age_daughter_1 <- ifelse(gender_loop==1,age_conj_orig-35,age-35)
      age_daughter_1 <- ifelse(age_daughter_1>age_max_daughter,-1,age_daughter_1)
      
      
      
      vec_mortality_woman_rp_2000          <- github_mortality_woman_rp_2000$mortality
      
      # daughter
      vec_mortality_daughter_rp_2000<- github_mortality_woman_rp_2000$mortality
      
      
      
      # man
      vec_mortality_man_rp_2000 <- github_mortality_man_rp_2000$mortality
      
      
      
      
      
        vec_mortality_woman_rp_2000_t1 <- vec_mortality_woman_rp_2000
        vec_mortality_daughter_rp_2000_t1 <- vec_mortality_daughter_rp_2000

      vecpar <- 0:120
      vecconj <- vecpar<=age_max_conjuge & vecpar>=age_conj_orig
      vecfil <-  vecpar<=age_max_daughter & vecpar>=age_daughter_1
      vecpar <- vecpar>=age_part
      
      
      vec_mortality_woman_rp_2000_t1 <-   vec_mortality_woman_rp_2000_t1[vecconj]
      vec_mortality_daughter_rp_2000_t1 <-   vec_mortality_daughter_rp_2000_t1[vecfil]
      vec_mortality_man_rp_2000_t1 <- vec_mortality_man_rp_2000[vecpar]
      
      maxtemp <- max(length(vec_mortality_man_rp_2000_t1),length(vec_mortality_woman_rp_2000_t1),length(vec_mortality_daughter_rp_2000_t1))
      
      
      matrixprobs <- matrix(1,nrow=maxtemp,ncol=3)
      matrixprobs[1:length(vec_mortality_man_rp_2000_t1),1] <- vec_mortality_man_rp_2000_t1
      matrixprobs[1:length(vec_mortality_woman_rp_2000_t1),2] <- vec_mortality_woman_rp_2000_t1
      matrixprobs[1:length(vec_mortality_daughter_rp_2000_t1),3] <- vec_mortality_daughter_rp_2000_t1

      
ks <- 40000000

a <- func_bernoulli_correl_fam(matrixprobs,
                                        ks,
                                        rho1,
                                        rho2,
                                        5)



probss <- -(a[,-1]-a[,-ncol(a)])/a[,-ncol(a)]
probss <- ifelse(is.na(probss),1,probss)
probss_full_mean <- apply(probss,2,mean)

probss_full_median <- apply(probss,2,median)
probss_full_q05 <- apply(probss,2,function(x)quantile(x,.05))
probss_full_q95 <- apply(probss,2,function(x)quantile(x,.95))


datatemp <- data.frame(prob = probss_full_mean,
                       probabilage_med = probss_full_median,
                       probabilage_05 = probss_full_q05,
                       probabilage_95 = probss_full_q95,
                       gender = gender_loop,
                       age = age,
                       genderconj = gender_conj_loop,
                       idadeconj = age_conj_orig,
                       rho1 = rho1,
                       rho2 = rho2,
                       probnormalman= matrixprobs[-1,1],
                       probnormalwoman = matrixprobs[-1,2]) %>%
  mutate(idadeplot = row_number()-1)

datatemp_full <- rbind(datatemp_full,datatemp)

  }
  

  t1 <- Sys.time()
  print(t1 - t0)
}


vs <- c(30,40,50,60,70,80,90)
vs <- paste0("Man ",vs," Woman ",vs-3," Daughter ",vs-38)
names(vs) <- c(30,40,50,60,70,80,90)

ggplot(datatemp_full %>% filter(age==40 &  (idadeplot+idadeconj+1)<100 & rho2!=.01)) +
  geom_line(aes(x=(idadeplot+idadeconj)+1,
                y=log(probabilage_med),
                col=rho1 %>% as.character(),
                linetype=rho2 %>% as.character()),
            size=1) +
  labs(col=expression(delta[1]),
       linetype=expression(delta[2])) +
  ylab("log mortality") +
  xlab("Woman age") +  
    scale_x_continuous(breaks = c(seq(37,100,by=5),100))+ 
  theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 1))



txjuros <- 2/100
aa <- datatemp_full %>% filter(idadeplot>=(idadeconj)) %>% 
  mutate(qts = exp(cumsum(log(1-probabilage_med)))) 

fff <- function(x){
  return(rev(cumsum(rev(exp(cumsum(log(1-x))))))[1])
}


ggg <- function(x,txjuros,idadeplot,idade_conj){
  return(sum(exp(cumsum(log(1-x)))*(1/(1+txjuros))^((1:length(x)-1))))
}


aa <-  datatemp_full %>% 
  group_by(rho1,rho2)%>%
  summarise(expects = fff(probabilage_med),
            fator = ggg(probabilage_med,txjuros))


res1 <- aa %>% select(rho1,rho2,expects)%>%
  mutate(expects = round(expects,2)) %>%
  spread(key = "rho1",value="expects")

print(xtable::xtable(res1),include.rownames = F)

