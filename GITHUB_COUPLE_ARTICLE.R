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
    age_part <- age - 1
    age_max_conjuge <- 120
    age_max_daughter <- 30
    gender_conj_loop <- 3-gender_loop 
    age_conj_orig <-ifelse(gender_loop==1,age-3,age+3)
    age_conj_orig <- ifelse(age_conj_orig>age_max_conjuge,-1,age_conj_orig-1)
    
    vec_mortality_woman_rp_2000          <- github_mortality_woman_rp_2000$mortality
    
    # man
    vec_mortality_man_rp_2000 <- github_mortality_man_rp_2000$mortality
    
    
    vec_mortality_woman_rp_2000_t1 <- vec_mortality_woman_rp_2000
    vecpar <- 0:120
    vecconj <- vecpar<=age_max_conjuge & vecpar>=age_conj_orig
    vecfil <-  vecpar<=age_max_daughter & vecpar>=age_daughter_1
    vecpar <- vecpar>=age_part
    
    
    vec_mortality_woman_rp_2000_t1 <-   vec_mortality_woman_rp_2000_t1[vecconj]
        vec_mortality_man_rp_2000_t1 <- vec_mortality_man_rp_2000[vecpar]
    
    maxtemp <- max(length(vec_mortality_man_rp_2000_t1),length(vec_mortality_woman_rp_2000_t1))
    
    
    matrixprobs <- matrix(1,nrow=maxtemp,ncol=2)
    matrixprobs[1:length(vec_mortality_man_rp_2000_t1),1] <- vec_mortality_man_rp_2000_t1
    matrixprobs[1:length(vec_mortality_woman_rp_2000_t1),2] <- vec_mortality_woman_rp_2000_t1
    
    ks <- 1000000

a <- func_bernoulli_correl_coup(matrixprobs,
                                         ks,
                                         rho1,
                                         100)




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
                       probnormalman= matrixprobs[-1,1],
                       probnormalwoman = matrixprobs[-1,2]) %>%
  mutate(idadeplot = row_number()-1)

datatemp_full <- rbind(datatemp_full,datatemp)
t1 <- Sys.time()
print(t1 - t0)
  }

vs <- c(60,70,80,90)
vs <- paste0("Man ",vs," Woman ",vs-3)
names(vs) <- c(60,70,80,90)
datatemp_full_filter <- datatemp_full %>% filter((idadeplot+idadeconj+1) <100)


ggplot(datatemp_full_filter) +
  geom_line(aes(x=idadeplot+idadeconj+1 ,
                y=log(prob),
                col=rho1 %>% as.character()),
            size=1) +
  geom_line(aes(x=idadeplot+idadeconj +1,
                y=log(probnormalman),
                col="Man"),
            size=1,
            linetype="dotted") +
  labs(col=expression(delta[1]),
       linetype="Man") +
  ylab("log mortality") +
  xlab("Woman age") + 
  scale_color_manual(values=c("red", "blue", "green","orange","gray","black"))+
  scale_x_continuous(breaks = c(seq(57,97,by=5),100)) + theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 1))

ggsave(file.path(dir_path,"resultados","REVISAO","modelo_last_survivor_60.pdf"),
       height = 4,width = 8)

ggsave(file.path(dir_path,"resultados","REVISAO","modelo_last_survivor_60.tiff"),
       height = 4,width = 8)


txjuros <- 2/100
aa <- datatemp_full  %>% 
  mutate(qts = exp(cumsum(log(1-probabilidade_med)))) 

fff <- function(x){
  return(rev(cumsum(rev(exp(cumsum(log(1-x))))))[1])
}


ggg <- function(x,txjuros,idadeplot,idade_conj){
  return(sum(exp(cumsum(log(1-x)))*(1/(1+txjuros))^((1:length(x)-1))))
}

aa <-  datatemp_full  %>% 
  group_by(age,idadeconj,rho1) %>%
  summarise(expects = fff(probabilage_med),
            fator = ggg(probabilage_med,txjuros))


res1 <- aa %>% select(age,idadeconj,rho1,expects)%>% 
  mutate(expects = round(expects,2),
         age = as.character(age),
         idadeconj = as.character(idadeconj)) %>%
  spread(key = "rho1",value="expects")

print(xtable::xtable(res1),include.rownames = F)


res2 <- aa %>% select(age,idadeconj,rho1,fator) %>% 
  mutate(fator = paste0("$",round(10^6/fator,2)),
         age = as.character(age),
         idadeconj = as.character(idadeconj)) %>%
  spread(key = "rho1",value="fator")

print(xtable::xtable(res2),include.rownames = F)


res3 <- aa %>% select(age,idadeconj,rho1,fator,expects) %>% filter(rho1 != .01) %>%
  mutate(valor = paste0("$",round(10^6/fator,2)," (",round(expects,2),")")) %>%
  select(-fator,-expects) %>%
  spread(key = "rho1",value="valor") %>%
  filter(age %in% c(60,70,80,90) )
print(xtable::xtable(res3),include.rownames = F)

