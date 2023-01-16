###
library(betapart)
library(adespatial)
library(vegan)
library(BAT)
library(car)
library(dplyr)
library(tidyverse)
library(plyr)
library(betaC)
library(rstatix)
library(ggplot2)

source("site_dist.R")

set.seed(1000)

dissim = "jaccard"
sample_pattern <- data.frame(A=c(20,10,5,10),B=c(20,10,15,30))
sample_label <- paste0("(",sample_pattern$A,",",sample_pattern$B,")")

data(BCI)
nm <- nullmodel(BCI,"r2dtable")
mat1 <- as.data.frame(simulate(nm,nsim=1))
mat2 <- as.data.frame(simulate(nm,nsim=1))
colnames(mat2) <- paste0(colnames(mat2),"_B")
mat <- rbind.fill(mat1,mat2)
mat[is.na(mat)] <- 0
env_data <- c(rep("A",50),rep("B",50))
predict_df <- result_df <- cor_df_subset_long_all <- NULL

avg_dis_matrix <- as.matrix(vegdist(mat,dissim))
#C <- C_target(mat,factor=1)
#beta_pairwise<- beta_stand(mat, func = list( "beta_C"), setsize=2,args = list(C=C),summarise=F)
#avg_dis_matrix <- matrix(NA,100,100)
#avg_dis_matrix[lower.tri(avg_dis_matrix,diag=F)] <- beta_pairwise
#avg_dis_matrix[upper.tri(avg_dis_matrix,diag=F)] <- t(avg_dis_matrix)[upper.tri(avg_dis_matrix,diag=F)]

diag(avg_dis_matrix) <- NA
avg_dis <- apply(avg_dis_matrix,2,function(x) mean(x,na.rm=T))

m_observed_ad_full <- lm(avg_dis~env_data)
summary(m_observed_ad_full)
Anova(m_observed_ad_full,white.adjust=T)
result_observed_ad_full <- Anova(m_observed_ad_full,white.adjust=T)
predict_observed_ad_full <- unique(predict(m_observed_ad_full,type="response"))

LCBD_full <- LCBD.comp(vegdist(mat,dissim))$LCBD
m_observed_LCBD_full <- lm(LCBD_full~env_data)
summary(m_observed_LCBD_full)
predict_observed_LCBD_full <- unique(predict(m_observed_LCBD_full,type="response"))

CD_full <- betadisper(as.dist(avg_dis_matrix),group = c(rep("A",nrow(mat))),type="centroid",sqrt.dist=T)$distance
m_observed_CD_full <- lm(CD_full~env_data)
summary(m_observed_CD_full)
predict_observed_CD_full <- unique(predict(m_observed_CD_full,type="response"))

env_data <- data.frame(env_data=env_data)
reg_df_full <- site_dist(formula=pair_dist~env_data,mat=mat,env_data=env_data)

m_niche_full <- lm(avg_dis~env_data,data=reg_df_full)
result_niche_full <- Anova(m_niche_full,white.adjust=T)
summary(m_niche_full)
predict_niche_full <- unique(predict(m_niche_full,type="response"))

u_df_full <- data.frame(uniche=reg_df_full$avg_dis,uobserved=avg_dis,LCBD_full=LCBD_full,CD_full=CD_full)
u_df_corr_full <- cor(u_df_full)
u_df_corr_full_long <- cor_gather(u_df_corr_full)

for (j in 1:nrow(sample_pattern)) {
  for (sim in 1: 100) {
    message(sim)
    
    sub_mat_seq <- c(sample(1:(nrow(mat)/2),sample_pattern[j,1]),sample((nrow(mat)/2+1):nrow(mat),sample_pattern[j,2]))
    sub_mat <- mat[sub_mat_seq,]
    
    avg_dis_matrix <- as.matrix(vegdist(sub_mat,dissim))

    diag(avg_dis_matrix) <- NA
    avg_dis <- apply(avg_dis_matrix,2,function(x) mean(x,na.rm=T))
    
    sub_env_data <- env_data[sub_mat_seq,]
    
    m_observed <- lm(avg_dis~sub_env_data)
    summary(m_observed)
    Anova(m_observed,white.adjust=T)
    p_observed <- Anova(m_observed,white.adjust=T)
    predict_observed <- unique(predict(m_observed,type="response"))
    
    LCBD <- LCBD.comp(vegdist(sub_mat,dissim))$LCBD
    m_observed_LCBD <- lm(LCBD~sub_env_data)
    summary(m_observed_LCBD)
    predict_observed_LCBD <- unique(predict(m_observed_LCBD,type="response"))
    
    CD <- betadisper(as.dist(avg_dis_matrix),group = c(rep("A",nrow(sub_mat))),type="centroid",sqrt.dist=T)$distance
    m_observed_CD <- lm(CD~sub_env_data)
    summary(m_observed_CD)
    predict_observed_CD <- unique(predict(m_observed_CD,type="response"))
    
    sub_env_data <- data.frame(env_data=sub_env_data)
    reg_df <- site_dist(formula=pair_dist~env_data,mat=sub_mat,env_data=sub_env_data)
    
    m_niche <- lm(avg_dis~env_data,data=reg_df)
    p_niche <- Anova(m_niche,white.adjust=T)
    summary(m_niche)
    predict_niche <- unique(predict(m_niche,type="response"))
    
    u_df <- data.frame(uniche=reg_df$avg_dis,uobserved=avg_dis,LCBD=LCBD,CD=CD)
    u_cor_subset <- cor(u_df)
    u_cor_subset_long <- cor_gather(u_cor_subset)
    u_cor_subset_long$sim <- sim
    u_cor_subset_long$Effort <- sample_label[[j]]
    
    cor_df_subset_long_all <- rbind(cor_df_subset_long_all,u_cor_subset_long) 
    
    Effort <- sample_label[[j]]
    
    predict_df <- rbind(predict_df,as.data.frame(rbind(cbind(predict_observed,env=c("A","B"),"Observed",sim,Effort),
                                                       cbind(predict_niche,env=c("A","B"),"Niche",sim,Effort),
                                                       cbind(predict_observed_LCBD,c("A","B"),"Observed_LCBD",sim,Effort),
                                                       cbind(predict_observed_CD,c("A","B"),"Observed_CD",sim,Effort))))
    
    temp_df <- data.frame(p_observed = summary(m_observed)$coefficients[2,4],coef_observed=summary(m_observed)$coefficients[2,1],
                          p_niche = summary(m_niche)$coefficients[2,4],coef_niche=summary(m_niche)$coefficients[2,1],
                          p_observed_LCBD = summary(m_observed_LCBD)$coefficients[2,4],coef_observed_LCBD=summary(m_observed_LCBD)$coefficients[2,1],
                          p_observed_CD = summary(m_observed_CD)$coefficients[2,4],coef_observed_CD=summary(m_observed_CD)$coefficients[2,1],
                          sample_pattern = Effort,sim=sim)
    
    result_df <- rbind(result_df,temp_df)
    
  }
}

balanced <- cor_df_subset_long_all[cor_df_subset_long_all$Effort == "(10,10)" | cor_df_subset_long_all$Effort == "(20,20)",]
unbalanced <- cor_df_subset_long_all[cor_df_subset_long_all$Effort != "(10,10)" & cor_df_subset_long_all$Effort != "(20,20)",]

unbalanced_obs <- subset(unbalanced,var1 != "uniche" & var2 != "uniche")
unbalanced_niche <- subset(unbalanced,var1 == "uniche" | var2 == "uniche")
unbalanced_niche <- unbalanced_niche[!(unbalanced_niche$var1 == "uniche" & unbalanced_niche$var2 == "uniche"),]
unbalanced_niche <- subset(unbalanced_niche,var1 == "uniche")

cor_metric_labs <- as_labeller(c(CD = "Distance~to~centroid",
                                 LCBD = "LCBD",
                                 uobserved = "U[observed]"),
                               default = label_parsed)

cor_p <- ggplot(data=unbalanced_niche)+
  geom_violin(aes(y=cor,x=Effort))+
  facet_wrap(~var2,scale="fixed",labeller= labeller(var2=cor_metric_labs))+
  ylab("Pearson's coefficient")+
  xlab("Sampling Effort (A,B)")+
  theme_classic()
plot(cor_p)

ggsave("FigS1.tiff",width=16.8,height=8.4,units="cm",compression="lzw")
colnames(predict_df) <- c("predicted_uniqueness","Cond","Metric","sim","Effort")
predict_df$predicted_uniqueness <- as.numeric(predict_df$predicted_uniqueness)
predict_df$id <- interaction(predict_df$Cond,predict_df$Metric)

average_set <- aggregate(predicted_uniqueness~Cond+Effort+Metric,mean,data=predict_df)

long_result <- data.frame(coef=-c(result_df$coef_observed_LCBD,result_df$coef_observed_CD,result_df$coef_observed,result_df$coef_niche),
                          Effort= rep(result_df$sample_pattern,4),
                          Metric=c(rep("LCBD",400),rep("CD",400),rep("Uobserved",400),rep("Uniche",400)),
                          Sig = c(rep(c(rep("Sig",100),rep("Insig",100)),6),rep("Insig",400)))

hline_result <- data.frame(coef=-c(coefficients(m_observed_CD_full)[[2]],coefficients(m_observed_LCBD_full)[[2]],coefficients(m_observed_ad_full)[[2]],coefficients(m_niche_full)[[2]]),
                           Metric=c("CD","LCBD","Uobserved","Uniche"))

long_result$Effort <-  factor(long_result$Effort, levels = c("(10,10)", "(5,15)","(20,20)","(10,30)"))
long_result$Metric <- factor(long_result$Metric, levels=c("LCBD","CD","Uobserved","Uniche"))
hline_result$Metric <- factor(hline_result$Metric, levels=c("LCBD","CD","Uobserved","Uniche"))

long_result <- subset(long_result,Metric == "Uobserved" | Metric == "Uniche")
long_result$Metric <- droplevels(long_result$Metric)
hline_result <- subset(hline_result, Metric == "Uobserved" | Metric == "Uniche")
hline_result$Metric <- droplevels(as.factor(hline_result$Metric))
hline_result$label <- c("(a)","(b)")

metric_labs <- as_labeller(c(Uobserved = "U[observed]",
                             Uniche = "U[niche]"),
                           default = label_parsed)

p5_ef<- ggplot(long_result)+
  geom_violin(aes(y=coef,x=Effort,fill=Metric),show.legend=F)+
  geom_hline(data=hline_result,aes(yintercept=coef,group=Metric),linetype=2)+
  geom_text(data=hline_result,aes(label=label),x=-Inf,y=Inf,hjust=-0.5,vjust=1)+
  ylab("Uniqueness differences across habitats (A-B)")+
  xlab("Sampling effort (A,B)")+
  facet_wrap(~Metric,labeller= labeller(Metric=metric_labs))+
  scale_fill_discrete(values=c("blue","red"),labels=c(expression(U[observed]),expression(U[niche])))+
  theme_classic()+
  theme(legend.position="bottom")

plot(p5_ef)

ggsave("p_factor_jaccard.tiff",dpi=800,width=12.6,height=12.6,compression="lzw",units="cm")

long_result_observed <- data.frame(coef=-c(result_df$coef_observed,result_df$coef_observed_LCBD,result_df$coef_observed_CD),
                                   Effort = rep(result_df$sample_pattern,3),
                                   Metric = c(rep("Observed",400),rep("LCBD",400),rep("CD",400)))
long_result_observed$Effort <-  factor(long_result_observed$Effort, levels = c("(10,10)", "(5,15)","(20,20)","(10,30)"))

long_result_ref <- data.frame(coef=-c(coefficients(m_observed_ad_full)[[2]],coefficients(m_observed_LCBD_full)[[2]],coefficients(m_observed_CD_full)[[2]]),
                              Metric = c("Observed","LCBD","CD"))

metric_labs <- as_labeller(c(Observed="U[observed]",
                             LCBD = "LCBD",
                             CD = "Distance to centroid",default = label_parsed))


t.test(result_df$coef_observed[result_df$sample_pattern == "(20,20)"]-coefficients(m_observed_ad_full)[[2]])
t.test(result_df$coef_niche[result_df$sample_pattern == "(20,20)"]-coefficients(m_observed_ad_full)[[2]])

t.test(result_df$coef_observed[result_df$sample_pattern == "(10,10)"]-coefficients(m_observed_ad_full)[[2]])
t.test(result_df$coef_niche[result_df$sample_pattern ==  "(10,10)"]-coefficients(m_observed_ad_full)[[2]])

t.test(result_df$coef_observed[result_df$sample_pattern == "(10,30)"]-coefficients(m_observed_ad_full)[[2]])
t.test(result_df$coef_niche[result_df$sample_pattern == "(10,30)"]-coefficients(m_observed_ad_full)[[2]])

t.test(result_df$coef_observed[result_df$sample_pattern == "(5,15)"]-coefficients(m_observed_ad_full)[[2]])
t.test(result_df$coef_niche[result_df$sample_pattern == "(5,15)"]-coefficients(m_observed_ad_full)[[2]])

length(which(merged_coef_df$p_niche[merged_coef_df$Effort == "(20,20)"] < 0.05))
length(which(merged_coef_df$p_niche[merged_coef_df$Effort == "(10,10)"] < 0.05))
length(which(merged_coef_df$p_niche[merged_coef_df$Effort == "(10,30)"] < 0.05))
length(which(merged_coef_df$p_niche[merged_coef_df$Effort == "(5,15)"] < 0.05))
