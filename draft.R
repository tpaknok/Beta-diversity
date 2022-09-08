#####
library(plyr)
library(BAT)
library(adespatial)
library(vegan)
library(glmmTMB)
library(ggeffects)
source("site_dist.R")
random_sample <- seq(6,20,2)
n_h1 <- 20
n_h2 <- 40
site_h1 <- 20
site_h2 <- 20
SR <- 5
dissim <- "jaccard"
env_dissim <- "gower"
resample_repeat <- 1
result_df <- temp_result_df<-NULL

set.seed(999)
for (j in 1:100) {
  message(j)
  h1_mat <- matrix(0,nrow=site_h1,ncol=n_h1)
  for (i in 1:site_h1){
    seq <- sample(1:n_h1,SR)
    h1_mat[i,seq] <- 1
  }
  h2_mat <- matrix(0,nrow=site_h2,ncol=n_h2)
  for (i in 1:site_h2){
    seq <- sample(1:n_h2,SR)
    h2_mat[i,seq] <- 1
  }
  
  colnames(h2_mat) <- paste0("sp",1:n_h2,"h2")
  
  h1_mat <- h1_mat[,colSums(h1_mat) > 0]
  h2_mat <- h2_mat[,colSums(h2_mat) > 0]
  
  mat <- rbind.fill(as.data.frame(h1_mat),as.data.frame(h2_mat))
  mat[is.na(mat)] <- 0
  group <- c(rep("Group1",site_h1),rep("Group2",site_h2))
  
  cent_dist_all  <- betadisper(vegdist(mat,"bray"),group=rep("a",nrow(mat)))
  cent_dist_all  <- cent_dist_all$distances
  cent_dist_all <- (cent_dist_all*(length(cent_dist_all)-1)+0.5)/length(cent_dist_all)
  cent_dist_within <- betadisper(vegdist(mat,"bray"),group=group)
  cent_dist_within <- cent_dist_within$distances
  cent_dist_within <- (cent_dist_within*(length(cent_dist_within)-1)+0.5)/length(cent_dist_within)
  
  m_balanced_all <- glmmTMB(cent_dist_all~group,family=beta_family())
  summary(m_balanced_all)
  m_balanced_within <- glmmTMB(cent_dist_within~group,family=beta_family())
  summary(m_balanced_within)
  
  pairwise_dist <- as.matrix(vegdist(mat,dissim))
  pairwise_dist[col(pairwise_dist)==row(pairwise_dist)] <- NA
  
  avg_dis_all <- colMeans(pairwise_dist,na.rm=T)
  avg_dis_all <- (avg_dis_all*(length(avg_dis_all)-1)+0.5)/length(avg_dis_all)
  avg_dis_within <- apply(pairwise_dist,2,function(x) mean(x[group== group[which(is.na(x))]],na.rm=T))
  avg_dis_within <- (avg_dis_within*(length(avg_dis_within)-1)+0.5)/length(avg_dis_within)
  
  m_balanced_all_ad <- glmmTMB(avg_dis_all~group,family=beta_family())
  summary(m_balanced_all_ad)
  m_balanced_within_ad <- glmmTMB(avg_dis_within~group,family=beta_family())
  summary(m_balanced_within_ad)
  
  for (sample in random_sample) {
    message(sample)
    temp_result_df <- NULL
    
    for (m in 1:resample_repeat) {  
      message("m = ",m)
      resample <- c(sample(1:20,sample),sample(21:40,20))
      subset_mat <- mat[resample,]
      subset_mat <- subset_mat[,colSums(subset_mat) > 0]
      subset_group <- group[resample]
      
      cent_dist_all  <- betadisper(vegdist(subset_mat,"bray"),group=rep("a",nrow(subset_mat)))
      cent_dist_all  <- cent_dist_all$distances
      cent_dist_all <- (cent_dist_all*(length(cent_dist_all)-1)+0.5)/length(cent_dist_all)
      cent_dist_within <- betadisper(vegdist(subset_mat,"bray"),group=subset_group)
      cent_dist_within <- cent_dist_within$distances
      cent_dist_within <- (cent_dist_within*(length(cent_dist_within)-1)+0.5)/length(cent_dist_within)
      
      m_unbalanced_all <- glmmTMB(cent_dist_all~subset_group,family=beta_family())
      summary(m_unbalanced_all)
      m_unbalanced_within <- glmmTMB(cent_dist_within~subset_group,family=beta_family())
      summary(m_unbalanced_within)
      
      pairwise_dist <- as.matrix(vegdist(subset_mat,dissim))
      pairwise_dist[col(pairwise_dist)==row(pairwise_dist)] <- NA
      
      avg_dis_all <- colMeans(pairwise_dist,na.rm=T)
      avg_dis_all <- (avg_dis_all*(length(avg_dis_all)-1)+0.5)/length(avg_dis_all)
      avg_dis_within <- apply(pairwise_dist,2,function(x) mean(x[subset_group== subset_group[which(is.na(x))]],na.rm=T))
      avg_dis_within <- (avg_dis_within*(length(avg_dis_within)-1)+0.5)/length(avg_dis_within)
      avg_dis_controlled_0 <- site_dist(subset_mat,subset_group,dissim,env_dissim,0,quadratic.term=F,trans=T)
      
      m_unbalanced_all_ad <- glmmTMB(avg_dis_all~subset_group,family=beta_family())
      summary(m_unbalanced_all_ad)
      m_unbalanced_within_ad <- glmmTMB(avg_dis_within~subset_group,family=beta_family())
      summary(m_unbalanced_within_ad)
      m_unbalanced_controlled_ad0 <- glmmTMB(avg_dis~env,family=beta_family(),data=avg_dis_controlled_0)
      summary(m_unbalanced_controlled_ad0)
      
      avg_dis_controlled_1 <- site_dist(subset_mat,subset_group,dissim,env_dissim,1,quadratic.term=F,trans=T)
      m_unbalanced_controlled_ad1 <- glmmTMB(avg_dis~env,family=beta_family(),data=avg_dis_controlled_1)
      summary(m_unbalanced_controlled_ad1)
      
      avg_dis_controlled_0.5 <- site_dist(subset_mat,subset_group,dissim,env_dissim,0.5,quadratic.term=F,trans=T)
      m_unbalanced_controlled_ad_0.5 <- glmmTMB(avg_dis~env,family=beta_family(),data=avg_dis_controlled_0.5)
      summary(m_unbalanced_controlled_ad_0.5)
      
      temp_result <- c(summary(m_balanced_all)$coefficients$cond[1,1],
                       summary(m_balanced_all)$coefficients$cond[2,c(1,4)],
                       summary(m_balanced_within)$coefficients$cond[1,1],
                       summary(m_balanced_within)$coefficients$cond[2,c(1,4)],
                       summary(m_balanced_all_ad)$coefficients$cond[1,1],
                       summary(m_balanced_all_ad)$coefficients$cond[2,c(1,4)],
                       summary(m_balanced_within_ad)$coefficients$cond[1,1],
                       summary(m_balanced_within_ad)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_all)$coefficients$cond[1,1],
                       summary(m_unbalanced_all)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_within)$coefficients$cond[1,1],
                       summary(m_unbalanced_within)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_all_ad)$coefficients$cond[1,1],
                       summary(m_unbalanced_all_ad)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_within_ad)$coefficients$cond[1,1],
                       summary(m_unbalanced_within_ad)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_controlled_ad0)$coefficients$cond[1,1],
                       summary(m_unbalanced_controlled_ad0)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_controlled_ad1)$coefficients$cond[1,1],
                       summary(m_unbalanced_controlled_ad1)$coefficients$cond[2,c(1,4)],
                       summary(m_unbalanced_controlled_ad_0.5)$coefficients$cond[1,1],
                       summary(m_unbalanced_controlled_ad_0.5)$coefficients$cond[2,c(1,4)])
      
      temp_result <- as.data.frame(t(temp_result))
      temp_result$sample <- sample
      temp_result$resample_repeat <- m
      temp_result$matrix <- j
      temp_result_df <- rbind(temp_result_df,temp_result)
    }
    result_df <- rbind(result_df,temp_result_df)
    
    #t.test(temp_result_df$balanced_eff,temp_result_df$env_controlled_eff,paired = T)
    #t.test(temp_result_df$result_df[,2],temp_result_df$result_df[,8],paired = T)
    #t.test(temp_result_df$result_df[,5],temp_result_df$result_df[,8],paired = T)
  }
}

library(ggplot2)
plot_df <- data.frame(Effect = c(result_df[,14],result_df[,17]), Group = c(rep("All",nrow(result_df)),rep("Within",nrow(result_df))),sample=rep(result_df$sample,2))
plot_df$id <- interaction(plot_df$Group,plot_df$sample)

p_dist <- ggplot(plot_df)+
  geom_boxplot(aes(y=Effect,x=sample,group=id,fill=Group),position = position_dodge(preserve = "single"))+
  ylab("Effects (B-A) based on distance to centroid")+
  xlab("Number of assemblages of Habitat A")+
  theme_classic()

plot(p_dist)

plot_df_ad <- data.frame(Effect = c(result_df[,20],result_df[,23],result_df[,26],result_df[,29],result_df[,32]), Group = c(rep("All",nrow(result_df)),rep("Within",nrow(result_df)),rep("ED = 0",nrow(result_df)),rep("ED = 1",nrow(result_df)),rep("ED = 0.5",nrow(result_df))),sample=rep(result_df$sample,5))
plot_df_ad$id <- interaction(plot_df_ad$Group,plot_df_ad$sample)

p_ad2 <- ggplot(plot_df_ad)+
  geom_boxplot(aes(y=Effect,x=sample,group=id,fill=Group),position = position_dodge(preserve = "single"))+
  ylab("Effects (B-A) based on average pairwise dissimilarity")+
  xlab("Number of assemblages of Habitat A")+
  theme_classic()

plot(p_ad2)

library(ggpubr)
p <- ggarrange(p_dist,p_ad,common.legend = T, legend="bottom")
plot(p)

ggsave("p1.tiff",dpi=800,width=16,height=8,compression="lzw")
  reg_df <- site_dist(subset_mat,subset_group,"bray","gower",0)
  m_env_controlled_ad <- lm(avg_dis~subset_group,data=reg_df)
  
###
  
  site_initial <- c(20,20,20)
  sp_richness <- 5
  dissim = "jaccard"
  
  result_df <- list()
  set.seed(999)
  for (i in 1:100) {
    message("i = ",i)
    site_mat_list <- list()
    
    mat <- matrix(0,sum(site_initial),40)
    
    for (j in 1:nrow(mat)) {
      if (j <= 20) {
        mat[j,sample(1:20,sp_richness)] <- 1
      }  else if (j >= 41){
        mat[j,sample(21:40,sp_richness)] <- 1
      } else {
        mat[j,sample(1:40,sp_richness)] <- 1
      }
    }
    
    env_data <- c(rep("A",20),rep("B",20),rep("C",20))
    env_data <- as.factor(env_data)
    
    avg_dis <- as.matrix(vegdist(mat,dissim))
    avg_dis[col(avg_dis)==row(avg_dis)] <- NA
    avg_dis <- colMeans(avg_dis,na.rm=T)
    
    env_dist <- as.matrix(gower(as.matrix(env_data)))
    env_dist[col(env_dist)==row(env_dist)] <- NA
    env_dist <- colMeans(env_dist,na.rm=T)
    
    plot(env_data,avg_dis)
    m_balanced_ad <- glmmTMB(avg_dis~env_data,family=beta_family())
    summary(m_balanced_ad)
    predict_observed <- predict(m_balanced_ad,newdata = data.frame(env_data=c("A","B","C")),type="response")
    
    reg_df <- site_dist(mat,env_data,"bray","gower",min(apply(as.matrix(gower(as.matrix(env_data))),2,max)),quadratic.term=F,trans = T)
    m_env_controlled_ad_max <- glmmTMB(avg_dis~env,data=reg_df,family=beta_family())
    summary(m_env_controlled_ad_max)
    plot(reg_df$env,reg_df$avg_dis)
    predict_max <- predict(m_env_controlled_ad_max,newdata = data.frame(env=c("A","B","C")),type="response")
    
    reg_df <- site_dist(mat,env_data,"bray","gower",0,quadratic.term=F,trans = T)
    m_env_controlled_ad_0 <- glmmTMB(avg_dis~env,data=reg_df,family=beta_family())
    summary(m_env_controlled_ad_0)
    plot(reg_df$env,reg_df$avg_dis)
    predict_0 <- predict(m_env_controlled_ad_0,newdata = data.frame(env=c("A","B","C")),type="response")
    
    result_df[[i]] <- as.data.frame(rbind(cbind(predict_observed,c("A","B","C"),"Observed",i),
                                    cbind(predict_max,c("A","B","C"),"1",i),
                                    cbind(predict_0,c("A","B","C"),"0",i)))
  }

  result_df_factor <- do.call(rbind,result_df)
  colnames(result_df_factor) <- c("predicted_beta_site","Habitat","ED_cond","i")
  result_df_factor$predicted_beta_site <- as.numeric(result_df_factor$predicted_beta_site)
  result_df_factor$id <- interaction(result_df_factor$Habitat,result_df_factor$ED_cond)
  
  p4 <- ggplot(result_df_factor)+
    geom_boxplot(aes(y=predicted_beta_site,x=ED_cond,group=id,fill=Habitat),position = position_dodge(preserve = "single"))+
    ylab(~ paste("Predicted ",beta[Site]))+
    xlab("ED value")+
    scale_fill_discrete(labels= c(A = expression("A ("~gamma~"=20 )"),
                        B = expression("B ("~gamma~"=40 )"),
                        C = expression("C ("~gamma~"=20 )")))+
    scale_x_discrete(labels= c(Max = 0.81))+
    theme_classic()+
    theme(legend.position="bottom")
  
  plot(p4)
  
  ggsave("p4.tiff",dpi=800,compression="lzw")
###
  
  site_initial <- c(20,20,10)
  sp_richness <- 5
  dissim = "jaccard"
  set.seed(999)
  gamma <- 40
  
  result_df <- list()
  for (i in 1:100) {
    message("i = ",i)
    site_mat_list <- list()
    
    mat <- matrix(0,sum(site_initial),gamma)
    
    for (j in 1:nrow(mat)) {
      if (j <= site_initial[[1]]) {
        mat[j,sample(1:(gamma/2),sp_richness)] <- 1
      }  else if (j > sum(site_initial)-site_initial[[3]]){
        mat[j,sample((gamma/2+1):gamma,sp_richness)] <- 1
      } else {
        mat[j,sample(1:gamma,sp_richness)] <- 1
      }
    }
    
    env_data <- c(rep(1,site_initial[[1]]),rep(2,site_initial[[2]]),rep(3,site_initial[[3]]))
    
    avg_dis <- as.matrix(vegdist(mat,dissim))
    avg_dis[col(avg_dis)==row(avg_dis)] <- NA
    avg_dis <- colMeans(avg_dis,na.rm=T)
    
    env_dist <- as.matrix(gower(as.matrix(env_data)))
    env_dist[col(env_dist)==row(env_dist)] <- NA
    env_dist <- colMeans(env_dist,na.rm=T)
    
    plot(env_data,avg_dis)
    m_balanced_ad <- glmmTMB(avg_dis~env_data+I(env_data^2),family=beta_family())
    summary(m_balanced_ad)
    predict_observed <- predict(m_balanced_ad,newdata = data.frame(env_data=c(1,2,3)),type="response")
    
    reg_df <- site_dist(mat,env_data,dissim,"gower",min(apply(gower(as.matrix(env_data)),2,max)),quadratic.term=F,trans = T)
    m_env_controlled_ad_max <- glmmTMB(avg_dis~env+I(env^2),data=reg_df,family=beta_family())
    summary(m_env_controlled_ad_max)
    plot(reg_df$env,reg_df$avg_dis)
    predict_max <- predict(m_env_controlled_ad_max,newdata = data.frame(env=c(1,2,3)),type="response")
    
    reg_df <- site_dist(mat,env_data,dissim,"gower",0,quadratic.term=F,trans = T)
    m_env_controlled_ad_0 <- glmmTMB(avg_dis~env+I(env^2),data=reg_df,family=beta_family())
    summary(m_env_controlled_ad_0)
    plot(reg_df$env,reg_df$avg_dis)
    predict_0 <- predict(m_env_controlled_ad_0,newdata = data.frame(env=c(1,2,3)),type="response")
    
    reg_df <- site_dist(mat,env_data,dissim,"gower",apply(gower(as.matrix(env_data)),2,max),quadratic.term=F,trans = T)
    #reg_df <- site_dist(mat,env_data,dissim,"gower",env_dist,quadratic.term=F,trans = T)
    m_uncontrolled_ED <- glmmTMB(avg_dis~env+I(env^2),data=reg_df,family=beta_family())
    summary(m_uncontrolled_ED)
    plot(reg_df$env,reg_df$avg_dis)
    predict_uncontrolled_ED <- predict(m_uncontrolled_ED ,newdata = data.frame(env=c(1,2,3)),type="response")
    
    result_df[[i]] <- as.data.frame(rbind(cbind(predict_observed,c(1,2,3),"Observed",i),
                                          cbind(predict_max,c(1,2,3),"Max",i),
                                          cbind(predict_0,c(1,2,3),"0",i),
                                    cbind(predict_uncontrolled_ED,c(1,2,3),"Habitat_Max",i)))
    }
  
  result_df_ordinal <- do.call(rbind,result_df)
  colnames(result_df_ordinal) <- c("predicted_beta_site","Habitat","ED_cond","i")
  result_df_ordinal$predicted_beta_site <- as.numeric(result_df_ordinal$predicted_beta_site)
  result_df_ordinal$id <- interaction(result_df_ordinal$Habitat,result_df_ordinal$ED_cond)
  
  library(ggplot2)
  p52 <- ggplot(result_df_ordinal)+
    geom_boxplot(aes(y=predicted_beta_site,x=ED_cond,group=id,fill=Habitat),position = position_dodge(preserve = "single"))+
    ylab(~ paste("Predicted ",beta[Site]))+
    xlab("ED value")+
    ylim(0.80,1)+
    scale_fill_discrete(labels= c("1" = expression("1 ("~gamma~"=20 )"),
                                  "2" = expression("2 ("~gamma~"=40 )"),
                                  "3" = expression("3 ("~gamma~"=20 )")))+
    scale_x_discrete(labels= c(Max = 0.5))+
    theme_classic()+
    theme(legend.position="bottom")
  
  plot(p52)
  
  plot(p5)
  
  library(ggpubr)
  p_overall <- ggarrange(p5,p52)
  ggsave("p5.tiff",dpi=800,compression="lzw")
  