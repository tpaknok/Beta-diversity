library(betapart)
library(adespatial)
library(vegan)
library(BAT)
merged_df <- NULL
source("site_dist.R")
set.seed(100)

for (sample_size in c(20,100)) {
  site_initial <- c(20,20,sample_size)
  sp_richness <- 10
  dissim = "jaccard"
  habitat_gamma <- c(20,40,20)
  tl_gamma <- habitat_gamma[[1]]+habitat_gamma[[3]]
  link <- "log"
  generalist <- 0
  
  result_df <- list()
  
  env_data <- c(rep("A",site_initial[[1]]),rep("B",site_initial[[2]]),rep("C",site_initial[[3]]))

  for (i in 1:100) {
    message("i = ",i)
    site_mat_list <- list()
    
    mat <- matrix(0,sum(site_initial),tl_gamma)
    
    for (j in 1:nrow(mat)) {
      if (j <= site_initial[[1]]) {
        mat[j,sample(1:habitat_gamma[[1]],sp_richness)] <- 1
      }  else if (j > sum(site_initial)-site_initial[[3]]){
        mat[j,sample((tl_gamma/2+1):tl_gamma,sp_richness)] <- 1
      } else {
        mat[j,sample(1:habitat_gamma[[2]],sp_richness)] <- 1
      }
    }
    
    avg_dis_matrix <- as.matrix(vegdist(mat,dissim))
    diag(avg_dis_matrix) <- NA
    avg_dis <- apply(avg_dis_matrix,2,function(x) mean(x,na.rm=T))

    env_dist <- as.matrix(gower(as.matrix(env_data)))
    env_dist[col(env_dist)==row(env_dist)] <- NA
    
    m_observed_ad <- lm(avg_dis~env_data)
    summary(m_observed_ad)
    predict_observed_ad <- unique(predict(m_observed_ad,type="response"))

    LCBD <- LCBD.comp(vegdist(mat,dissim))$LCBD
    m_observed_LCBD <- lm(LCBD~env_data)
    summary(m_observed_LCBD)
    predict_observed_LCBD <- unique(predict(m_observed_LCBD,type="response"))
    
    CD <- LCBD.comp(vegdist(mat,dissim))$LCBD*LCBD.comp(vegdist(mat,dissim))$beta[2]
    m_observed_CD <- lm(CD~env_data)
    summary(m_observed_CD)
    predict_observed_CD <- unique(predict(m_observed_CD,type="response"))
    
    w <- 1/table(env_data)
    w <- w/min(w)
    w <- c(rep(w[[1]],site_initial[[1]]),rep(w[[2]],site_initial[[2]]),rep(w[[3]],site_initial[[3]]))
    
    w_avg_dis <- apply(avg_dis_matrix,2,function(x) weighted.mean(x,w=w,na.rm=T))
    m_weighted_ad <- lm(w_avg_dis~env_data)
    summary(m_weighted_ad)
    predict_weighted_ad <- unique(predict(m_weighted_ad,type="response"))
    
    reg_df <- site_dist(mat,env_data)
    
    within_beta <- subset(reg_df,type=="Within")
    m_within <- lm(avg_dis~env_data,data=within_beta)
    summary(m_within)
    predict_within <- unique(predict(m_within,type="response"))
    
    between_beta <- subset(reg_df,type=="Between")
    m_between <- lm(avg_dis~env_data,data=between_beta)
    summary(m_between)
    predict_between <- unique(predict(m_between,type="response"))
    
    total_beta <- subset(reg_df,type=="Total")
    m_total <- lm(avg_dis~env_data,data=total_beta)
    summary(m_total)
    predict_total <- unique(predict(m_total,type="response"))
    
    #weight_LCBD <- c(rep(site_initial[[1]],site_initial[[1]]),rep(site_initial[[2]],site_initial[[2]]),rep(site_initial[[3]],site_initial[[3]]))
    
    #w_LCBD <- W_LCBD(mat,env,weight_LCBD,adjust.bias=F)
    #m_weighted_LCBD <- lm(w_LCBD~numeric_env+I(numeric_env^2))
    #summary(m_weighted_LCBD)
    #predict_weighted_LCBD <- unique(predict(m_weighted_LCBD,type="response"))
    
    Effort <- ifelse(length(unique(site_initial)) == 1, "Balanced","Unbalanced")
    
    result_df[[i]] <- as.data.frame(rbind(cbind(predict_observed_ad,c("A","B","C"),"Observed_ad",i),
                                          cbind(predict_observed_LCBD,c("A","B","C"),"Observed_LCBD",i),
                                          cbind(predict_observed_CD,c("A","B","C"),"Observed_CD",i),
                                          cbind(predict_within,c("A","B","C"),"Beta_Within",i),
                                          cbind(predict_between,c("A","B","C"),"Beta_Between",i),
                                          cbind(predict_total,c("A","B","C"),"Beta_Total",i)))
    
    
  }
  
  result_df_ordinal <- do.call(rbind,result_df)
  colnames(result_df_ordinal) <- c("predicted_beta_site","Habitat","Metric","i")
  result_df_ordinal$predicted_beta_site <- as.numeric(result_df_ordinal$predicted_beta_site)
  result_df_ordinal$id <- interaction(result_df_ordinal$Habitat,result_df_ordinal$Metric)
  result_df_ordinal$type <- ifelse(length(unique(site_initial)) == 1, "Balanced","Unbalanced")
  merged_df <- rbind(merged_df,result_df_ordinal)
}

library(ggplot2)
library(ggh4x)
#merged_df <- subset(merged_df,Metric != "Beta_Total")
metric_labs <- as_labeller(c(Beta_Between="beta[between]",
                             Beta_Within="beta[within]",
                             Beta_Total="beta[total]",
                             Observed_ad="Mean~beta[pair]",
                             Observed_CD="Distance~to~centroid",
                             Observed_LCBD= "LCBD"),default = label_parsed)

p5_factor<- ggplot(merged_df)+
  geom_boxplot(aes(y=predicted_beta_site,x=Habitat,fill=Habitat))+
  ylab("Value")+
  xlab("Habitat")+
  #ylim(0.75,1)+
  scale_fill_discrete(labels= c("A" = expression("A ("~gamma~"=20 )"),
                                "B" = expression("B ("~gamma~"=40 )"),
                                "C" = expression("C ("~gamma~"=20 )")))+
  theme(legend.position="bottom")+
  facet_wrap(type~Metric,scales="free",ncol=5,labeller= labeller(Metric=metric_labs))+
  theme_classic()+
  theme(legend.position = "bottom")

plot(p5_factor)
ggsave("revised_p1.tiff",width=10,heigh=5,dpi=800,compression="lzw")
###
habitat_gamma1 <- 20
habitat_gamma2 <- 20
site_num <- 99+2
sp_richness <- 10
distribution <- c("skewed_normal","uniform")
merged_df <- merged_df_env <- NULL
simulation <- 100
dissim <- "gower"
set.seed(999)

library(fGarch)
library(truncnorm)

  for (j in distribution) {
    
    if (j == "uniform") {
    numeric_env <- seq(0,1,0.01)
    } else {
      numeric_env <- rtruncnorm(site_num-2, a=0, b=1, mean = 0.3, sd = 0.3)
      numeric_env <- c(numeric_env,0,1)
      #numeric_env <- ifelse(numeric_env < 0, 0, ifelse(numeric_env > 1,1,numeric_env))
      plot(density(numeric_env))
    }
  
    result_df <- NULL
    density_df <- NULL
    for (sim in 1:simulation) {
      message(sim)
      
      mat <- as.data.frame(matrix(0,site_num,habitat_gamma1+habitat_gamma2))
      
      for (i in 1:nrow(mat)) {
      p1 <- rep(numeric_env[[i]],habitat_gamma1)
      p2 <- rep(1-numeric_env[[i]],habitat_gamma2)
      
      p <- c(p1,p2)
      mat[i,sample(1:40,sp_richness,prob=p)] <- 1
      }
    
      avg_dis_matrix <- as.matrix(vegdist(mat,dissim))
      diag(avg_dis_matrix) <- NA
      avg_dis <- apply(avg_dis_matrix,2,function(x) mean(x,na.rm=T))
      
      predict_num_env <- data.frame(numeric_env=numeric_env)
      m_observed_ad <- lm(avg_dis~numeric_env+I(numeric_env^2))
      summary(m_observed_ad)
      predict_observed_ad <- predict(m_observed_ad, predict_num_env,type="response")
      
      LCBD <- LCBD.comp(vegdist(mat,dissim))$LCBD
      m_observed_LCBD <- lm(LCBD~numeric_env+I(numeric_env^2))
      summary(m_observed_LCBD)
      predict_observed_LCBD <- unique(predict(m_observed_LCBD,predict_num_env,type="response"))
      
      CD <- LCBD.comp(vegdist(mat,dissim))$LCBD*LCBD.comp(vegdist(mat,dissim))$beta[2]
      m_observed_CD <- lm(CD~numeric_env+I(numeric_env^2))
      summary(m_observed_CD)
      predict_observed_CD <- unique(predict(m_observed_CD,predict_num_env,type="response"))
      
      reg_df <- site_dist(mat,numeric_env,splines=10)
      
      reg_df$env <- as.numeric(reg_df$env)
      within_beta <- subset(reg_df,type=="Within")
      m_within <- lm(avg_dis~env+I(env^2),data=within_beta)
      summary(m_within)
      predict_within <- predict(m_within,data.frame(env=unique(reg_df$env)),type="response")
      
      between_beta <- subset(reg_df,type=="Between")
      m_between <- lm(avg_dis~env+I(env^2),data=between_beta)
      summary(m_between)
      predict_between <- predict(m_between,data.frame(env=unique(reg_df$env)),type="response")
      
      total_beta <- subset(reg_df,type=="Total")
      m_total <- lm(avg_dis~env+I(env^2),data=total_beta)
      summary(m_total)
      predict_total <- predict(m_total,data.frame(env=unique(reg_df$env)),type="response")
      
      result_df[[sim]] <- as.data.frame(rbind(cbind(predict_observed_ad,unique(reg_df$env),"Observed_ad",sim,j),
                                            cbind(predict_observed_LCBD,unique(reg_df$env),"Observed_LCBD",sim,j),
                                            cbind(predict_observed_CD,unique(reg_df$env),"Observed_CD",sim,j),
                                            cbind(predict_within,unique(reg_df$env),"Beta_Within",sim,j),
                                            cbind(predict_between,unique(reg_df$env),"Beta_Between",sim,j),
                                            cbind(predict_total,unique(reg_df$env),"Beta_Total",sim,j)))
      
      density_df[[sim]] <- data.frame(env=numeric_env,sim=sim,distribution=j)
    }
    
    result_df_cont <- do.call(rbind,result_df)
    colnames(result_df_cont) <- c("predicted_beta_site","Cond","Metric","sim","distribution")
    result_df_cont$predicted_beta_site <- as.numeric(result_df_cont$predicted_beta_site)
    result_df_cont$Cond <- as.numeric(result_df_cont$Cond)
    
    result_df_cont$id <- interaction(result_df_cont$distribution,result_df_cont$Metric,result_df_cont$sim)
    merged_df <- rbind(merged_df,result_df_cont)
    merged_df_env <- rbind(merged_df_env,do.call(rbind,density_df))
    }
  
  #merged_df <- subset(merged_df,Metric != "Beta_Total")
  average_set <- aggregate(predicted_beta_site~Cond+distribution+Metric,mean,data=merged_df)
  
  library(ggplot2)
  library(ggh4x)
  #merged_df$Metric <- factor(merged_df$Metric,levels=c("Observed_ad","Observed_CD","Observed_LCBD","Within","Between"))
  metric_labs <- as_labeller(c(Beta_Between="beta[between]",
                               Beta_Within="beta[within]",
                               Beta_Total="beta[Total]",
                               Observed_ad="Mean~beta[pair]",
                               Observed_CD="Distance~to~centroid",
                               Observed_LCBD= "LCBD"),default = label_parsed,)
  distribution_labs <- as_labeller(c(skewed_normal = "Skewed normal",
                                     uniform = "Unifrom"))
  
  p5_cont<- ggplot(merged_df)+
    geom_line(aes(y=predicted_beta_site,x=Cond,group=id),alpha=0.05)+
    geom_line(data=average_set,aes(y=predicted_beta_site,x=Cond,group=distribution),alpha=1)+
    ylab("Value")+
    xlab("Env")+
    theme(legend.position="bottom")+
    facet_wrap(~Metric,scales="free",ncol=6,labeller= labeller(Metric=metric_labs,distribution=distribution_labs))+
    theme_classic()
  
  plot(p5_cont)
  
  density_plot_df <- data.frame(env=rep(seq(0,1,0.01),2),distribution=c(rep("truncated-normal",101),rep("uniform",101)))
  density_plot_df$y <- 1
  density_plot_df$y[density_plot_df$distribution == "truncated-normal"]  <- dtruncnorm(seq(0,1,0.01), a=0, b=1, mean = 0.3, sd = 0.3)
  
  density_ggplot <- ggplot(density_plot_df)+
    geom_line(aes(y=y,x=env,group=distribution))+
    xlim(0,1)+
    ylab("Density")+
    xlab("Env")+
    facet_wrap(~distribution,nrow=2,scale="free")+
    theme_classic()
  
  plot(density_ggplot)  
  
  library(ggpubr)
  overall <- ggarrange(density_ggplot,p5_cont,widths=c(1/6,1))
  plot(overall)  
  
  ggsave("revised_p2.tiff",width=12,heigh=6,dpi=800,compression="lzw")