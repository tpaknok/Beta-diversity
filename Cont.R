### For continuous variable simulation
library(truncnorm)
library(dplyr)
library(tidyverse)
library(plyr)
library(vegan)
library(betaC)
library(reshape2)

source("site_dist.R")
set.seed(999) 

gamma_sp <- 200 #200 species for all sites
gamma_abundance <- 100 #100 abundance for each species
gradient <- seq(0:100) #0-100 gradient. Integer
simulation <- 100 #100 simulations

sp_pool <- matrix(0,length(gradient),gamma_sp) #make a species pool matrix first
dissim = "jaccard"

for (i in 1:101) {
    sp_pool[i,i:(i+99)] <- gamma_abundance #replace species gradually along the gradient...per unit increase one species will be replaced by another (same abundance)
}

### check species pool
library(ggplot2)
longData <- melt(sp_pool)
longData[longData$value!=0,"value"] <- "100"
longData[longData$value==0,"value"] <- "0"

sp_pool_p <-ggplot(longData, aes(x = Var2, y = Var1)) + 
          geom_raster(aes(fill=value)) + 
          scale_fill_manual(values=c("white","black"))+
          labs(x="Species ID", y="Site ID (Gradient location)")+
          theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(sp_pool_p)

###

distribution <- c("gridded","biased") #gridded = evenly-spaced, biased = uneven
site_num <- 26 
local_abundance <- 100 #100 individuals in the local community (not the species pools)

result_df <- predict_df <- NULL
  for (j in distribution) { #one sampling scneario
    for (sim in 1:simulation) { # simulation number
      message(sim)
    if (j == "gridded") {
          numeric_env <- seq(0,100,by=4)
        } else {
          numeric_env <- c(sample(0:10,12,replace=T),sample(11:40,6,replace=T),sample(41:70,4,replace=T),sample(71:100,2,replace=T))
          numeric_env <- c(numeric_env,0,100)
          
          density_ggplot <- ggplot(as.data.frame(numeric_env))+
            geom_bar(aes(numeric_env))+
            scale_x_binned(n.breaks=10,labels=c("10","","30","","50","","70","","90"))+
            ylab("Count")+
            xlab("Environmental gradient")+
            theme_classic()
          
          plot(density_ggplot) #check density plot for each simulation
          
        }
        
        mat <- matrix(0,site_num,ncol(sp_pool)) #create a 26 site * 200 species matrix 
        
        for (k in 1: nrow(mat)){
          sp <- sample(1:ncol(sp_pool),local_abundance,replace=T,prob=sp_pool[numeric_env[[k]]+1,]) #sample 100 individuals from the species pool for each environmental conditons, with replacement
          comm_abundance <- count(sp)
          mat[k,comm_abundance$x] <- comm_abundance$freq
          }
        
          mat <- mat[,colSums(mat)>0] #remove species with zero abundance
          avg_dis_matrix <- as.matrix(vegdist(mat,dissim))
          diag(avg_dis_matrix) <- NA
          avg_dis <- apply(avg_dis_matrix,2,function(x) mean(x,na.rm=T)) #get Uobs
          
          m_observed <- lm(avg_dis~numeric_env+I(numeric_env^2))
          summary(m_observed)
          predict_observed <- predict(m_observed, data.frame(numeric_env=seq(0,100)),type="response")
          
          niche_df <- site_dist(pair_dist~s(numeric_env,k=3),mat=mat,env_data=as.data.frame(numeric_env)) #get Uniche
          
          m_niche <- lm(avg_dis~numeric_env+I(numeric_env^2),data=niche_df)
          summary(m_niche)
          predict_niche <- predict(m_niche,data.frame(numeric_env=seq(0,100)),type="response")
          
          predict_df <- rbind(predict_df,as.data.frame(rbind(cbind(predict_observed,numeric_env=seq(0,100),"Observed",sim,j),
                                                  cbind(predict_niche,numeric_env=seq(0,100),"Niche",sim,j))))
          
          temp_df <- data.frame(p_observed_linear = summary(m_observed)$coefficients[2,4],coef_observed_linear=summary(m_observed)$coefficients[2,1],
                                p_observed_quadratic = summary(m_observed)$coefficients[2,4],coef_observed_quadratic=summary(m_observed)$coefficients[3,1],
                                p_niche_linear = summary(m_niche)$coefficients[2,4],coef_niche_linear=summary(m_niche)$coefficients[2,1],
                                p_niche_quadratic = summary(m_observed)$coefficients[2,4],coef_niche_quadratic=summary(m_observed)$coefficients[3,1],
                                j) #summarize result
          
          result_df <- rbind(result_df,temp_df)
          

    }
  }
    
  t.test(result_df[result_df$j == "gridded","coef_observed_linear"]) #linear term estimates significant across simulations?
  t.test(result_df[result_df$j == "biased","coef_observed_linear"])
  t.test(result_df[result_df$j == "gridded","coef_observed_quadratic"])
  t.test(result_df[result_df$j == "biased","coef_observed_quadratic"])
  
  t.test(result_df[result_df$j == "gridded","coef_niche_linear"]) #non-linear term estimates significant across simulations?
  t.test(result_df[result_df$j == "biased","coef_niche_linear"])
  t.test(result_df[result_df$j == "gridded","coef_niche_quadratic"])
  t.test(result_df[result_df$j == "biased","coef_niche_quadratic"])

### Figure S1
  colnames(predict_df) <- c("predicted_uniqueness","Cond","Metric","sim","Distribution")
  predict_df$predicted_uniqueness <- as.numeric(predict_df$predicted_uniqueness)
  predict_df$Cond <- as.numeric(predict_df$Cond)
  predict_df$id <- interaction(predict_df$Distribution,predict_df$Metric,predict_df$sim)
  
  average_set <- aggregate(predicted_uniqueness~Cond+Distribution+Metric,mean,data=predict_df)
  
  Distribution_labs <- as_labeller(c(biased = "Biased",
                                     gridded = "Gridded"))
  metric_labs <- as_labeller(c(Niche="U[niche]",
                               Observed="U[observed]",default = label_parsed))

  p5_cont<- ggplot(predict_df)+
    geom_line(aes(y=predicted_uniqueness,x=Cond,group=id,colour=Metric),alpha=0.05)+
    geom_line(data=average_set,aes(y=predicted_uniqueness,x=Cond,colour=Metric))+
    ylab("Uniqueness")+
    xlab("Environmental gradient")+
    scale_color_manual(values=c("red","blue"))+
    facet_wrap(Distribution~Metric,ncol=2,scales="fixed",labeller= labeller(Metric=metric_labs,Distribution=Distribution_labs))+
    theme_classic()+
    theme(legend.position="none")
  
  plot(p5_cont)
  
  ggsave("FigS1.tiff",width=16,height=16,unit="cm",dpi=600,compression="lzw")
  