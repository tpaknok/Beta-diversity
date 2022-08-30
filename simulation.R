library(CommEcol)
library(plyr)
library(vegan)
coo<- c(seq(10,100,10))


g1 <- compas(S=50, dims=1, am=1, beta.R=1, clump=1,coords=coo, 
       n.quanti=0.01,  n.quali=0.01, add1=0.05)

g2 <- compas(S=50, dims=1, am=1, beta.R=1, clump=1,coords=coo, 
       n.quanti=0.01,  n.quali=0.01, add1=0.05)

library(BAT)
pa_g1 <- decostand(g1,"pa")
colnames(pa_g1) <- paste0(colnames(pa_g1),"_g1")
pa_g2 <- decostand(g2,"pa")
colnames(pa_g2) <- paste0(colnames(pa_g2),"_g2")

binded_df <- rbind.fill(as.data.frame(pa_g1),as.data.frame(pa_g2))
binded_df$group <- c(rep("Group1",10),rep("Group2",10))
binded_df[is.na(binded_df)] <- 0

LCBD_balanced <- LCBD.comp(vegdist(binded_df[,colnames(binded_df) != "group"],"jaccard"))$LCBD
m_balanced <- lm(LCBD_balanced~group,data=binded_df)
summary(m_balanced)

binded_df <- binded_df[c(sample(1:10,3),sample(11:20,10)),]
LCBD <- LCBD.comp(vegdist(binded_df[,colnames(binded_df) != "group"],"jaccard"))
LCBD_df <- data.frame(LCBD=LCBD$LCBD,group=binded_df$group)

LCBD_env <- LCBD.comp(gower(binded_df$group))$LCBD

LCBD_null <- NULL
coef_null <- NULL
for (i in 1:100){
  message(i)
  nm1 <- nullmodel(binded_df[binded_df$group == "Group1",colnames(binded_df) != "group"],"quasiswap")
  sm1 <- as.data.frame(simulate(nm1,1))
  nm2 <- nullmodel(binded_df[binded_df$group == "Group2",colnames(binded_df) != "group"],"quasiswap")
  sm2 <- as.data.frame(simulate(nm2,1))

  sm <- rbind(sm1,sm2)
  sm <- sm[,colSums(sm) > 0]
  #nm <- nullmodel(binded_df[,colnames(binded_df) != "group"],"c0")
  #sm <- as.data.frame(simulate(nm,1))
  
  LCBD_null <- cbind(LCBD_null,LCBD.comp(vegdist(sm,"jaccard"))$LCBD)
  temp_null <- LCBD.comp(vegdist(sm,"jaccard"))$LCBD
  
  m_iter<-glmmTMB(temp_null~group,data=binded_df)
  summary(m_iter)
  coef_null <- c(coef_null,summary(m_iter)$coefficients$cond[2,1])
}


library(car)
m_unbalanced<-glmmTMB(LCBD~group,data=LCBD_df)
summary(m_unbalanced)
summary(m_unbalanced)$coefficients$cond[2,1] > coef_null
Anova(m,white.adjust=T)
predict(m,newdata=data.frame(group=c("Group1","Group2")),type="response")


###
m2 <- glm(binded_df$sp.1_g1~group,data=binded_df,family="binomial")
predict(m2,type="response")

#####
library(plyr)
library(BAT)
library(adespatial)
library(vegan)
library(CommEcol)
source("C:/Users/pakno/OneDrive/Desktop/LCBD/site_dist.R")
random_sample <- seq(6,20,2)
n_h1 <- 20
n_h2 <- 20
site_h1 <- 20
site_h2 <- 20
SR <- 5
dissim <- "bray"
resample_repeat <- 1
result_df <- temp_result_df<-NULL

  for (j in 1:100) {
    message(j)
  h1_mat <- matrix(0,nrow=site_h1,ncol=n_h1)
  for (i in 1:site_h1){
    seq <- sample(1:n_h1,SR)
    h1_mat[i,seq] <- 1
  }
  
  colnames(h1_mat) <- paste0("sp",1:n_h1,"h1")
  
  
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
  
  LCBD_balanced <- LCBD.comp(vegdist(mat,dissim),sqrt.D=T)$LCBD*LCBD.comp(vegdist(mat,dissim),sqrt.D=T)$beta[[1]]
  avg_dis <- as.matrix(vegdist(mat,dissim))
  avg_dis[col(avg_dis)==row(avg_dis)] <- NA
  avg_dis <- colMeans(avg_dis,na.rm=T)
  
  #LCBD_balanced <- LCBD.comp(dis.nness(mat),sqrt.D=T)$LCBD
  #LCBD_balanced <- LCBD.comp(ESS(mat))$LCBD
  m_balanced <- lm(LCBD_balanced~group)
  summary(m_balanced)
  m_balanced_ad <- lm(avg_dis~group)
  summary(m_balanced_ad)
  
  for (sample in random_sample) {
  message(sample)
  temp_result_df <- NULL
  
  for (m in 1:resample_repeat) {  
    message("m = ",m)
  resample <- c(sample(1:20,sample),sample(21:40,20))
  subset_mat <- mat[resample,]
  subset_mat <- subset_mat[,colSums(subset_mat) > 0]
  subset_group <- group[resample]
  LCBD_env <- LCBD.comp(gower(subset_group),sqrt.D=T)$LCBD
  
  LCBD <- LCBD.comp(vegdist(subset_mat[order(as.numeric(rownames(subset_mat))),],dissim),sqrt.D=T)$LCBD*LCBD.comp(vegdist(mat,dissim),sqrt.D=T)$beta[[1]]
  avg_dis <- as.matrix(vegdist(subset_mat,dissim))
  avg_dis[col(avg_dis)==row(avg_dis)] <- NA
  
  #LCBD <- LCBD.comp(dis.nness(subset_mat))
  #LCBD <- LCBD.comp(ESS(subset_mat))
  
  LCBD_df <- data.frame(LCBD=LCBD,avg_dis=colMeans(avg_dis,na.rm=T),group=subset_group )
  
  ####
  m_unbalanced <- lm(LCBD~subset_group,data=LCBD_df)
  summary(m_unbalanced)
  
  m_unbalanced_ad <- lm(avg_dis~subset_group,data=LCBD_df)
  summary(m_unbalanced_ad)
  ####
  n <- rowSums(as.matrix(1-gower(subset_group)))+1
  
  rep <- sqrt(LCBD_env/min(LCBD_env))*min(n)

  new_mat <- NULL
  for (i in 1:nrow(subset_mat)) {
    new_mat <- rbind(new_mat,subset_mat[rep(i, rep[[i]]), ])
  }
  
  dist<- sqrt(vegdist(new_mat,dissim))
  #dist <- sqrt(ESS(new_mat))
  #dist <- sqrt(dis.nness(new_mat))
  ord <- wcmdscale(dist,eig=T)
  unique.ord.point <- ord$point
  unique.ord.negaxes <- ord$negaxes
  #unique.ord.point <- unique.ord.point[!grepl("\\.",rownames(unique.ord.point)),]
  #unique.ord.negaxes <- unique.ord.negaxes[!grepl("\\.",rownames(unique.ord.negaxes)),]

  cleaned.dist <- sqrt(dist(unique.ord.point)^2)
  #cleaned.dist <- sqrt(dist(unique.ord.point)^2 - dist(unique.ord.negaxes)^2)
  cleaned.dist[is.nan(cleaned.dist)] <- 0
  
  lower.tri.dist <- as.matrix(cleaned.dist)[lower.tri(cleaned.dist)]
  
  SSTotal <- sum(lower.tri.dist^2)/nrow(unique.ord.point)
  
  centroid_dist <- betadisper(cleaned.dist,group=rep("a",nrow(unique.ord.point)),type="centroid")
  centroid_dist <- centroid_dist$distances[!grepl("\\.",rownames(unique.ord.point))]
  #w <-apply(as.matrix(n),2,function(x) sqrt(x/(x-1)))
  w <- 1
  LCBD <-(centroid_dist*w)^2/sum((centroid_dist*w)^2)*sum((centroid_dist*w)^2)
  
  m_env_controlled <- lm(LCBD~subset_group)
  summary(m_env_controlled)

  #new_LCBD <- centroid_dist$distances^2/SSTotal #trial
  #sum(centroid_dist$distances^2/SSTotal)
  
  avg_dis <- list()
  for (i in 1:length(subset_group)){
    w <- sqrt(LCBD.comp(gower(subset_group[-i]))$LCBD/min(LCBD.comp(gower(subset_group[-i]))$LCBD))
    avg_dis[[i]] <- weighted.mean(as.matrix(vegdist(subset_mat,dissim))[-i,i],w=w)
  }

  avg_dis <- do.call(rbind,avg_dis)
  m_env_controlled_ad <- lm(avg_dis~subset_group)
  summary(m_env_controlled_ad)
  
  #reg_df <- site_dist(subset_mat,subset_group,"bray","gower",0.5)
  #m_env_controlled_ad <- lm(avg_dis~subset_group,data=reg_df)
  
  temp_result <- c(summary(m_balanced)$coefficients[1,1],
                   summary(m_balanced)$coefficients[2,c(1,4)],
                   summary(m_unbalanced)$coefficients[1,1],
                   summary(m_unbalanced)$coefficients[2,c(1,4)],
                   summary(m_env_controlled)$coefficients[1,1],
                   summary(m_env_controlled)$coefficients[2,c(1,4)],
                   summary(m_balanced_ad)$coefficients[1,1],
                   summary(m_balanced_ad)$coefficients[2,c(1,4)],
                   summary(m_unbalanced_ad)$coefficients[1,1],
                   summary(m_unbalanced_ad)$coefficients[2,c(1,4)],
                   summary(m_env_controlled_ad)$coefficients[1,1],
                   summary(m_env_controlled_ad)$coefficients[2,c(1,4)])
  
  temp_result <- as.data.frame(t(temp_result))
  temp_result$balanced_eff <- temp_result[,2]/temp_result[,1]*100
  temp_result$unbalanced_eff <- temp_result[,5]/temp_result[,4]*100
  temp_result$env_controlled_eff <- temp_result[,8]/temp_result[,7]*100
  temp_result$sample <- sample
  temp_result$resample_repeat <- m
  temp_result$matrix <- j
  temp_result_df <- rbind(temp_result_df,temp_result)
  }

#boxplot(temp_result_df$balanced_eff,temp_result_df$unbalanced_eff,temp_result_df$env_controlled_eff)
#boxplot(temp_result_df[,2],temp_result_df[,5],temp_result_df[,8])

result_df <- rbind(result_df,temp_result_df)

#t.test(temp_result_df$balanced_eff,temp_result_df$env_controlled_eff,paired = T)
#t.test(temp_result_df$result_df[,2],temp_result_df$result_df[,8],paired = T)
#t.test(temp_result_df$result_df[,5],temp_result_df$result_df[,8],paired = T)
  }
  }

library(ggplot2)
plot_df <- data.frame(Effect = c(result_df[,2],result_df[,5],result_df[,8]), Group = c(rep("Equal effort",nrow(result_df)),rep("Unequal effort - Unweighted",nrow(result_df)),rep("Unequal effort - Weighted",nrow(result_df))),sample=rep(result_df$sample,3))
plot_df$id <- interaction(plot_df$Group,plot_df$sample)

p_dist <- ggplot(plot_df)+
  geom_boxplot(aes(y=Effect,x=sample,Group=id,fill=Group))+
  ylab("Effect based on distance to centroid")+
  xlab("Number of assemblages of Habitat A")+
  theme_classic()

plot(p_dist)

library(ggplot2)
plot_df_ad <- data.frame(Effect = c(result_df[,11],result_df[,14],result_df[,17]), Group = c(rep("Equal effort",nrow(result_df)),rep("Unequal effort - Unweighted",nrow(result_df)),rep("Unequal effort - Weighted",nrow(result_df))),sample=rep(result_df$sample,3))
plot_df_ad$id <- interaction(plot_df_ad$Group,plot_df_ad$sample)

p <- ggplot(plot_df_ad)+
  geom_boxplot(aes(y=Effect,x=sample,group=id,fill=Group))+
  ylab("Effect based on average pairwise dissimilarity")+
  xlab("Number of assemblages of Habitat A")+
  theme_classic()

plot(p)

library(ggplot2)
plot_df_ad <- data.frame(Effect = c(result_df[,8],result_df[,11],result_df[,14]), Group = c(rep("Equal effort",nrow(result_df)),rep("Unequal effort - Unweighted",nrow(result_df)),rep("Unequal effort - Weighted",nrow(result_df))),sample=rep(result_df$sample,3))
plot_df_ad$id <- interaction(plot_df_ad$group,plot_df_ad$sample)

p_ad <- ggplot(plot_df_ad)+
  geom_boxplot(aes(y=Effect,x=sample,group=id,fill=group))+
  ylab("Effect based on average pairwise dissimilarity")+
  xlab("Number of assemblages of Habitat A")+
  theme_classic()

plot(p_ad)

library(ggpubr)
p <- ggarrange(p_dist,p,common.legend = T, legend="bottom")
plot(p)
ggsave("p3.tiff",dpi=800,,compression="lzw")
### nestedness

env <- seq(20,30,by=2)
site_initial <- 10
sp_richness <- 10

for (i in 1:100) {
  site_mat_list <- list()
  for (j in 1:length(env)){
    sp_pool <- 2*env[[j]]-20
    sp_pool <- round(sp_pool)
    site_mat <- matrix(0,nrow=site_initial,ncol=sp_pool)
    for (k in 1:nrow(site_mat)) {
    site_mat[k,sample(1:sp_pool,sp_richness)] <- 1 
    }
    colnames(site_mat) <- paste0("sp",1:sp_pool)
    site_mat_list[[j]] <- site_mat
  }
  
  mat <- do.call(rbind.fill,lapply(site_mat_list,as.data.frame))
  mat[is.na(mat)]<-0
  env_data <- rep(env,each=site_initial)
  
  avg_dis <- as.matrix(vegdist(mat,dissim))
  avg_dis[col(avg_dis)==row(avg_dis)] <- NA
  avg_dis <- colMeans(avg_dis,na.rm=T)
  
  m_balanced_ad <- lm(avg_dis~env_data)
  summary(m_balanced_ad)
  
  reg_df <- site_dist(mat,env_data,"bray","gower",0)
  m_env_controlled_ad <- lm(avg_dis~env,data=reg_df)
  summary(m_env_controlled_ad)
}

### mid-domain

env <- seq(20,30,5)
site_initial <- 20
sp_richness <- 10
dissim = "bray"

result_df <- list()
for (i in 1:100) {
  message("i = ",i)
  site_mat_list <- list()
  
  mat <- matrix(0,site_initial*length(env),40)
  
  for (j in 1:nrow(mat)) {
    if (j <= nrow(mat)/3) {
    mat[j,sample(1:20,sp_richness)] <- 1
    }  else if (j > nrow(mat)/3*2){
      mat[j,sample(21:40,sp_richness)] <- 1
    } else {
      mat[j,sample(1:40,sp_richness)] <- 1
    }
  }
  
  env_data <- rep(env,each=site_initial)
  
  avg_dis <- as.matrix(vegdist(mat,dissim))
  avg_dis[col(avg_dis)==row(avg_dis)] <- NA
  avg_dis <- colMeans(avg_dis,na.rm=T)
  
  env_dist <- as.matrix(gower(env_data))
  env_dist[col(env_dist)==row(env_dist)] <- NA
  env_dist <- colMeans(env_dist,na.rm=T)
  
  plot(env_data,avg_dis)
  m_balanced_ad <- lm(avg_dis~env_data+I(env_data^2))
  summary(m_balanced_ad)

  reg_df <- site_dist(mat,env_data,"bray","gower",0)
  m_env_controlled_ad <- lm(avg_dis~env+I(env^2),data=reg_df)
  summary(m_env_controlled_ad)
  plot(reg_df$env,reg_df$avg_dis)
  
  result_df[[i]] <- c(summary(m_balanced_ad)$coefficients[1,1],
                       c(summary(m_balanced_ad)$coefficients[2:3,c(1,4)]),
                       summary(m_env_controlled_ad)$coefficients[1,1],
                       c(summary(m_env_controlled_ad)$coefficients[2:3,c(1,4)]))
}

result_df <- do.call(rbind,result_df)
