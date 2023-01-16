library(BiodiversityR)
library(vegan)
source("site_dist.R")
data(BCI)
data(BCI.env)

Start <- Sys.time()
predict_pb <- site_dist(pair_dist ~ s(elevation,k=3) + s(convex,k=3) + s(slope,k=3),
                        spatial_var = c("UTM.EW","UTM.NS"),
                        correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),
                        mat=BCI,env_data=BCI.env,
                        length.cont=10,control=lmeControl(opt = 'optim',maxIter = 1e8, msMaxIter = 1e8))

End <- Sys.time()
End-Start

m_niche <- gls(avg_dis~elevation+convex+slope,correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),data=predict_pb)
summary(m_niche)
#vif(m_niche)

library(ggeffects)
ele_niche <- ggemmeans(m_niche,"elevation")
convex_niche <- ggemmeans(m_niche,"convex")
slope_niche <- ggemmeans(m_niche,"slope")

niche_pre <- rbind(ele_niche,convex_niche,slope_niche)
niche_pre$env <- c(rep("elevation",nrow(ele_niche)),
                 rep("convex",nrow(convex_niche)),
                 rep("slope",nrow(slope_niche)))
niche_pre$Significance <- c(rep("Insignificant",nrow(ele_niche)),
                 rep("Insignificant",nrow(convex_niche)),
                 rep("Significant",nrow(slope_niche)))
niche_pre$Model <- "Niche"
####

v_dist <- as.matrix(vegdist(BCI,"jaccard"))
diag(v_dist) <- NA
v_distm <- apply(v_dist,2,function(x) mean(x,na.rm=T))

analysis_BCI <- cbind(BCI.env,avg_dis=v_distm)

m_obs <- gls(v_distm~elevation+convex+slope,correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),data=analysis_BCI)
summary(m_obs)

library(adespatial)
LCBD <- LCBD.comp(vegdist(BCI,"jaccard"))
analysis_BCI$LCBD <- LCBD$LCBD
m_LCBD <- gls(LCBD~elevation+convex+slope,correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),data=analysis_BCI)
summary(m_LCBD)

CD <- betadisper(as.dist(vegdist(BCI,"jaccard")),group = c(rep("A",nrow(BCI))),type="centroid",sqrt.dist=T)$distance
analysis_BCI$CD <- CD
m_CD <- gls(CD~elevation+convex+slope,correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),data=analysis_BCI)
summary(m_CD)

result_table <- rbind(summary(m_niche)$tTable,summary(m_obs)$tTable,summary(m_LCBD)$tTable,summary(m_CD)$tTable)
result_table <- round(result_table,3)
write.csv(result_table,"SI_table.csv")

####
library(ggeffects)
ele_obs <- ggemmeans(m_obs,"elevation")
convex_obs <- ggemmeans(m_obs,"convex")
slope_obs <- ggemmeans(m_obs,"slope")

obs_pre <- rbind(ele_obs,convex_obs,slope_obs)
obs_pre$env <- c(rep("elevation",nrow(ele_obs)),
                 rep("convex",nrow(convex_obs)),
                 rep("slope",nrow(slope_obs)))
obs_pre$Significance <- c(rep("Insignificant",nrow(ele_obs)),
                 rep("Insignificant",nrow(convex_obs)),
                 rep("Insignificant",nrow(slope_obs)))
obs_pre$Model <- "Observed"

overall_pre <- rbind(niche_pre,obs_pre)

########################################
library(tidyverse)
predict_pb_long <- predict_pb[,-5:-7] %>% pivot_longer(!avg_dis)
obs_long <- predict_pb[,-5:-7] %>% pivot_longer(!avg_dis)
predict_pb_long$Model <- "Niche"

analysis_BCI_long <- analysis_BCI[,c(3:5,8)] %>% pivot_longer(!avg_dis)
analysis_BCI_long$Model <- "Observed"

Overall_long_df <- rbind(predict_pb_long,analysis_BCI_long)
colnames(Overall_long_df)[[2]] <- "env"
###
env.labs <- c("Elevation (m)", "Convexity (m)", "Slope (°)")
names(env.labs) <- c("elevation", "convex", "slope")

Overall_long_df <- Overall_long_df[Overall_long_df$env != "UTM.EW" & Overall_long_df$env != "UTM.NS",]
analysis_BCI_long$label <- "(a)"
analysis_BCI_long$label[analysis_BCI_long$name == "elevation"] <- "(b)"
analysis_BCI_long$label[analysis_BCI_long$name == "slope"] <- "(c)"

label_df <- data.frame(name=unique(sort(analysis_BCI_long$name)),label=sort(unique(analysis_BCI_long$label)))

library(ggplot2)
density_p <- ggplot(analysis_BCI_long,aes(x=value))+
  geom_density()+
  geom_text(data=label_df,aes(label=label),x=-Inf,y=Inf,hjust=-0.5,vjust=1)+
  facet_wrap(~name,scale="free",nrow=5,labeller = labeller(name=env.labs))+
  ylab("Density")+
  xlab("Environmental conditions")+
  theme_classic()
plot(density_p)

label_df <- data.frame(env=unique(sort(analysis_BCI_long$name)),label=c("(d)","(e)","(f)"))

library(ggplot2)
empirical_p <- ggplot(overall_pre,aes(x=x,y=predicted))+
  geom_point(data=Overall_long_df,aes(x=value,y=avg_dis,colour=Model),alpha=0.5)+
  geom_line(aes(linetype=Significance,colour=Model),size=1.2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high,fill=Model),alpha=0.25)+
  geom_text(data=label_df,aes(label=label),x=-Inf,y=Inf,hjust=-0.5,vjust=1)+
  facet_wrap(~env,scale="free_x",nrow=3,labeller = labeller(env=env.labs))+
  ylab("Uniqueness")+
  xlab("Environmental conditions")+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_colour_manual(values=c("red","blue"),labels=c(expression(U[niche]),expression(U[observed])))+
  scale_fill_manual(values=c("red","blue"),labels=c(expression(U[niche]),expression(U[observed])))+
  scale_x_continuous(n.breaks=4)+
  theme_classic()+
  theme(legend.position="bottom")+
  guides(linetype = "none")
plot(empirical_p)

library(ggpubr)
p <- ggarrange(density_p,empirical_p,ncol=2,common.legend=T,legend="bottom")
plot(p)

ggsave("empirical_revised.tiff",width=7,height=7,dpi=800,compression="lzw")

