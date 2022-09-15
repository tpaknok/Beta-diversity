site_dist <- function(mat,env,dissim_mat,dissim_env,ref,quadratic.term=T,trans = T,link="log",weight=NULL) {
  require(BAT)
  require(vegan)
  require(ks)
  require(drc)
  require(aomisc)
  avg_dis <- list()
  dis_matrix <- vegdist(mat,dissim)
  env <- as.numeric(env)
  #if (trans == T) {
    #dis_matrix <- apply(as.matrix(dis_matrix),2,function(x) (x*(length(x)-1)+0.5)/length(x))
  #}
  
  if (length(ref) == 1){
    ref <- rep(ref,nrow(mat))
  }
  
  if (is.null(weight)) {
    w <- rep(1,nrow(mat))
  } else {
    w <- weight
  }
  
  #st.val <- c(0.005, -1e-05)
  for (k in 1:nrow(dis_matrix)){
    total <- t(apply(mat,1,function(x) colSums(rbind(x,mat[k,]))))
    total <- specnumber(total)[-k]
    pair_dist <- as.matrix(dis_matrix)[-k,k]
    #pair_dist <- 1-pair_dist
    if (dissim_env == "gower") {
    env_dist <- as.matrix(gower(as.matrix(env)))[-k,k]
    } else {
      env_dist <- as.matrix(vegdist(env,dissim_env))[-k,k]
    }
    
    pair_df <- as.data.frame(cbind(pair_dist,env_dist,env=env[-k],total=total,id=seq(1:(nrow(mat)-1))))
    #reg_df$scaled_dist <- c(scale(reg_df$pair_dist,center=min(reg_df[,1]),scale = max(reg_df[,1])-min(reg_df[,1])))
    #reg_df$scaled_env_dist <- c(scale(reg_df$env_dist,center=min(reg_df[,2]),scale = max(reg_df[,2])-min(reg_df[,2])))
    #scaled_ref <- (ref-min(reg_df[,2]))/(max(reg_df[,2])-min(reg_df[,2]))
    
    if (quadratic.term == T) {
    #dis_controlled_m <- glmmTMB(pair_dist~env_dist+I(env_dist^2),family=betabinomial,reg_df,weights=w[-k])
      #dis_controlled_m <- glm(pair_dist~env_dist+I(env_dist^2),family=binomial(link),reg_df,weights = total)
    #dis_controlled_m <-  glm(pair_dist~env_dist+I(env_dist^2),data=reg_df,family=binomial(link="log"),weights=w[-k])
    } else {
      #dis_controlled_m <- glm(pair_dist~env_dist,family = gaussian(link = "log"),data=reg_df,weights=w[-k],start=st.val)
      #dis_controlled_m <-  randomForest(pair_dist~env_dist,data=reg_df,weights=w[-k])
      #dis_controlled_m <-neuralnet(scaled_dist~scaled_env_dist,data=reg_df)
      #dis_controlled_m <- nls(pair_dist ~ NLS.expoDecay(env_dist, a, k),data = pair_df)
      dis_controlled_m <- drm(pair_dist ~ env_dist,fct = DRC.expoDecay(),data=pair_df)
      
    }
    #avg_dis[[k]] <- unique(predict(dis_controlled_m,newdata = data.frame(env_dist=rep(ref[[k]],length(env_pair_w))),type="response"))
    avg_dis[[k]] <- unique(predict(dis_controlled_m,newdata = data.frame(env_dist=rep(ref[[k]],nrow(mat)-1)),type="response"))
    #avg_dis[[k]] <- unique(compute(dis_controlled_m,covariate = data.frame(scaled_env_dist=rep(scaled_ref[[k]],nrow(mat)-1)))$net.result)
    #avg_dis[[k]] <- avg_dis[[k]]*(max(reg_df[,1])-min(reg_df[,1]))+min(reg_df[,1])
    avg_dis[[k]] <- ifelse(avg_dis[[k]] > 1, 1,ifelse(avg_dis[[k]] < 0,0,avg_dis[[k]]))
  }
  result_df <- do.call(rbind,avg_dis)
  result_df <- data.frame(avg_dis=as.numeric(result_df),env)
  return(result_df)
}