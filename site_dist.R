site_dist <- function(mat,env,dissim_mat,dissim_env,ref,quadratic.term=T,trans = T) {
  require(BAT)
  require(vegan)
  avg_dis <- list()
  dis_matrix <- vegdist(mat,dissim)
  
  if (trans == T) {
    dis_matrix <- apply(as.matrix(dis_matrix),2,function(x) (x*(length(x)-1)+0.5)/length(x))
  }
  
  if (length(ref) == 1){
    ref <- rep(ref,nrow(mat))
  }
  
  for (k in 1:nrow(dis_matrix)){
    pair_dist <- as.matrix(dis_matrix)[-k,k]
    if (dissim_env == "gower") {
    env_dist <- as.matrix(gower(as.matrix(env)))[-k,k]
    } else {
      env_dist <- as.matrix(vegdist(env,dissim_env))[-k,k]
    }
    reg_df <- as.data.frame(cbind(pair_dist,env_dist,env[-k]))

    if (quadratic.term == T) {
    dis_controlled_m <- glmmTMB(pair_dist~env_dist+I(env_dist^2),family=beta_family(),reg_df)
    } else {
      dis_controlled_m <-  glmmTMB(pair_dist~env_dist,family=beta_family(),reg_df)
    }
    avg_dis[[k]] <- predict(dis_controlled_m,newdata = data.frame(env_dist=ref[[k]]),type="response")
  }
  result_df <- do.call(rbind,avg_dis)
  result_df <- data.frame(avg_dis=as.numeric(result_df),env)
  return(result_df)
}