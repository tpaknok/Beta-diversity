site_dist <- function(mat,env,dissim_mat,dissim_env,ref,quadratic.term=T) {
  require(BAT)
  require(vegan)
  avg_dis <- list()
  dis_matrix <- vegdist(mat,dissim)
  for (k in 1:nrow(dis_matrix)){
    pair_dist <- as.matrix(dis_matrix)[-k,k]
    if (dissim_env == "gower") {
    env_dist <- as.matrix(gower(env))[-k,k]
    } else {
      env_dist <- as.matrix(vegdist(env,dissim_env))[-k,k]
    }
    reg_df <- as.data.frame(cbind(pair_dist,env_dist))
    if (quadratic.term == T) {
    dis_controlled_m <- lm(pair_dist~env_dist+I(env_dist^2),reg_df)
    } else {
      dis_controlled_m <- lm(pair_dist~env_dist,reg_df)
    }
    avg_dis[[k]] <- predict(dis_controlled_m,newdata = data.frame(env_dist=ref))
  }
  result_df <- do.call(rbind,avg_dis)
  result_df <- data.frame(avg_dis=as.numeric(result_df),env)
  return(result_df)
}