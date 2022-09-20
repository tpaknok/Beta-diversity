site_dist <- function(mat,env,dissim_mat,dissim_env,ref,quadratic.term=T,trans = T,link="log",weight=NULL) {
  require(BAT)
  require(vegan)
  require(mgcv)
  avg_dis <- list()
  dis_matrix <- vegdist(mat,dissim)
  env <- as.data.frame(env)
  #if (trans == T) {
  #dis_matrix <- apply(as.matrix(dis_matrix),2,function(x) (x*(length(x)-1)+0.5)/length(x))
  #}
  
  #if (length(ref) == 1){
    #ref <- rep(ref,nrow(mat))
  #}
  
  if (is.null(weight)) {
    w <- rep(1,nrow(mat))
  } else {
    w <- weight
  }
  
  env_list <- list()
  for (env_var in colnames(env)) {
    if (is.numeric(env[,env_var])) {
      numeric_env <- env[,env_var]
      max_num_env <- max(numeric_env)
      min_num_env <- min(numeric_env)
      range_env <- rbind(min_num_env,max_num_env)
  
      env_list[[env_var]]<- apply(range_env,2,function (x) seq(x[[1]],range_env[[2]],length.out=11))
    } else {
        env_factor <- env[,env_var]
        env_list[[env_var]] <- unique(env_factor)
    }
    }
  
  env_space <- expand.grid(env_list)
  env_space[,sapply(env_space,is.ordered)] <- c(apply(as.data.frame(env_space[,sapply(env_space,is.ordered)]),2,as.numeric))
  
  #st.val <- c(0.005, -1e-05)
  for (k in 1:nrow(dis_matrix)){
    total <- t(apply(mat,1,function(x) colSums(rbind(x,mat[k,]))))
    total <- specnumber(total)[-k]
    pair_dist <- as.matrix(dis_matrix)[-k,k]
    #pair_dist <- 1-pair_dist
    #if (dissim_env == "gower") {
      #env_dist <- as.matrix(gower(as.matrix(env)))[-k,k]
    #} else {
      #env_dist <- as.matrix(vegdist(env,dissim_env))[-k,k]
    #}
    
    #env_diff <- SED <- NULL
    #for (env_var in colnames(env)) {
      #temp_SED <- min(apply(vegdist(env[,env_var],"manhattan"),2,max))
      #if (is.factor(env[,env_var])) {      
        #temp_env_diff <-as.matrix(gower(env[,env_var]))
        #diag(temp_env_diff) <- NA
        #temp_env_diff <- na.omit(temp_env_diff[,k])
        
      #} else {
        #temp_env_diff <- env[,env_var]-env[k,env_var]
        #temp_env_diff[[k]] <- NA
        #temp_env_diff <- na.omit(temp_env_diff)
      #}
      #env_diff <- cbind(env_diff,temp_env_diff)
    #}
    #env_diff <- as.data.frame(env_diff)
    #colnames(env_diff) <- paste0(colnames(env),"_diff")

    pair_df <- data.frame(pair_dist,env=as.data.frame(env)[-k,])
    
    pair_df[,sapply(pair_df,is.ordered)] <- c(apply(as.data.frame(pair_df[,sapply(pair_df,is.ordered)]),2,as.numeric))
    #predict_between <- unique(pair_df[abs(pair_df$env_diff)==max(abs(env_diff)),env_var])

    #reg_df$scaled_dist <- c(scale(reg_df$pair_dist,center=min(reg_df[,1]),scale = max(reg_df[,1])-min(reg_df[,1])))
    #reg_df$scaled_env_dist <- c(scale(reg_df$env_dist,center=min(reg_df[,2]),scale = max(reg_df[,2])-min(reg_df[,2])))
    #scaled_ref <- (ref-min(reg_df[,2]))/(max(reg_df[,2])-min(reg_df[,2]))
    
   # if (quadratic.term == T) {
      #dis_controlled_m <- glmmTMB(pair_dist~env_dist+I(env_dist^2),family=betabinomial,reg_df,weights=w[-k])
      #dis_controlled_m <- glm(pair_dist~env_dist+I(env_dist^2),family=binomial(link),reg_df,weights = total)
      #dis_controlled_m <-  glm(pair_dist~env_dist+I(env_dist^2),data=reg_df,family=binomial(link="log"),weights=w[-k])
    #} else {
      #dis_controlled_m <- glm(pair_dist~env_dist,family = gaussian(link = "log"),data=reg_df,weights=w[-k],start=st.val)
      #dis_controlled_m <-  randomForest(pair_dist~env_dist,data=reg_df,weights=w[-k])
      #dis_controlled_m <-neuralnet(scaled_dist~scaled_env_dist,data=reg_df)
      #dis_controlled_m <- nls(pair_dist ~ NLS.expoDecay(env_dist, a, k),data = pair_df)
      #dis_controlled_m <- drm(pair_dist ~ env_dist,fct = DRC.expoDecay(),data=pair_df)
      dis_controlled_m <- gam(pair_dist ~ s(env,k=3),data=pair_df)
      summary(dis_controlled_m)
    #}
    #avg_dis[[k]] <- unique(predict(dis_controlled_m,newdata = data.frame(env_dist=rep(ref[[k]],length(env_pair_w))),type="response"))
    #avg_dis[[k]] <- unique(predict(dis_controlled_m,newdata = data.frame(env_dist=rep(ref[[k]],nrow(mat)-1)),type="response"))
    #avg_dis[[k]] <- unique(compute(dis_controlled_m,covariate = data.frame(scaled_env_dist=rep(scaled_ref[[k]],nrow(mat)-1)))$net.result)
    #avg_dis[[k]] <- avg_dis[[k]]*(max(reg_df[,1])-min(reg_df[,1]))+min(reg_df[,1])
    predict_within <- data.frame(ifelse(is.ordered(env),as.numeric(env[k,]),env[k,]))
    colnames(predict_within) <- colnames(env)
    within_env <- predict(dis_controlled_m,newdata=predict_within)
    within_env <- ifelse(within_env > 1, 1, ifelse(within_env < 0, 0, within_env))
    
    predict_net<- as.data.frame(env_space)
    colnames(predict_net) <- colnames(env)
    net_env <- c(predict(dis_controlled_m,newdata=predict_net))
    net_env <- ifelse(net_env > 1, 1, ifelse(net_env < 0, 0, net_env))
    total_env <- sum(net_env)/length(net_env)
    
    between_env <- sum(net_env-c(within_env))/(length(net_env)-1) #doesn't matter including within-env in sum - always zero. -1 to remove within env
    avg_dis[[k]] <- data.frame(c(within_env,between_env,total_env),c("Within","Between","Total"),env[k,],k)
    names(avg_dis[[k]]) <- c("avg_dis","type",colnames(env),"Site")
    
    #avg_dis[[k]] <- ifelse(avg_dis[[k]] > 1, 1,ifelse(avg_dis[[k]] < 0,0,avg_dis[[k]]))
  }
  result_df <- do.call(rbind,avg_dis)
  #result_df <- data.frame(avg_dis=as.numeric(result_df),env)
  return(result_df)
}
