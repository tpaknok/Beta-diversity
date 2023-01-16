site_dist <- function(formula,random=NULL,correlation=NULL,spatial_var = NULL,
                      mat=NULL,env_data,dissim="jaccard",C=0.5,dissim_mat = NULL,
                      env_space = NULL,
                      family="gaussian",
                      weights=NULL,
                      length.cont=25,
                      ...) {
  require(BAT)
  require(vegan)
  require(mgcv)
  require(betaC)
  
  avg_dis <- list()
  env <- as.data.frame(env_data)
  
  #####DEVELOPMENTAL PURPOSES
  if (dissim == "beta_C" & is.null(dissim_mat)) {
    message("rarefying....")
    beta_pairwise<- beta_stand(mat, func = list( "beta_C"), setsize=2,args = list(C=C),summarise=F)
    dis_matrix <- matrix(NA,nrow(mat),nrow(mat))
    dis_matrix[lower.tri(avg_dis_matrix,diag=F)] <- beta_pairwise
    dis_matrix[upper.tri(avg_dis_matrix,diag=F)] <- t(avg_dis_matrix)[upper.tri(avg_dis_matrix,diag=F)]
  } 
  ######  
  
  if (dissim != "beta_C" & is.null(dissim_mat)) {
    dis_matrix <- vegdist(mat,dissim)
  } 
  
  if (!is.null(dissim_mat)){
    dis_matrix <- dissim_mat
  }
  
  if (is.null(weights)) {
    w = NULL
  } else {
    w = weights
  }
  
  library(stringr)
  
  seq_env <-apply(as.data.frame(colnames(env)),1,function(x) grepl(x,as.character(formula)[[3]]))
  
  env <- as.data.frame(env[,seq_env])
  colnames(env) <- colnames(env_data)[seq_env]
  
  env_list <- list()
  for (env_var in colnames(env)) {
    if (is.numeric(env[,env_var])) {
      numeric_env <- env[,env_var]
      max_num_env <- max(numeric_env)
      min_num_env <- min(numeric_env)
      range_env <- rbind(min_num_env,max_num_env)
  
      env_list[[env_var]]<- apply(range_env,2,function (x) seq(x[[1]],range_env[[2]],length.out=length.cont))
    } else {
        env_factor <- env[,env_var]
        env_list[[env_var]] <- unique(env_factor)
        
        if (is.ordered(env_list[[env_var]])) {
          env_list[[env_var]] <- as.numeric(env_list[[env_var]])
        }
      }
    }
  
  env_space <- expand.grid(env_list)
  message("finish creating environmental space")
  
  if (length(which(sapply(env_data,is.ordered) == T)) > 1) {
    message("Detected more than one ordinal factor")
    env_data[,sapply(env_data,is.ordered)] <- sapply(as.data.frame(env_data[,sapply(env_data,is.ordered)]),as.numeric)
  }
  
  
  if (length(which(sapply(env_data,is.ordered) == T)) == 1) {
    message("Detected one ordinal factor")
    env_data[,sapply(env_data,is.ordered)] <- as.vector(sapply(as.data.frame(env_data[,sapply(env_data,is.ordered)]),as.numeric))
  }
  
  for (k in 1:nrow(dis_matrix)){
    message(k,"/",nrow(dis_matrix))
    pair_dist <- as.matrix(dis_matrix)[-k,k]

    pair_df <- data.frame(pair_dist,env_data[-k,])
    colnames(pair_df)[-1] <- colnames(env_data)
    
    env_name <- colnames(pair_df)[-1]
    num_seq <- which(sapply(pair_df,is.numeric))[-1] - 1
    
    if (is.null(random) & is.null(correlation) & is.null(w)) {
      dis_controlled_m <- gam(formula=formula,data=pair_df,family=family,...)
      } 
    
    if (is.null(random) & is.null(correlation) & !is.null(w)) {
      pair_df$w_site <- w[-k]
      dis_controlled_m <- gam(formula=formula,data=pair_df,family=family,weights=w_site,...)
      } 
    
    if (!(is.null(random) & is.null(correlation)) & is.null(w)) {
      dis_controlled_m <- gamm(formula=formula,data=pair_df,family=family,...)
    } 
    
    if (!(is.null(random) & is.null(correlation)) & !is.null(w)) {
      pair_df$w_site <- w[-k]
      dis_controlled_m <- gamm(formula=formula,data=pair_df,family=family,weights=w_site,...)
      } 
    
    predict_net<- env_space
    colnames(predict_net) <- colnames(env)
    
    if (is.null(random) & is.null(correlation)) {
      net_env <- c(predict(dis_controlled_m,newdata=predict_net,type="response"))
    } else {
      net_env <- c(predict(dis_controlled_m$gam,newdata=predict_net,type="response"))
    }
    
    total_env <- sum(net_env)/length(net_env)

        avg_dis[[k]] <- data.frame(total_env,c(env[k,]),k)
    colnames(avg_dis[[k]]) <- c("avg_dis",colnames(env),"Site")
    if (ncol(env) == 1) {
      colnames(avg_dis)[[2]] <- colnames(env_var)
    }
  }
  result_df <- do.call(rbind,avg_dis)
  result_df <- cbind(result_df,env_data[,spatial_var])
}
