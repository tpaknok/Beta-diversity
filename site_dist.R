site_dist <- function(formula,random=NULL,correlation=NULL,spatial_var = NULL,
                      mat=NULL,env_data,dissim="jaccard",C=0.5,dissim_mat = NULL,
                      env_space = NULL,
                      family="gaussian",
                      weights=NULL,
                      length.cont=25,
                      ...) {
  
convert_env_df <- function(env,length.cont=length.cont) {
    if (is.numeric(env)) { #for numeric variables
      numeric_env <- env
      max_num_env <- max(numeric_env) #get max
      min_num_env <- min(numeric_env) # get min
      range_env <- rbind(min_num_env,max_num_env)
      
      env_list<- apply(range_env,2,function (x) seq(x[[1]],range_env[[2]],length.out=length.cont)) #generate an evenly-spcaed gradient 
    } else { #for non-numeric variables....
      env_factor <- env 
      env_list <- unique(env_factor) #get unique factors 
      
      if (is.ordered(env_list)) { #for ranked factor (or ordinal variables)
        env_list <- as.numeric(env_list) #change them to numbers
      }
    }
    return(env_list)
  }
  
  require(BAT)
  require(vegan)
  require(mgcv)
  require(betaC)
  
  avg_dis <- list()
  env <- as.data.frame(env_data)
  
  #####DEVELOPMENTAL PURPOSES. This is for rarefaction beta diversity. Need more test.
  if (dissim == "beta_C" & is.null(dissim_mat)) {
    message("rarefying....")
    beta_pairwise<- beta_stand(mat, func = list( "beta_C"), setsize=2,args = list(C=C),summarise=F)
    dis_matrix <- matrix(NA,nrow(mat),nrow(mat))
    dis_matrix[lower.tri(avg_dis_matrix,diag=F)] <- beta_pairwise
    dis_matrix[upper.tri(avg_dis_matrix,diag=F)] <- t(avg_dis_matrix)[upper.tri(avg_dis_matrix,diag=F)]
  } 
  ######  
  
  if (dissim != "beta_C" & is.null(dissim_mat)) {
    dis_matrix <- vegdist(mat,dissim) #calculate dissimilarity
  } 
  
  if (!is.null(dissim_mat)){
    dis_matrix <- dissim_mat 
  }
  
  if (is.null(weights)) {
    w = NULL 
  } else {
    w = weights #do sites have different weight?
  }
  
  library(stringr)
  
  seq_env <-apply(as.data.frame(colnames(env)),1,function(x) grepl(x,as.character(formula)[[3]]))
  
  env <- as.data.frame(env[,seq_env])
  colnames(env) <- colnames(env_data)[seq_env]
  
  env_list <- lapply(env,function(x) convert_env_df(x,length.cont=length.cont)) #create the niche space for projection
  
  env_space <- expand.grid(env_list)
  message("finish creating environmental space")
  
  if (length(which(sapply(env_data,is.ordered) == T)) > 1) {
    message("Detected more than one ordinal factor")
    env_data[,sapply(env_data,is.ordered)] <- sapply(as.data.frame(env_data[,sapply(env_data,is.ordered)]),as.numeric) #turn ranked factor to numeric variable 
  }
  
  
  if (length(which(sapply(env_data,is.ordered) == T)) == 1) {
    message("Detected one ordinal factor")
    env_data[,sapply(env_data,is.ordered)] <- as.vector(sapply(as.data.frame(env_data[,sapply(env_data,is.ordered)]),as.numeric)) #turn ranked factor to numeric variable 
  }
  
  for (k in 1:nrow(dis_matrix)){ #generate Uniche of each site
    message(k,"/",nrow(dis_matrix)) 
    pair_dist <- as.matrix(dis_matrix)[-k,k]

    pair_df <- data.frame(pair_dist,env_data[-k,])
    colnames(pair_df)[-1] <- colnames(env_data)
    
    env_name <- colnames(pair_df)[-1]
    num_seq <- which(sapply(pair_df,is.numeric))[-1] - 1
    
    if (is.null(random) & is.null(correlation) & is.null(w)) { #gam if no random effect & correlation structure
      dis_controlled_m <- gam(formula=formula,data=pair_df,family=family,...)
      } 
    
    if (is.null(random) & is.null(correlation) & !is.null(w)) { #gam with weight
      pair_df$w_site <- w[-k]
      dis_controlled_m <- gam(formula=formula,data=pair_df,family=family,weights=w_site,...)
      } 
    
    if (!(is.null(random) & is.null(correlation)) & is.null(w)) { #gamm without weight
      dis_controlled_m <- gamm(formula=formula,data=pair_df,family=family,...)
    } 
    
    if (!(is.null(random) & is.null(correlation)) & !is.null(w)) { #gamm with weight
      pair_df$w_site <- w[-k]
      dis_controlled_m <- gamm(formula=formula,data=pair_df,family=family,weights=w_site,...)
      } 
    
    predict_net<- env_space 
    colnames(predict_net) <- colnames(env)
    
    if (is.null(random) & is.null(correlation)) {
      net_env <- c(predict(dis_controlled_m,newdata=predict_net,type="response")) #project the model using the environmental space
    } else {
      net_env <- c(predict(dis_controlled_m$gam,newdata=predict_net,type="response"))
    }
    
    net_env[net_env > 1] <- 1 #cap projected dissimilarity between 0-1
    net_env[net_env < 0] <- 0
    
    total_env <- sum(net_env)/length(net_env) #average projected dissimilarity across the environmental space, thus obtain Uniche

    avg_dis[[k]] <- data.frame(total_env,c(env[k,]),k)
    colnames(avg_dis[[k]]) <- c("avg_dis",colnames(env),"Site")
  }
  result_df <- do.call(rbind,avg_dis)
  result_df <- cbind(result_df,env_data[,spatial_var])
}
