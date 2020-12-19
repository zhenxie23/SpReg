DisVecAVSim <- function(DatR, VarR, DatY, VarY, ncore, cutoffR, cutoffYR, ...){
  
  library(sp);
  library(gstat);
  library(geosphere);
  library(mvtnorm);
  
  # extract data from dot 
  dots <- list(...);
  beta0 = dots$beta;
  psill0 = dots$Grouped_theta1_mle;
  range0 = dots$Grouped_theta2_mle;
  nuget0 = dots$Grouped_nuget_mle;
  Grouped_vgm_model = dots$Grouped_vgm_model;
  V_QMLE = dots$V_QMLE;
  
  num_of_group <- length(DatR);
  num_of_beta <- length(DatY)/length(DatR);
  
  for(g in 1:1){
  
    # computing the covariance matrix of R and Y
    dist_R <- spDists(DatR[[g]]);
    num_R <- length(DatR[[g]][[VarR]]);
    m_R <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR[[g]][[VarR]])));
    K_R <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){
      
      value <- Semivariogram_function(Grouped_vgm_model[g], dist_R[i,j], 
                                      psill0[g], range0[g], nuget0[g], cov = TRUE);
      return(value);
      
    })});
    
    K_Y <- list();
    m_Y <- list();
    
    for(v in 1:num_of_beta){
      
      g_Y <- g+(v-1)*num_of_group;
      num_Y <- length(DatY[[g_Y]][[VarY]]);
      dist_Y <- spDists(DatY[[g_Y]]);
      m_Y[[v]] <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatY[[g_Y]][[VarY]] )));
      
      dist_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){distGeo(DatY[[g_Y]]@coords[i, ],DatR[[g]]@coords[j, ])/1000})});
      m1 <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatR[[g]][[VarR]] )));
      m2 <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR[[g]][[VarR]] )));
      K1 <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){
        value <- Semivariogram_function(Grouped_vgm_model[g], dist_YR[i,j], 
                                        psill0[g], range0[g], nuget0[g], cov = TRUE);
        return(value);
      })});
      K1 <- matrix(nrow = num_Y, ncol = num_R, K1);
      K2 <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){
        value <- Semivariogram_function(Grouped_vgm_model[g], dist_R[i,j], 
                                        psill0[g], range0[g], nuget0[g], cov = TRUE);
        return(value);
      })});
      
      # Got the predicted rainfall
      DatY[[g_Y]][[VarR]]  <- c(m1+K1%*%solve(K2, DatR[[g]][[VarR]]-m2));
      res <- DatY[[g_Y]][[VarY]]-beta0[v]*DatY[[g_Y]][[VarR]];
      u_res <- var(res);
      
      K_Y[[v]] <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){
        
        value <- (beta0[v]^2)*Semivariogram_function(Grouped_vgm_model[g], dist_Y[i,j], 
                                                   psill0[g], range0[g], nuget0[g], cov = TRUE)+
                  u_res;
        return(value);
      })});
    
    }
    
    J = 100;
    gamma_s <- NULL;

    # R_s, Y_s are independent draws sampling from the Gaussian Random Field
    for(j in 1:J){
      
      DatR[[g]]$R_s <- t(rmvnorm(1, mean = c(t(m_R)), sigma = K_R));
      g_R_s <- gstat(NULL, "R_s", R_s~1, DatR[[g]]);
      one_gamma_s <- NULL;
      
      for(v in 1:num_of_beta){
        
        g_Y <- g+(v-1)*num_of_group;
        # sample_cov_YR is the sample cross covariogram of Y and R
        DatY[[g_Y]]$Y_s <- t(rmvnorm(1, mean = c(t(m_Y[[v]])), sigma = K_Y[[v]]));
        g_YR_s <- gstat(g_R_s, "Y_s", Y_s~1, DatY[[g_Y]]);
        sample_cov_YR_s <- variogram(g_YR_s, covariogram = TRUE, cutoff = cutoffYR);
        
        #sample_cov_YR <- variogram(g_YR, covariogram = TRUE);
        Id_YR <- "R_s.Y_s";
        sample_cov_YR_s <- sample_cov_YR_s[sample_cov_YR_s$id == Id_YR,];
        sample_cov_YR_s <- sample_cov_YR_s[order(sample_cov_YR_s$dist),];
        
        one_gamma_s <- c(one_gamma_s, sample_cov_YR_s$gamma);
        
        rm(sample_cov_YR_s);
        
      }
      
      # sample_cov_R is the sample covariogram of R
      sample_cov_R_s <- variogram(g_R_s, covariogram = TRUE, cutoff = cutoffR);
      index_distance_0 <- which(sample_cov_R_s$dist == 0);
      index_new <- setdiff(1:nrow(sample_cov_R_s), index_distance_0);
      index_new <- c(index_new, index_distance_0[1]);
      sample_cov_R_s <- sample_cov_R_s[index_new, ];
      sample_cov_R_s <- sample_cov_R_s[order(sample_cov_R_s$dist), ];
      one_gamma_s <- c(one_gamma_s, sample_cov_R_s$gamma);
    
      gamma_s <- cbind(gamma_s, one_gamma_s);
      rm(sample_cov_R_s);
      rm(one_gamma_s);
      
    }

  # Asymptotic matrix
  num_of_gamma <- nrow(gamma_s);
  one_B <- matrix(nrow = num_of_gamma, ncol = num_of_gamma, 0);
  for(i in 1:num_of_gamma){
    for(j in 1:num_of_gamma){
      one_B[i,j] <- cov(gamma_s[i,],gamma_s[j,]);
    }
  }
  
  }
  
  return(one_B);
  
}