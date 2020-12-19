ABCMHEstimation <- function(DatR, VarR, DatY, VarY, epsilon, W_mat, cutoffR, cutoffYR, ...){
  
  library(Formula);
  library(geosphere);
  library(gstat);
  library(sp);
  library(rgdal);
  library(rARPACK);
  library(mixlm);
  library(MASS);
  library(georob);
  library(mvtnorm);
  library(optimx);
  library(spatstat);
  library(mvQuad);
  #library(MinDistFunctions);
  library(doParallel);
  library(distributions3);
  
  # extract data from dot 
  dots <- list(...);
  TSBS_point_estimation = dots$TSBS_point_estimation;
  TSBS_sd = dots$TSBS_sd;
  Grouped_theta1_mle = dots$Grouped_theta1_mle;
  Grouped_theta2_mle = dots$Grouped_theta2_mle;
  Grouped_nuget_mle = dots$Grouped_nuget_mle;
  Grouped_vgm_model = dots$Grouped_vgm_model;
  V_QMLE = dots$V_QMLE;
  
  ##### 1 Prepare Dist, Gamma, B, starting values 
  
  #### 2.1.1 Dist, Gamma -> from gstat nonparametric variogram value
  
  Formula_R <- as.formula(paste(VarR,"~","1"));
  Formula_Y <- as.formula(paste(VarY,"~","1"));
  
  num_of_group <- length(DatR);
  num_of_beta <- length(DatY)/length(DatR);
  
  Dist <- list();
  Gamma <- list();
  Num_of_bins <- list();
  
  for(g in 1:num_of_group){
    
    g_R <- gstat(NULL, VarR, Formula_R, DatR[[g]]);
    one_dist <- NULL;
    one_gamma <- NULL;
    one_num_of_bins <- NULL;
    
    for(v in 1:num_of_beta){
      
      # sample_cov_YR is the sample cross covariogram of Y and R
      g_YR <- gstat(g_R, VarY, Formula_Y, DatY[[g+(v-1)*num_of_group]]);
      sample_cov_YR <- variogram(g_YR, covariogram = TRUE, cutoff = cutoffYR);
      #sample_cov_YR <- variogram(g_YR, covariogram = TRUE);
      Id_YR <- paste(VarR,".",VarY,sep = "");
      sample_cov_YR <- sample_cov_YR[sample_cov_YR$id == Id_YR,];
      sample_cov_YR <- sample_cov_YR[order(sample_cov_YR$dist),];
      SAMPLE_COV_YR <<- sample_cov_YR;
      
      one_dist <- c(one_dist, sample_cov_YR$dist);
      one_gamma <- c(one_gamma, sample_cov_YR$gamma);
      one_num_of_bins <- c(one_num_of_bins, nrow(sample_cov_YR));
      
      rm(sample_cov_YR);
      
    }
    
    # sample_cov_R is the sample covariogram of R
    sample_cov_R <- variogram(g_R, covariogram = TRUE, cutoff = cutoffR);
    index_distance_0 <- which(sample_cov_R$dist == 0);
    index_new <- setdiff(1:nrow(sample_cov_R), index_distance_0);
    index_new <- c(index_new, index_distance_0[1]);
    sample_cov_R <- sample_cov_R[index_new, ];
    SAMPLE_COV_R <<- sample_cov_R;
    #sample_cov_R <- variogram(g_R, covariogram = TRUE);
    sample_cov_R <- sample_cov_R[order(sample_cov_R$dist), ];
    
    one_dist <- c(one_dist, sample_cov_R$dist);
    one_gamma <- c(one_gamma, sample_cov_R$gamma);
    one_num_of_bins <- c(one_num_of_bins, nrow(sample_cov_R));
    
    # Dist is the default lags of the sample covariogram/cross variogram 
    Dist[[g]] <- one_dist;
    # Gamma_YR is the value of sample covariogram/cross variogram
    Gamma[[g]] <- one_gamma;
    # Num_of_bins() stores the number of bins for each group
    Num_of_bins[[g]] <- one_num_of_bins;
    
    rm(one_dist, one_gamma, one_num_of_bins);
    
  }
  
  
  AV <- W_mat;
  
  J = 100;
  
  t1 <- Sys.time();
  
  #epsilon = 0.25;
  value_list <- NULL;
  
  beta_group <- NULL;
  psill_group <- NULL;
  range_group <- NULL;
  Gamma_list <- list();
  prob_list <- NULL;
  prob_new_list <- NULL;
  
  for(j0 in 1:J){
    
    print(j0);
    
    # sample beta and theta
    if(num_of_beta > 1){
      new_beta <- mvrnorm(1, mu = c(TSBS_point_estimation), Sigma = diag(TSBS_sd));
    }else{
      new_beta <- mvrnorm(1, mu = TSBS_point_estimation, Sigma = (TSBS_sd^2));
    }
    
    #print(new_beta);
    
    new_psill <- NULL;
    new_range <- NULL;
    for(g in 1:num_of_group){
      
      psill_g <- -1.0;## positive number
      range_g <- -1.0;
      
      while(psill_g < 0 || range_g < 0){
        
        theta_g <- mvrnorm(1, mu = c(Grouped_theta1_mle[g],Grouped_theta2_mle[g]), Sigma = V_QMLE[[g]]);
        psill_g <- theta_g[1]; range_g <- theta_g[2];
      }
      
      new_psill <- c(new_psill, theta_g[1]);
      new_range <- c(new_range, theta_g[2]);
      
    }
    
    #print(new_psill);
    #print(new_range);
    new_Gamma <- list();

    for(g in 1:num_of_group){
      
      one_gamma <- NULL;
      index <- 0;
      
      for(v in 1:num_of_beta){
        for(i in 1: Num_of_bins[[g]][v]){
          
          index <- index + 1;
          gamma_i <- new_beta[v]*Semivariogram_function(Grouped_vgm_model[g], Dist[[g]][index], 
                                                        new_psill, new_range, 0, cov = TRUE);
          one_gamma <- c(one_gamma, gamma_i);
          
        }
      }
      
      for(i in 1: Num_of_bins[[g]][num_of_beta+1]){
        index <- index + 1;
        gamma_i <- Semivariogram_function(Grouped_vgm_model[g], Dist[[g]][index], 
                                          new_psill, new_range, 0, cov = TRUE);
        one_gamma <- c(one_gamma, gamma_i);
        
      }
      new_Gamma[[g]] <- one_gamma;
    }
    
    
    value <- 0;
    for(g in 1:num_of_group){
      
      value <- value + c(t(Gamma[[g]]-new_Gamma[[g]])%*%AV%*%(Gamma[[g]]-new_Gamma[[g]]));
      
    }
    
    Gamma_list[[j0]] <- new_Gamma[[1]];
    
    value_list <- c(value_list, value);
    psill_group <- c(psill_group, new_psill[g]);
    range_group <- c(range_group, new_range[g]);
    #Gamma_list[[j0]] <- new_Gamma[[g]];
    
    if(j0 == 1){
      beta_group <- matrix(nrow = num_of_beta, ncol = 1, new_beta);
      if(num_of_beta > 1){
        prob_par <- dmvnorm(new_beta, mean = c(TSBS_point_estimation), sigma = diag(TSBS_sd^2));
      }else{
        prob_par <- dnorm(new_beta, mean = TSBS_point_estimation, sd = TSBS_sd);
        prob_par_new <- prob_par;
      }
      for(g in 1:num_of_group){
        
        prob_theta <- dmvnorm(c(new_psill[g],new_range[g]), mean = c(Grouped_theta1_mle[g],Grouped_theta2_mle[g]), sigma = V_QMLE[[g]])
        prob_par <- prob_par * prob_theta;
        
      }
    }else{
      if(value < epsilon){
        
        # sample u ~ U[0,1]
        u <- runif(1, min = 0, max = 1);
        # compute the proposal distribution q(theta*|theta)
        if(num_of_beta > 1){
          prob_par_new <- dmvnorm(new_beta, mean = c(TSBS_point_estimation), sigma = diag(TSBS_sd^2));
        }else{
          prob_par_new <- dnorm(new_beta, mean = TSBS_point_estimation, sd = TSBS_sd);
        }
        for(g in 1:num_of_group){
          
          prob_theta_new <- dmvnorm(c(new_psill[g],new_range[g]), mean = c(Grouped_theta1_mle[g],Grouped_theta2_mle[g]), sigma = V_QMLE[[g]])
          prob_par_new <- prob_par_new * prob_theta_new;
          
        }
        # compute the acceptance ratio
        ar <- prob_par/prob_par_new;
        #print(ar);
        
        if(u < ar){
          beta_group <- cbind(beta_group, new_beta);
          prob_par <- prob_par_new;
        }else{
          beta_j <- beta_group[,j0-1];
          beta_group <- cbind(beta_group, beta_j);
        }
      }else{
        beta_j <- beta_group[,j0-1];
        beta_group <- cbind(beta_group, beta_j);
      }
    }
    
    #prob_list <- c(prob_list, prob_par);
    #prob_new_list <- c(prob_new_list, prob_par_new);
    
  }
  
  abc_mean <- apply(beta_group,1,mean);
  
  t2 <- Sys.time();
  print(t2-t1);
  
  Output <- list(coefficient = abc_mean,
                 beta_group = beta_group,
                 value_list = value_list);
  
  return(Output);
  
}
