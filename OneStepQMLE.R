OneStepQMLE <- function(DatR, VarR, DatY, VarY, Variogram_Model, 
                        TSBS_Starting_Value, is_cov_misspecified, is_den_misspecified, 
                        cutoffR, cutoffYR, cutoffY,
                        Plot = TRUE, Start_Value_Method = 2, projected = FALSE, ...){
  # DatR: a spatial object, see function coordinates()
  # VarR: name of variable R
  # DatY: a spatial object, see function coordinates()
  # VarY: name of variable Y
  # Variogram_Model: type of variogram model, model type, e.g. "Exp", "Sph", "Gau", "Mat". Calling vgm() 
  #            without a model argument returns a data.frame with available models, see function vgm()
  # is_cov_misspecified: logical; if TRUE, the covariance function is misspecified
  # is_den_misspecified: logical; if TRUE, the density function is not Gaussian distribution.
  # Plot_Start_Value: logical, plot the variogram and the fitted variogram curve corresponding to starting values
  # Start_Value_Method: fitting method, used by gstat, see function fit.variogram()
  # projected: logical; if FALSE, data are assumed to be unprojected, meaning decimal longitude/latitude. 
  #            For projected data, Euclidian distances are computed, for unprojected great circle distances(km)
  
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
  library(matrixcalc);
  library(MisalignFuncR)
  
  if(projected == FALSE){
    proj4string(DatR) =  "+proj=longlat +datum=WGS84";
    proj4string(DatY) = "+proj=longlat +datum=WGS84";
  }
  
  # read the starting values
  psill0 <- TSBS_Starting_Value[1];
  range0 <- TSBS_Starting_Value[2];
  beta0 <- TSBS_Starting_Value[3];
  beta0_std <- TSBS_Starting_Value[4];
  starting_values <- c(psill = unname(psill0), 
                       range = unname(range0), 
                       nuget = 0, 
                       beta = unname(beta0));
  
  print(starting_values);
  
  # extract data from dot 
  dots <- list(...);
  DatY$residuals <- dots$Residual;
  
  # distance matrix
  ## need to divide into several groups
  Dist_R <- spDists(DatR);
  Dist_Y <- spDists(DatY);
  num_R <- length(DatR@data[[VarR]]);
  num_Y <- length(DatY@data[[VarY]]);
  Dist_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){distGeo(DatY@coords[i, ],DatR@coords[j, ])/1000})});
  
  
  # Rainfall
  ## need to divide into several groups
  Rainfall <- list(matrix(nrow = length(DatR@data[[VarR]]), ncol = 1, c(DatR@data[[VarR]])));
  
  Formula_R <- as.formula(paste(VarR,"~","1"));
  Formula_Y <- as.formula(paste(VarY,"~","1"));
  
  # sample_cov_R is the sample covariogram of R
  g_R <- gstat(NULL, VarR, Formula_R, DatR);
  g_Y <- gstat(NULL, VarY, Formula_Y, DatY);
  g_YR <- gstat(g_R, VarY, Formula_Y, DatY);
  
  Vsample <- list();
  Vsample[[1]] <- matrix(nrow = (num_Y+num_R), ncol = (num_Y+num_R), 0);
  
  Cov_sample_RR <- variogram(g_R, cutoff = cutoffR, covariogram = TRUE);
  Cov_sample_RR <- Cov_sample_RR[order(Cov_sample_RR$dist),];
  Vsample[[1]][((num_Y+1):(num_Y+num_R)),((num_Y+1):(num_Y+num_R))] <- SampleCovarianceMatrix(Cov_sample_RR, Dist_R);
  
  COV_RR <<- Cov_sample_RR;
  
  Cov_sample_YY <- variogram(g_Y, cutoff = cutoffY, covariogram = TRUE);
  Cov_sample_YY <- Cov_sample_YY[order(Cov_sample_YY$dist),];
  Vsample[[1]][1:num_Y,1:num_Y] <- SampleCovarianceMatrix(Cov_sample_YY, Dist_Y);
  
  COV_YY <<- Cov_sample_YY;
  
  # Cov_sample_YR is the sample cross covariogram of Y and R
  Cov_sample_YR <- variogram(g_YR, cutoff = cutoffYR, covariogram = TRUE);
  Id_YR <- paste(VarR,".",VarY,sep = "");
  Cov_sample_YR <- Cov_sample_YR[Cov_sample_YR$id == Id_YR,];
  Cov_sample_YR <- Cov_sample_YR[order(Cov_sample_YR$dist),];
  Vsample[[1]][1:num_Y,(num_Y+1):(num_R+num_Y)] <- SampleCovarianceMatrix(Cov_sample_YR, Dist_YR);
  Vsample[[1]][(num_Y+1):(num_R+num_Y),1:num_Y] <- t(Vsample[[1]][1:num_Y,(num_Y+1):(num_R+num_Y)]);
  
  COV_YR <<- Cov_sample_YR;
  
  V_sample <<- Vsample[[1]];
  
  ############################
  # Modeling the covariance function of error term Sigma_res
  m1 <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatR@data[[VarR]] )));
  m2 <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR@data[[VarR]] )));
  K1 <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){psill0*exp(-Dist_YR[i,j]/range0)})});
  K2 <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill0*exp(-Dist_R[i,j]/range0)})});
  
  # Got the predicted rainfall
  #DatY$R  <- c(m1+K1%*%solve(K2, DatR[[VarR]]-m2));
  
  #DatY$residuals <- DatY$Y-beta0*DatY$R;
  
  Sigma_res <- diag(rep(sd(DatY$residuals), num_Y));
  
  #SIGMA_RES <<- Sigma_res;
  
  #vgm_res <- variogram(residuals~1, DatY);
  #mfit_res <- fit.variogram(vgm_res, vgm(Variogram_Model), fit.method = Start_Value_Method);
  #psill_res <- mfit_res$psill[2]; range_res = mfit_res$range[2];
  #nugget_res <- mfit_res$psill[1];
  #Sigma_res <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){psill_res*exp(-Dist_Y[i,j]/range_res)})});
  #diag(Sigma_res) <- rep((psill_res+nugget_res),num_Y);
  #SIGMA_RES <<- Sigma_res;
  #print(nugget_res);
  #P <<- plot(vgm_res,vgm(psill_res, "Exp", range_res, nugget_res));
  
  # Object Function
  if(isTRUE(Variogram_Model == "Exp")){
    
    function_lists <- list(
      
      function(par, gs){
        
        psill_var <- par[1];
        range_var <- par[2];
        nuget_var <- par[3];
        beta_var <- par[4];
        
        #Vsample <- gs$Vsample; 
        Rainfall <- gs$Rainfall;
        Rainfall <- matrix(nrow = (num_R+num_Y), ncol = 1, Rainfall);
        m <- gs$m;
        m <- matrix(nrow = (num_R+num_Y), ncol = 1, m);
        Dist_R <- gs$Dist_R;
        Dist_Y <- gs$Dist_Y; Dist_YR <- gs$Dist_YR;
        Sigma_res <- gs$Sigma_res;
        
        Omega <- matrix(nrow = (num_Y+num_R), ncol = (num_Y+num_R), 0);
        Omega[((num_Y+1):(num_Y+num_R)),((num_Y+1):(num_Y+num_R))] <-  sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 
          
          value <- Semivariogram_function(1, Dist_R[i,j], psill_var, range_var, nuget_var, cov = TRUE);
          return(value);
          
        })});
        Omega[1:num_Y,1:num_Y] <-  sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 
          
          value <- Semivariogram_function(1, Dist_Y[i,j], psill_var, range_var, nuget_var, cov = TRUE);
          value <- (beta_var^2)*value;
          return(value);
          
        })});
        Omega[1:num_Y,1:num_Y] <-  Omega[1:num_Y,1:num_Y] + Sigma_res;
        Omega[(1:num_Y),((num_Y+1):(num_R+num_Y))] <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ 
          
          value <- Semivariogram_function(1, Dist_YR[i,j], psill_var, range_var,nuget_var, cov = TRUE);
          value <- beta_var*value;
          return(value);
        
        })});
        Omega[((num_Y+1):(num_R+num_Y)),(1:num_Y)] <- t(Omega[1:num_Y,(num_Y+1):(num_R+num_Y)]);
        
        Omega_eig <- eigen(Omega)$values;
        Omega_det <- sum(log(abs(Omega_eig)));
        #value <- 0.5*Omega_det + 0.5*matrix.trace(solve(Omega, Vsample)); 
        value <- Omega_det + t(Rainfall-m)%*%solve(Omega,(Rainfall-m));
        
        print(psill_var);
        print(range_var);
        print(nuget_var);
        print(beta_var);
        print(value);
        print("-------------------------");
        
        return(value);
      }
    )
    
  }
  
  
  #gs_dat <- list(Vsample = Vsample[[1]], Dist_R = Dist_R, 
  #               Dist_Y = Dist_Y, Dist_YR = Dist_YR, Sigma_res = Sigma_res);
  gs_dat <- list(Rainfall = c(DatY[[VarY]],DatR[[VarR]]), 
                 m = c(rep(mean(DatY[[VarY]]), num_Y),
                       rep(mean(DatR[[VarR]]), num_R)),
                 Dist_R = Dist_R, Dist_Y = Dist_Y, Dist_YR = Dist_YR, Sigma_res = Sigma_res);
 
  lower_values <- c(psill = 0.0001, range = 0.0001, nuget = -Inf, beta = beta0-2*beta0_std);
  upper_values <- c(psill = 10, range = 100, nuget = Inf, beta = beta0+2*beta0_std);
  #lower_values <- c(beta0-3*beta0_std);
  #upper_values <- c(beta0+3*beta0_std);
 
  # minimization of Kullback-Leibler Distance
  mle_result <- opm(starting_values, fn = function_lists[[1]], 
                    lower = lower_values, upper = upper_values, 
                    method = c("L-BFGS-B"),
                    #method = c("BFGS"),
                    gs = gs_dat);
  
  # Point Estimation
  psill_mle <- mle_result$psill;
  range_mle <- mle_result$range;
  nuget_mle <- mle_result$nuget;
  beta_mle <- mle_result$beta;
  
  print(psill_mle);
  print(range_mle);
  print(beta_mle);
  
  ### Computation of asymtotic variance of Quasi ML estimator
  W <- list();
  m_R <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR[[VarR]])));
  m_Y <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatY[[VarY]])));
  w_R <- matrix(nrow = num_R, ncol = 1, c(DatR[[VarR]]));
  w_Y <- matrix(nrow = num_Y, ncol = 1, c(DatY[[VarY]]));
  W[[1]] <- rbind((w_Y-m_Y), (w_R-m_R));
  
  V <- list();
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 
    
    value <- Semivariogram_function(1, Dist_R[i,j], psill_mle, range_mle, nuget_mle, cov = TRUE);
    return(value);
  
  })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 
    
    value <- Semivariogram_function(1, Dist_Y[i,j], psill_mle, range_mle, nuget_mle, cov = TRUE);
    value <- (beta_mle^2)*value;
    return(value);
    
  })});
  V_Y  <- V_Y + Sigma_res;
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ 
  
    value <- Semivariogram_function(1, Dist_YR[i,j], psill_mle, range_mle, nuget_mle, cov = TRUE);
    value <- beta_mle*value;
    return(value);
    
  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  V[[1]] <- rbind(V1, V2); 
  #VV1 <<- V[[1]];
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  FDV <- list();
  FDV[[1]] <- list();
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ exp(-Dist_R[i,j]/range_mle)})});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ (beta_mle^2)*exp(-Dist_Y[i,j]/range_mle)})});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ beta_mle*exp(-Dist_YR[i,j]/range_mle)  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  FDV[[1]][[1]] <- rbind(V1, V2); 
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ psill_mle*exp(-Dist_R[i,j]/range_mle)*(Dist_R[i,j]/(range_mle^2))  })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ (beta_mle^2)*psill_mle*exp(-Dist_Y[i,j]/range_mle)*(Dist_Y[i,j]/(range_mle^2))  })});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ beta_mle*psill_mle*exp(-Dist_YR[i,j]/range_mle)*(Dist_YR[i,j]/(range_mle^2))  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  FDV[[1]][[2]] <- rbind(V1, V2); 
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 0 })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 2*beta_mle*psill_mle*exp(-Dist_Y[i,j]/range_mle)})});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ psill_mle*exp(-Dist_YR[i,j]/range_mle)  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  FDV[[1]][[3]] <- rbind(V1, V2); 
  #VV2 <<- FDV[[1]][[1]];
  
  
  SDV <- list();
  SDV[[1]] <- list();
  SDV[[1]][[1]] <- list();
  SDV[[1]][[2]] <- list();
  SDV[[1]][[3]] <- list();
  
  SDV[[1]][[1]][[1]] <- sapply(1:(num_R+num_Y), function(i){sapply(1:(num_R+num_Y), function(j){ 0 })});
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ exp(-Dist_R[i,j]/range_mle)*(Dist_R[i,j]/(range_mle^2))  })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ (beta_mle^2)*exp(-Dist_Y[i,j]/range_mle)*(Dist_Y[i,j]/(range_mle^2))  })});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ beta_mle*exp(-Dist_YR[i,j]/range_mle)*(Dist_YR[i,j]/(range_mle^2))  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  SDV[[1]][[1]][[2]] <- rbind(V1, V2);
  SDV[[1]][[2]][[1]] <- SDV[[1]][[1]][[2]];
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 0 })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 2*beta_mle*exp(-Dist_Y[i,j]/range_mle)})});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ exp(-Dist_YR[i,j]/range_mle)  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  SDV[[1]][[1]][[3]] <- rbind(V1, V2);
  SDV[[1]][[3]][[1]] <- SDV[[1]][[1]][[3]];
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ psill_mle*exp(-Dist_R[i,j]/range_mle)*(((Dist_R[i,j]^2)/(range_mle^4))-(2*Dist_R[i,j]/(range_mle^3)))} )});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ (beta_mle)^2*psill_mle*exp(-Dist_Y[i,j]/range_mle)*(((Dist_Y[i,j]^2)/(range_mle^4))-(2*Dist_Y[i,j]/(range_mle^3)))  })});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ beta_mle*psill_mle*exp(-Dist_YR[i,j]/range_mle)*(((Dist_YR[i,j]^2)/(range_mle^4))-(2*Dist_YR[i,j]/(range_mle^3)))} )});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  SDV[[1]][[2]][[2]] <- rbind(V1, V2);
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 0 })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 2*psill_mle*beta_mle*exp(-Dist_Y[i,j]/range_mle)*(Dist_Y[i,j]/(range_mle^2))})});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ psill_mle*exp(-Dist_YR[i,j]/range_mle)*(Dist_YR[i,j]/(range_mle^2))  })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  SDV[[1]][[2]][[3]] <- rbind(V1, V2);
  SDV[[1]][[3]][[2]] <- SDV[[1]][[2]][[3]];
  
  rm(V_R);rm(V_Y);rm(V_YR);rm(V1);rm(V2);
  V_R  <- sapply(1:num_R, function(j){sapply(1:num_R, function(i){ 0 })});
  V_Y  <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ 2*psill_mle*exp(-Dist_Y[i,j]/range_mle) })});
  V_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ 0 })});
  V1 <- cbind(V_Y, V_YR);
  V2 <- cbind(t(V_YR), V_R);
  SDV[[1]][[3]][[3]] <- rbind(V1, V2);
  

  V_QMLE_1 <- QMLE_Variance(W, V, FDV, SDV, Vsample, is_cov_misspecified, is_den_misspecified);
  V_QMLE_2 <- QMLE_Variance(W, V, FDV, SDV, Vsample, TRUE, FALSE);
  
  Output <- list(
    
    theta_mean = c(psill = psill_mle, range = range_mle, beta = beta_mle),
    theta_cov_1 = V_QMLE_1[[1]],
    theta_cov_2 = V_QMLE_2[[1]]);
  
  return(Output);
}