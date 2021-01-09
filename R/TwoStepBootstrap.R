#' Estimation and Inference for Two-Step Estimators
#'
#' @description This function implements the esmation and inference for two-step estimators including Krig-and-regress(OLS), Krig-and-regress(GLS), two-step bootstrap.
#'
#' @usage TwoStepBootstrap(DatR, VarR, DatY, VarY, variogram.model,
#'                    is.cov.misspecified, is.den.misspecified,
#'                    plot.start.value = TRUE, cutoff.R, cutoff.res,
#'                    start.value.method = 2, projected = FALSE)
#'
#' @param DatR explanatory variable R, a spatial object, see \link[sp:coordinates]{coordintes()}
#' @param VarR name of variable R
#' @param DatY outcome variable Y, a spatial object, see \link[sp:coordinates]{coordintes()}
#' @param VarY name of variable Y
#' @param variogram.model variogram model type, e.g. "Exp", "Sph", "Gau", "Mat"
#' @param is.cov.misspecified logical; if TRUE, the covariance function is misspecified
#' @param is.den.misspecified logical; if TRUE, the density function is not Gaussian distribution.
#' @param plot.start.value logical, if TRUE, plot the variogram and the fitted variogram curve corresponding to starting values
#' @param cutoff.R cutoff for sample variogram of variable R
#' @param cutoff.res cutoff for sample variogram of regression residuals
#' @param start.value.method fitting method, see \link[gstat:fit.variogram]{fit.variogram()}
#' @param projected logical; if FALSE, data are assumed to be unprojected, meaning decimal longitude/latitude. For projected data, Euclidian distances are computed, for unprojected great circle distances(km) are computed.
#' @return \item{\code{vario.par.point.est}}{point estimates for variogram parameters(psill, range)}
#'         \item{\code{vario.par.var.mat}}{estimated variance-covariance matrix for variogram parameters(psill, range)}
#'         \item{\code{ols.point.est}}{point estimates for Krig-and-OLS estimator}
#'         \item{\code{ols.var.mat}}{estimated variance-covariance matrix for Krig-and-OLS estimator}
#'         \item{\code{gls.point.est}}{point estimates for Krig-and-GLS estimator}
#'         \item{\code{gls.var.mat}}{estimated variance-covariance matrix for Krig-and-GLS estimator}
#'         \item{\code{tsbs.point.est}}{point estimates for Two-Step Bootstrap estimator}
#'         \item{\code{tsbs.var.mat}}{estimated variance-covariance matrix for Two-Step Bootstrap estimator}
#' @examples
#' library(SpReg)
#' rivers <- read.csv(system.file("extdata", "rivers.csv", package = "SpReg"))
#'
#' rivers <- rivers[which(rivers$FOR_NLCD<100),]
#' train_set <- sample(1:558,277);
#' test_set <- setdiff(1:558,train_set);
#' rivers_train <- rivers[train_set, ];
#' rivers_test <- rivers[test_set, ];
#'
#' DatR <- rivers_train[,c("FOR_NLCD", "LAT_DD", "LON_DD")];
#' DatR$X <- log((DatR$FOR_NLCD)/(100-DatR$FOR_NLCD));
#' sp::coordinates(DatR) <- ~LON_DD+LAT_DD;
#' sp::proj4string(DatR) =  "+proj=longlat +datum=WGS84";
#' DatY <- rivers_test[, c("CL","LAT_DD","LON_DD")];
#' DatY$Y <- log(DatY$CL);
#' sp::coordinates(DatY) <- ~LON_DD+LAT_DD;
#' sp::proj4string(DatY) =  "+proj=longlat +datum=WGS84";
#'
#' TwoStep_Results <- TwoStepBootstrap(DatR, "X", DatY, "Y", "Exp", FALSE, FALSE,
#'                                     cutoff.R = 295, cutoff.res = 40);
#'
#' @seealso \link{sp}, \link{gstat}
#' @export

TwoStepBootstrap <- function(DatR, VarR, DatY, VarY, variogram.model,
                             is.cov.misspecified, is.den.misspecified,
                             plot.start.value= TRUE, cutoff.R, cutoff.res,
                             start.value.method = 2, projected = FALSE){

  if(projected == FALSE){
    sp::proj4string(DatR) =  "+proj=longlat +datum=WGS84";
    sp::proj4string(DatY) = "+proj=longlat +datum=WGS84";
  }


  ## estimation of theta obtained via classic geostatistical methods, as the starting point of optimization
  vgm1_formula <- as.formula(paste(VarR,"~","1"));
  vgm1 <- gstat::variogram(vgm1_formula, DatR, cutoff = cutoff.R);
  #vgm1 <- variogram(vgm1_formula, DatR);
  model1 <- gstat::vgm(psill = 4, variogram.model, range = 10);
  mfit <- gstat::fit.variogram(vgm1, model1, fit.method = start.value.method);

  # plot the staring points
  if(plot.start.value == TRUE){
    plot(vgm1, model = mfit);
  }

  starting_values <- c(psill = mfit$psill[1], range = mfit$range[1]);
  #starting_values <- c(psill = 10, range = 10);

  # distance matrix
  ## need to divide into several groups
  Dist <- list(sp::spDists(DatR));

  # Rainfall
  ## need to divide into several groups
  Rainfall <- list(matrix(nrow = length(DatR[[VarR]]), ncol = 1, c(DatR[[VarR]])));

  Vsample <- list();
  Formula_R <- as.formula(paste(VarR,"~","1"));
  g_R = gstat::gstat(NULL, VarR, Formula_R, DatR);
  Cov_sample_RR <- gstat::variogram(g_R, covariogram = TRUE, cutoff = cutoff.R);
  #Cov_sample_RR <- variogram(g_R, covariogram = TRUE);
  Cov_sample_RR <- Cov_sample_RR[order(Cov_sample_RR$dist),];
  Vsample[[1]] <- SampleCovarianceMatrix(Cov_sample_RR, Dist[[1]]);

  #VSAMPLE <<- Vsample[[1]]

  if(isTRUE(variogram.model == "Exp")){

    function_lists <- list(

      function(par, gs){

        psill <- par[1];range <- par[2];
        #Rainfall <- gs$Rainfall;
        Dist <- gs$Dist;
        Vsample <- gs$Vsample;
        #M <- nrow(Rainfall);
        M <- nrow(Vsample);
        #m <- matrix(nrow = M, ncol = 1, rep(0,M, mean(Rainfall)));
        #m <- matrix(nrow = M, ncol = 1, 0);
        Omega <- matrix(nrow = M, ncol = M, 0);
        Omega <- sapply(1:M, function(i){sapply(1:M, function(j){psill*exp(-Dist[i,j]/range)})});

        Omega_eig <- eigen(Omega)$values;
        Omega_det <- sum(log(abs(Omega_eig)));
        value <- 0.5*Omega_det + 0.5*(matrixcalc::matrix.trace(solve(Omega, Vsample)));
        #value <- Omega_det + t(Rainfall-m)%*%solve(Omega,Rainfall-m);

        print(psill);
        print(range);
        print(value);
        print("------------------------");

        return(value);
      }

    )

  }

  gs_dat <- list(Vsample = Vsample[[1]], Dist = Dist[[1]]);
  #gs_dat <- list(Rainfall = Rainfall[[1]], Dist = Dist[[1]]);

  lower_values <- c(psill = 0.0001, range = 0.0001);
  upper_values <- c(psill = Inf, range = Inf);

  # minimization of Kullback-Leibler Distance
  mle_theta_result <- optimx::optimx(starting_values, fn = function_lists[[1]],
                                  #lower = lower_values, upper = upper_values,
                                  #method = c("L-BFGS-B"),
                                  method = c("BFGS"),
                                  control = list(maximize = FALSE, fnscale = 1.0),
                                  gs = gs_dat);

  # Point Estimation
  psill_mle <-mle_theta_result$psill;
  range_mle <- mle_theta_result$range;

  if(plot.start.value == TRUE){

    plot(vgm1,model = gstat::vgm(psill_mle, "Exp", range_mle));

  }

  num_R <- length(DatR@data[[VarR]]);
  num_Y <- length(DatY@data[[VarY]]);
  dist_train_test <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){geosphere::distGeo(DatY@coords[i, ],DatR@coords[j, ])/1000})});
  #dist_train_test <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ sqrt((DatR@coords[j,1]-DatY@coords[i,1])^2+(DatR@coords[j,2]-DatY@coords[i,2])^2)  })});
  m1 <- matrix(nrow = num_Y, ncol = 1, rep(mean(Rainfall[[1]])));
  m2 <- matrix(nrow = num_R, ncol = 1, rep(mean(Rainfall[[1]])));
  #m1 <- matrix(nrow = num_Y, ncol = 1, 0);
  #m2 <- matrix(nrow = num_R, ncol = 1, 0);
  K1 <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){psill_mle*exp(-dist_train_test[i,j]/range_mle)})});
  K2 <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_mle*exp(-Dist[[1]][i,j]/range_mle)})});
  R_hat <- m1+K1%*%solve(K2)%*%(Rainfall[[1]]-m2);

  beta_mle <- lm(DatY[[VarY]] ~ R_hat)$coefficients;


  # OLS analysis
  ols_model <- lm(DatY[[VarY]] ~ R_hat);

  # GLS analysis
  DatY$residuals <- ols_model$residuals;
  vgm_res <- gstat::variogram(residuals~1, DatY, cutoff = cutoff.res);
  #vgm_res <- variogram(residuals~1, DatY);
  mfit_res <- gstat::fit.variogram(vgm_res, gstat::vgm(variogram.model),
                                   fit.method = start.value.method);
  psill_res <- mfit_res$psill[2]; range_res = mfit_res$range[2];
  nug_res <- mfit_res$psill[1];
  dist_test <- sp::spDists(DatY);

  #vgm_res <- variogram(residuals~1, DatY, cutoff = cutoff.res, covariogram = TRUE);
  #Cov_sample_RR <- variogram(g_R, covariogram = TRUE, cutoff = cutoff.res);
  #vgm_res <- vgm_res[order(vgm_res$dist),];
  #K3 <- SampleCovarianceMatrix(vgm_res, dist_test);

  K3 <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){ nug_res+psill_res*exp(-dist_test[i,j]/range_res)})});
  my_weights <- solve(K3);
  #gls_model <- lm.gls(DatY$Y ~ R_hat, W = my_weights);
  X_hat <- cbind(rep(1,length(R_hat)), R_hat);
  gls_estimate <- solve(t(X_hat)%*%my_weights%*%X_hat)%*%(t(X_hat)%*%my_weights%*%DatY[[VarY]]);
  gls_sd <- solve(t(X_hat)%*%my_weights%*%X_hat);

  ### Computation of asymtotic variance of Quasi ML estimator

  W <- list();
  #m <- matrix(nrow = num_R, ncol = 1, rep(mean(Rainfall[[1]])));
  W[[1]] <- Rainfall[[1]] - m2;
  V <- list();
  V[[1]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_mle*exp(-Dist[[1]][i,j]/range_mle)})});
  FDV <- list();
  FDV[[1]] <- list();
  FDV[[1]][[1]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){exp(-Dist[[1]][i,j]/range_mle)})});
  FDV[[1]][[2]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_mle*exp(-Dist[[1]][i,j]/range_mle)*(Dist[[1]][i,j]/range_mle^2)})});

  SDV <- list();
  SDV[[1]] <- list();
  SDV[[1]][[1]] <- list();
  SDV[[1]][[2]] <- list();
  SDV[[1]][[1]][[1]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){ 0 })});
  SDV[[1]][[1]][[2]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){exp(-Dist[[1]][i,j]/range_mle)*(Dist[[1]][i,j]/(range_mle^2))})});
  SDV[[1]][[2]][[1]] <- SDV[[1]][[1]][[2]]
  SDV[[1]][[2]][[2]] <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_mle*exp(-Dist[[1]][i,j]/range_mle)*(((Dist[[1]][i,j]^2)/(range_mle^4))-(2*Dist[[1]][i,j]/(range_mle^3)))})});

  V_QMLE <- QMLE_Variance(W, V, FDV, SDV, Vsample, is.cov.misspecified, is.den.misspecified);

  J = 100;

  beta_bootstrap <- NULL;

  for(b_j in 1:J){

    theta_j <- MASS::mvrnorm(1, mu = c(psill_mle,range_mle), Sigma = V_QMLE[[1]]);

    psill_j <- theta_j[1];
    range_j <- theta_j[2];

    if(psill_j <= 0 || range_j <= 0){
      next;
    }

    K1 <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){psill_j*exp(-dist_train_test[i,j]/range_j)})});
    K2 <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_j*exp(-Dist[[1]][i,j]/range_j)})});

    if(class(try(solve(K2)))[1] == "try-error"){
      next;
    }

    R_j_hat <- m1+K1%*%solve(K2)%*%(Rainfall[[1]]-m2);

    bootstrap_index <- sample(1:num_Y, num_Y, replace = TRUE, prob = NULL);
    Y_B <- DatY@data[[VarY]][bootstrap_index];
    R_B <- R_j_hat[bootstrap_index];

    beta_bootstrap <- cbind(beta_bootstrap, lm(Y_B ~ R_B)$coefficients);

  }

  bootstrap_mean <- matrix(nrow = 1, ncol = nrow(beta_bootstrap), 0);
  bootstrap_cov <- matrix(nrow = nrow(beta_bootstrap), ncol = nrow(beta_bootstrap), 0);

  for(i in 1:nrow(beta_bootstrap)){

    bootstrap_mean[i] <- mean(beta_bootstrap[i, ]);

  }

  for(i in 1:nrow(beta_bootstrap)){
    for(j in 1:nrow(beta_bootstrap)){
      bootstrap_cov[i,j] <- cov(beta_bootstrap[i, ],beta_bootstrap[j, ]);
    }
  }

  # formatting the outputs

  var_string <- names(ols_model$coefficients);
  bootstrap_mean <- c(bootstrap_mean);
  gls_estimate <- c(gls_estimate);
  names(bootstrap_mean) <- var_string;
  names(gls_estimate) <- var_string;

  rownames(gls_sd) <- var_string;
  colnames(gls_sd) <- var_string;
  rownames(bootstrap_cov) <- var_string;
  colnames(bootstrap_cov) <- var_string;

  vario_par_var_mat <- V_QMLE[[1]];
  rownames(vario_par_var_mat) <- c("psill", "range");
  colnames(vario_par_var_mat) <- c("psill", "range");

  Output = list(
    # variogram_R = starting_point_variogram,
    vario.par.point.est = c(psill = psill_mle, range = range_mle),
    vario.par.var.mat = vario_par_var_mat,
    ols.point.est = ols_model$coefficients,
    ols.var.mat = vcov(ols_model),
    gls.point.est = gls_estimate,
    gls.var.mat = gls_sd,
    tsbs.point.est = bootstrap_mean,
    tsbs.var.mat = bootstrap_cov
  );

  return(Output);

}
