#' Estimation and Inference for Minimum-Distance Estimator
#'
#' @description This function implements the esmation and large sample inference for
#'              Minimum-Distance Estimator.
#'
#' @usage MinimumDistance(DatR, VarR, DatY, VarY, variogram.model,
#'                        MD.starting.value, cutoff.R, cutoff.YR,
#'                        start.value.method = 2, projected = FALSE)
#'
#' @param DatR explanatory variable R, a spatial object, see \link[sp:coordinates]{coordintes()}
#' @param VarR name of variable R
#' @param DatY outcome variable Y, a spatial object, see \link[sp:coordinates]{coordintes()}
#' @param VarY name of variable Y
#' @param variogram.model variogram model type, e.g. "Exp", "Sph", "Gau", "Mat"
#' @param MD.starting.value the starting point of parameters
#' @param cutoff.R cutoff for sample variogram of variable R
#' @param cutoff.YR cutoff for sample cross variogram of variable R and Y
#' @param start.value.method fitting method, see \link[gstat:fit.variogram]{fit.variogram()}
#' @param projected logical; if FALSE, data are assumed to be unprojected, meaning decimal longitude/latitude. For projected data, Euclidian distances are computed, for unprojected great circle distances(km) are computed.
#' @return \item{\code{num.obs}}{the number of observations}
#'         \item{\code{vario.par.point.est}}{point estimates for variogram parameters(psill, range)}
#'         \item{\code{vario.par.var.mat}}{estimated variance-covariance matrix for variogram parameters(psill, range)}
#'         \item{\code{md.point.est}}{point estimates for Min-Dist estimator}
#'         \item{\code{md.var.mat}}{estimated aymptotic variance for Min-Dist estimator}
#'
#' @seealso \link{sp}, \link{gstat}
#' @export

MinimumDistance <- function(DatR, VarR, DatY, VarY, variogram.model,
                            MD.starting.value, cutoff.R, cutoff.YR,
                            start.value.method = 2, projected = FALSE){

  if(projected == FALSE){
    sp::proj4string(DatR) =  "+proj=longlat +datum=WGS84";
    sp::proj4string(DatY) = "+proj=longlat +datum=WGS84";
  }

  Formula_R <- as.formula(paste(VarR,"~","1"));
  Formula_Y <- as.formula(paste(VarY,"~","1"));

  # sample_cov_R is the sample covariogram of R
  g_R = gstat::gstat(NULL, VarR, Formula_R, DatR);
  g_YR = gstat::gstat(g_R, VarY, Formula_Y, DatY);
  sample_cov_R <- gstat::variogram(g_R, covariogram = TRUE, cutoff = cutoff.R);
  #sample_cov_R <- gstat::variogram(g_R, covariogram = TRUE);
  sample_cov_R <- sample_cov_R[order(sample_cov_R$dist),];

  # sample_cov_YR is the sample cross covariogram of Y and R
  sample_cov_YR <- gstat::variogram(g_YR, covariogram = TRUE, cutoff = cutoff.YR);
  #sample_cov_YR <- gstat::variogram(g_YR, covariogram = TRUE);
  Id_YR <- paste(VarR,".",VarY,sep = "");
  sample_cov_YR <- sample_cov_YR[sample_cov_YR$id == Id_YR,];
  sample_cov_YR <- sample_cov_YR[order(sample_cov_YR$dist),];

  # Dist_YR is the default lags of the sample cross covariogram of Y and R
  Dist_YR <- list();
  Dist_YR[[1]] <- sample_cov_YR$dist;

  # Gamma_YR is the value of sample cross covariogram of Y and R
  Gamma_YR <- list();
  Gamma_YR[[1]] <- sample_cov_YR$gamma;
  num_of_bin_YR <- nrow(sample_cov_YR);

  # Dist_R is the default lags of the sample covariogram of R
  Dist_R <- list();
  Dist_R[[1]] <- sample_cov_R$dist;

  # Gamma_R is the value of sample covariogram of R
  Gamma_R <- list();
  Gamma_R[[1]] <- sample_cov_R$gamma;
  num_of_bin_R <- nrow(sample_cov_R);

  # User-defined distance function, which we aim to minimize
  if(isTRUE(variogram.model == "Exp")){

    function_lists <- list(

        function(par, gs){

            psill = par[1]; range = par[2]; beta = par[3];
            Dist_YR = gs$Dist_YR;Gamma_YR = gs$Gamma_YR;B_YR = gs$B_YR;
            Dist_R = gs$Dist_R; Gamma_R = gs$Gamma_R; B_R = gs$B_R;

            n <- length(Dist_YR);
            m <- length(Dist_R);

            g_YR <- sapply(1:n, function(i){

              Gamma_YR[i]-beta*psill*exp(-Dist_YR[i]/range)

            });

            g_R <- sapply(1:m, function(i){

              Gamma_R[i]-psill*exp(-Dist_R[i]/range)

            });

            g <- c(g_YR, g_R);

            value <- c(t(g_YR)%*%B_YR%*%g_YR+t(g_R)%*%B_R%*%g_R);
            #value <- t(g)%*%W2%*%g;

            #print(beta);
            #print(value);
            #print("-----------");

            return(value);

        }
    )

  }

  # num_R, num_Y are the number of R and Y, respectively
  num_R <- length(DatR@data[[VarR]]);
  num_Y <- length(DatY@data[[VarY]]);

  # psill_mle, rangle_mle is the starting_value
  psill_mle <- MD.starting.value[1];
  range_mle <- MD.starting.value[2];


  ##### B_YR, B_R are identity matrix
  # B_YR is the positive-definite weighted matrix
  #B_YR <- list();
  #B_YR[[1]] <- W2;
  #B_YR[[1]] <- GAMMA_YR_B;
  #B_YR[[1]] <- diag(sample_cov_YR$np);
  #B_YR[[1]] <- W_YR[[1]];

  #B_R <- list();
  #B_R[[1]] <- W1;
  #B_R[[1]] <- diag(sample_cov_R$np);
  #B_R[[1]] <- GAMMA_R_B;
  #B_R[[1]] <- W_R[[1]];

  ## B_YR, B_R below are the inverse of variance of gamma

  # run a simulation to compute the covariance matrix of g_n(theta), g_m(theta)
  dist_R <- spDists(DatR);
  dist_Y <- spDists(DatY);
  m_Y <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatY@data[[VarY]] )));
  m_R <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR@data[[VarR]] )));
  K_Y <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){psill_mle*exp(-dist_Y[i,j]/range_mle)})});
  K_R <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){psill_mle*exp(-dist_R[i,j]/range_mle)})});

  J = 300;
  gamma_R_s <- NULL;
  gamma_YR_s <- NULL;
  # R_s, Y_s are independent draws sampling from the Gaussian Random Field
  for(j in 1:J){

    DatR$R_s <- t(mvtnorm::rmvnorm(1, mean = c(t(m_R)), sigma = K_R));
    DatY$Y_s <- t(mvtnorm::rmvnorm(1, mean = c(t(m_Y)), sigma = K_Y));

    g_R_s <- gstat::gstat(NULL, "R_s", R_s~1, DatR);
    g_YR_s <- gstat::gstat(g_R_s, "Y_s", Y_s~1, DatY);

    sample_v_R_s <- gstat::variogram(g_R_s, cutoff = cutoff.R, covariogram = TRUE);
    #sample_v_R_s <- variogram(g_R_s, covariogram = TRUE);
    sample_v_R_s <- sample_v_R_s[order(sample_v_R_s$dist),];


    sample_v_YR_s <- gstat::variogram(g_YR_s, cutoff = cutoff.YR, covariogram = TRUE);
    #sample_v_YR_s <- variogram(g_YR_s, covariogram = TRUE);
    Id_YR <- "R_s.Y_s";
    sample_v_YR_s <- sample_v_YR_s[sample_v_YR_s$id == Id_YR,];
    sample_v_YR_s <- sample_v_YR_s[order(sample_v_YR_s$dist),];

    gamma_R_s <- cbind(gamma_R_s, sample_v_R_s$gamma);
    gamma_YR_s <- cbind(gamma_YR_s, sample_v_YR_s$gamma);

  }
  # B_R is the weighted matrix, which is the inverse of the covariance of g_m(theta)
  B_R <- list();
  gamma_R_sd <- matrix(nrow = num_of_bin_R, ncol = num_of_bin_R, 0)
  for(i in 1:num_of_bin_R){
    for(j in 1:num_of_bin_R){
      gamma_R_sd[i,j] <- cov(gamma_R_s[i, ], gamma_R_s[j, ]);
    }
  }
  B_R[[1]] <- solve(gamma_R_sd);
  ##B_R[[1]] <- diag(sample_cov_R$np);

  ## B_YR is the weighted matrix, which is the inverse of the covariance of g_n(theta)
  B_YR <- list();
  gamma_YR_sd <- matrix(nrow = num_of_bin_YR, ncol = num_of_bin_YR, 0);
  for(i in 1:num_of_bin_YR){
    for(j in 1:num_of_bin_YR){
      gamma_YR_sd[i,j] <- cov(gamma_YR_s[i, ], gamma_YR_s[j, ]);
    }
  }
  B_YR[[1]] <- solve(gamma_YR_sd);
  ##B_YR[[1]] <- diag(sample_cov_YR$np);

  # starting_values are the starting values of minimization of distance function
  starting_values <- MD.starting.value[1:3];
  gs_dat <- list(Dist_YR = Dist_YR[[1]], Gamma_YR = Gamma_YR[[1]], B_YR = B_YR[[1]],
                 Dist_R = Dist_R[[1]], Gamma_R = Gamma_R[[1]], B_R = B_R[[1]]);
  #lower_values <- c(psill = 3, range = 0.2, beta = -1);
  #upper_values <- c(psill = 5, range = 100, beta = -0.0001);

  #lower_values <- c(psill = 1, range = 1, beta = -5);
  #upper_values <- c(psill = 20, range = 25, beta = 5);

  ## Unconstrained global optimization using BFGS
  MD_Results <- optimx::optimx(starting_values, fn = function_lists[[1]],
                #lower = lower_values, upper = upper_values,
                #method = c("L-BFGS-B", "nlminb","bobyqa"),
                #method = c("nlminb"),
                method = c("BFGS"),
                gs = gs_dat);

  # get the point estimate theta0
  psill0 <- MD_Results$psill;
  range0 <- MD_Results$range;
  beta0 <- MD_Results$beta;

  #####################################################################
  ###########################Inference#################################
  #####################################################################

  #Gama_n_theta0 is the partial derivative matrix of g_n(theta) located at theta0
  Gamma_n_theta0 <- matrix(nrow = num_of_bin_YR, ncol = 3, 0);
  for(i in 1:(num_of_bin_YR)){

    Gamma_n_theta0[i,1] <- psill0*exp(-Dist_YR[[1]][i]/range0);
    Gamma_n_theta0[i,2] <- beta0*exp(-Dist_YR[[1]][i]/range0);
    Gamma_n_theta0[i,3] <- beta0*psill0*exp(-Dist_YR[[1]][i]/range0)*(Dist_YR[[1]][i]/(range0^2));

  }

  #GAMMA1 <<- Gamma_n_theta0;

  # Gama_m_theta0 is the partial derivative matrix of g_m(theta) located at theta0
  Gamma_m_theta0 <- matrix(nrow = num_of_bin_R, ncol = 3, 0);
  for(i in 1:(num_of_bin_R)){

    Gamma_m_theta0[i,1] <- 0;
    Gamma_m_theta0[i,2] <- exp(-Dist_R[[1]][i]/range0);
    Gamma_m_theta0[i,3] <- psill0*exp(-Dist_R[[1]][i]/range0)*(Dist_R[[1]][i]/(range0^2));

  }

  #GAMMA2 <<- Gamma_m_theta0;

  A_theta_0 <- t(Gamma_n_theta0)%*%B_YR[[1]]%*%(Gamma_n_theta0) + t(Gamma_m_theta0)%*%B_R[[1]]%*%(Gamma_m_theta0);
  A_theta_0 <- solve(A_theta_0);

 # print(A_theta_0);

  ################################################

  ##### Computation of Sigma_RR
  # to compute the integral of f(x) over R_0, f(x) is the spatial distribution of observation R
  #Sampling_Region_R <<- owin(xrange=c(min(DatR@coords[,1])-1, max(DatR@coords[,1])+1),
  #                          yrange=c(min(DatR@coords[,2])-1, max(DatR@coords[,2])+1));
  #Location_R <<- ppp(x =DatR@coords[,1], y =DatR@coords[,2], window = Sampling_Region_R);
  #Density_R <<- density(Location_R);
  #Prob_Density_R <- Density_R$v/sum(Density_R$v);
  #Integral_R <- sum(Prob_Density_R^2);
  #print(Integral_R);
  Location_x <- c(DatR@coords[,1], DatY@coords[,1]);
  Location_y <- c(DatR@coords[,2], DatY@coords[,2]);

  Sampling_Region <- spatstat::owin(xrange=c(min(Location_x)-1, max(Location_x)+1),
                                    yrange=c(min(Location_y)-1, max(Location_y)+1));
  Location <- spatstat::ppp(x =Location_x, y = Location_y, window = Sampling_Region);
  Density <- density(Location);
  Prob_Density <- Density$v/sum(Density$v);
  Integral <- sum(Prob_Density^2);

  #unit_l1 <- distGeo(c(Density$xcol[1],Density$yrow[1]), c(Density$xcol[1],Density$yrow[2]))/1000;
  #unit_l2 <- distGeo(c(Density$xcol[1],Density$yrow[1]), c(Density$xcol[2],Density$yrow[1]))/1000;
  unit_area <- geosphere::distGeo(c(0,Density$xstep), c(0,Density$ystep))/1000;
  Integral <- Integral/unit_area;

  num_obs <- num_R + num_Y;

  # to compute Sigma_RR
  Sigma_RR <- sapply(1:(num_of_bin_R), function(j){sapply(1:(num_of_bin_R), function(i){

    Grid_R <- mvQuad::createNIGrid(dim=2, type="nHe", level=25);
    fun_RR <- function(x){

      D <- abs(Dist_R[[1]][i]-Dist_R[[1]][j])/2;
      2*psill0^2*exp(-1/range0*(sqrt((x[,1]-D)^2+x[,2]^2)+sqrt((x[,1]+D)^2+x[,2]^2)))

    }
    Int_RR <- mvQuad::quadrature(fun_RR, Grid_R);

    (psill0*exp(-Dist_R[[1]][i])/range0)*(psill0*exp(-Dist_R[[1]][j])/range0)+
      psill0*(psill0*exp(-abs(Dist_R[[1]][i]-Dist_R[[1]][j]))/range0)+
      num_obs*Integral*Int_RR
      #Int_RR/(range0^2)

  })});

  #print(Sigma_RR);

  ## Computation of Sigma_YRYR

  # to compute the integral of f(x) over R_0, f(x) is the spatial distribution of observation Y
  #Sampling_Region_Y <- owin(xrange=c(min(DatY@coords[,1])-1, max(DatY@coords[,1])+1),
  #                          yrange=c(min(DatY@coords[,2])-1, max(DatY@coords[,2])+1));
  #Location_Y <- ppp(x =DatY@coords[,1], y =DatY@coords[,2], window = Sampling_Region_Y);
  #Density_Y <- density(Location_Y);
  #Prob_Density_Y <- Density_Y$v/sum(Density_Y$v);
  #Integral_Y <- sum(Prob_Density_Y^2);

  # Modeling the covariance function of error term Sigma_res
  dist_R <- sp::spDists(DatR);
  dist_Y <- sp::spDists(DatY);
  #dist_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){ sqrt((DatR@coords[j,1]-DatY@coords[i,1])^2+(DatR@coords[j,2]-DatY@coords[i,2])^2)  })});
  dist_YR <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){
             geosphere::distGeo(DatY@coords[i, ],DatR@coords[j, ])/1000
             })});
  m1 <- matrix(nrow = num_Y, ncol = 1, rep(mean(DatR@data[[VarR]] )));
  m2 <- matrix(nrow = num_R, ncol = 1, rep(mean(DatR@data[[VarR]] )));
  K1 <- sapply(1:num_R, function(j){sapply(1:num_Y, function(i){
                psill0*exp(-dist_YR[i,j]/range0
        )})});
  K2 <- sapply(1:num_R, function(i){sapply(1:num_R, function(j){
                psill0*exp(-dist_R[i,j]/range0)
        })});

  # Got the predicted rainfall
  DatY$R  <- c(m1+K1%*%solve(K2, DatR[[VarR]]-m2));

  DatY$residuals <- DatY$Y-beta0*DatY$R;
  vgm_res <- gstat::variogram(residuals~1, DatY, cutoff = 100);
  mfit_res <- gstat::fit.variogram(vgm_res, gstat::vgm(psill = 4, variogram.model, range = 10), fit.method = start.value.method);
  psill_res <- mfit_res$psill[1];
  range_res <- mfit_res$range[1];
  #Sigma_res <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){10*exp(-dist_Y[i,j]/10)})});
  Sigma_res <- sapply(1:num_Y, function(j){sapply(1:num_Y, function(i){psill_res*exp(-dist_Y[i,j]/range_res)})});

  # to compute Sigma_YRYR
  Sigma_YRYR <- sapply(1:(num_of_bin_YR), function(j){sapply(1:(num_of_bin_YR), function(i){

    Grid_Y <- mvQuad::createNIGrid(dim=2, type="nHe", level=25);
    fun_YY <- function(x){

      D <- abs(Dist_YR[[1]][i]-Dist_YR[[1]][j])/2;
      2*beta0^2*psill0^2*exp(-1/range0*(sqrt((x[,1]-D)^2+x[,2]^2)+sqrt((x[,1]+D)^2+x[,2]^2)))+
      psill0*psill_res*exp(-1/range0*(sqrt((x[,1])^2+(x[,2])^2))-1/range_res*(sqrt((x[,1]-2*D)^2+(x[,2])^2)))

    }
    Int_YY <- mvQuad::quadrature(fun_YY, Grid_Y);

    (psill0*exp(-Dist_YR[[1]][i])/range0)*(psill0*exp(-Dist_YR[[1]][j])/range0)+
    psill0*(beta0^2*psill0*exp(-abs(Dist_YR[[1]][i]-Dist_YR[[1]][j]))/range0)+
    psill0*(psill_res*exp(-abs(Dist_YR[[1]][i]-Dist_YR[[1]][j]))/range_res)+
    Int_YY*num_obs*Integral
    #Int_YY/(range0^2)

  })});

 # #print(Sigma_YRYR);

  # To comput Sigma_YRRR

  Sigma_YRRR <- sapply(1:(num_of_bin_R), function(j){sapply(1:(num_of_bin_YR), function(i){

    Grid_YR <- mvQuad::createNIGrid(dim=2, type="nHe", level=25);
    fun_YR <- function(x){

      D <- abs(Dist_YR[[1]][i]-Dist_R[[1]][j])/2;
      2*beta0*psill0^2*exp(-1/range0*(sqrt((x[,1]-D)^2+x[,2]^2)+sqrt((x[,1]+D)^2+x[,2]^2)))

    }
    Int_YR <- mvQuad::quadrature(fun_YR, Grid_YR);

    (beta0*psill0*exp(-Dist_YR[[1]][i])/range0)*(psill0*exp(-Dist_R[[1]][j])/range0)+
    psill0*(beta0*psill0*exp(-abs(Dist_YR[[1]][i]-Dist_R[[1]][j]))/range0)+
    Int_YR*num_obs*Integral

  })});

  #print(Sigma_YRRR);

  #num_pair <- max(max(sample_cov_YR$np), max(sample_cov_R$np));

  Sigma <- (A_theta_0%*%t(Gamma_n_theta0)%*%B_YR[[1]]%*%Sigma_YRYR%*%B_YR[[1]]%*%Gamma_n_theta0%*%A_theta_0)/num_obs+
           (A_theta_0%*%t(Gamma_m_theta0)%*%B_R[[1]]%*%t(Sigma_YRRR)%*%B_YR[[1]]%*%Gamma_n_theta0%*%A_theta_0)/num_obs+
           (A_theta_0%*%t(Gamma_n_theta0)%*%B_YR[[1]]%*%Sigma_YRRR%*%B_R[[1]]%*%Gamma_m_theta0%*%A_theta_0)/num_obs+
           (A_theta_0%*%t(Gamma_m_theta0)%*%B_R[[1]]%*%Sigma_RR%*%B_R[[1]]%*%Gamma_m_theta0%*%A_theta_0)/num_obs;

  Output = list(
    # variogram_R = starting_point_variogram,
    num.obs = num_obs,
    vario.par.point.est = c(psill = unname(psill0), range = unname(range0)),
    vario.par.var.mat = Sigma[2:3,2:3],
    md.point.est = unname(beta0),
    md.var.mat = Sigma[1]
  );

  Output$call <- match.call();
  class(Output) <- "MinimumDistance";

  return(Output);

}




#'@describeIn  MinimumDistance \code{summary} method for class "\code{MinimumDistance}".
#'@param object class \code{MinimumDistance} objects.
#'@export
summary.MinimumDistance <- function(object, ...) {

  x <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }

  ### print output
  cat(paste(rep("=", 20+ 4 + 10 + 10 + 25), collapse=""));
  cat("\n")

  cat(format(" "            , width= 20))
  cat(format(" "            , width= 4))
  cat(format("Point"        , width= 10,  justify="right"))
  cat(format("Std."         , width= 10,  justify="right"))
  cat(format(" "            , width= 25,  justify="centre"))
  cat("\n")

  cat(format(" "        , width=20))
  cat(format("n"        , width=4 ,  justify="left"))
  cat(format("Est."     , width=10,  justify="right"))
  cat(format("Error"    , width=10,  justify="right"))
  cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep=""), width=25, justify="centre"))
  cat("\n")

  cat(paste(rep("=", 20+ 4 + 10 + 10 + 25), collapse=""));
  cat("\n")

  cat(format("Minimum Distance", justify = "left"));
  cat("\n")

    CI_l <- x$md.point.est-qnorm(1-alpha/2)*sqrt(x$md.var.mat);
    CI_r <- x$md.point.est+qnorm(1-alpha/2)*sqrt(x$md.var.mat);
    cat(format("R"  , width= 20,   justify = "left"))
    cat(format(sprintf("%3.0f", x$num.obs),
               width= 4,   justify = "left"))
    cat(format(sprintf("%3.3f", x$md.point.est),
               width= 10,  justify="right"))
    cat(format(sprintf("%3.3f", sqrt(x$md.var.mat)),
               width= 10,  justify="right"))
    cat(format(paste("[", sprintf("%3.3f", CI_l), " , ", sep=""),
               width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_r), "]", sep=""),
               width=11, justify="left"))
    cat("\n")

  cat(paste(rep("=", 20+ 4 + 10 + 10 + 25), collapse=""));

}
