library(gstat)
library(sp)
library(rgdal)
library(rARPACK)
library(geosphere)
library(bbmle)
library(npsp)
library(MisalignFuncR)

TSBS_Isserlis <- list();
#TSBS_MLE <- list();
#Min_Dist <- list();
#Naive_KNN <- list();
True_Reg <- list();
Min_ABC_MH <- list();
Min_ABC_MH_CH <- list();
#OSQMLE_Isserlis <- list();
#OSMLE <- list();
#KIR <- list();

setwd("C:/Users/Administrator/Desktop/CRAN_misaligned_regression/MinDistance");
rivers <- read.csv("Rivers.csv");
rivers <- rivers[which(rivers$FOR_NLCD<100),]

#cl <- makePSOCKcluster(8);
#registerDoParallel(cl);

#Output <- foreach(SI = 1:8, 
#                  .packages = c("MASS","MisalignFuncR","gstat","sp"),
#                  .combine = c)  %dopar%  {
  
  #print(SI);

  # loading rivers data
  
  train_set <- sample(1:558,277);
  test_set <- setdiff(1:558,train_set);
  rivers_train <- rivers[train_set, ];
  rivers_test <- rivers[test_set, ];

  DatR <- rivers_train[,c("FOR_NLCD", "LAT_DD", "LON_DD")];
  DatR$X <- log((DatR$FOR_NLCD)/(100-DatR$FOR_NLCD));
  coordinates(DatR) <- ~LON_DD+LAT_DD;
  proj4string(DatR) =  "+proj=longlat +datum=WGS84";

  DatY <- rivers_test[, c("CL","LAT_DD","LON_DD")];
  DatY$Y <- log(DatY$CL);
  coordinates(DatY) <- ~LON_DD+LAT_DD;
  proj4string(DatY) =  "+proj=longlat +datum=WGS84";
  
  #kNN_Results <- kNNandRegress(DatR, "X", DatY, "Y");
  #Naive_KNN[[SI]] <- kNN_Results;
  
  Bootstrap_Results <- TwoStepBootstrap(DatR, "X", DatY, "Y", "Exp", FALSE, FALSE,
                                        cutoff = 295, cutoff_u = 40);
  
  MD_Starting_Value <- c(psill = Bootstrap_Results$theta_mean[1], 
                         range = Bootstrap_Results$theta_mean[2], 
                         beta = Bootstrap_Results$point_estimate[2],
                         beta_sd = sd(Bootstrap_Results$bootstrap_results[2,]));
  
  t1 <- Sys.time();
  OSMLE_Results <- OneStepQMLE(DatR, "X", DatY, "Y", "Exp", 
                               TSBS_Starting_Value = MD_Starting_Value,
                               FALSE, FALSE,
                               cutoffR = 40, cutoffYR = 40, cutoffY = 50,
                               Residual = Bootstrap_Results$ols_model$residuals);
  t2 <- Sys.time();
  print(t2-t1);
  #OSMLE[[SI]] <- OSMLE_Results;
  
  #OSMLE_Results <- OneStepQMLE(DatR, "X", DatY, "Y", "Exp", 
  #                             TSBS_Starting_Value = MD_Starting_Value, 
  #                             TRUE, FALSE);
  
  #OSQMLE_Isserlis[[SI]] <- OSMLE_Results;
  
  #TSBS_Isserlis[[SI]] <- Bootstrap_Results;
  #save(DatR, DatY, Bootstrap_Results, V, file = "ABCpre.Rdata");
  
  #V <- DisVecAVSim(DatR = list(DatR), VarR = "X", 
  #                 DatY = list(DatY), VarY = "Y", 
  #                 cutoffR = 150, cutoffYR = 70,
  #                 beta = c(Bootstrap_Results$point[2]),
  #                 Grouped_vgm_model = c(1),
  #                 Grouped_theta1_mle = c(Bootstrap_Results$theta_mean[1]),
  #                 Grouped_theta2_mle = c(Bootstrap_Results$theta_mean[2]),
  #                 Grouped_nuget_mle = c(0),
  #                 V_QMLE = Bootstrap_Results$theta_cov);
  
  
  
  #ABC_Results_1 <- ABCMHEstimation(DatR = list(DatR), VarR = "X", 
  #                               DatY = list(DatY), VarY = "Y", 
  #                               epsilon = 70, W_mat = solve(V),
  #                               cutoffR = 150, cutoffYR = 70,
  #                               TSBS_point_estimation = c(Bootstrap_Results$point[2]), 
  #                               TSBS_sd = sd(Bootstrap_Results$bootstrap_results[2,]),
  #                               Grouped_vgm_model = c(1),
  #                               Grouped_theta1_mle = c(Bootstrap_Results$theta_mean[1]),
  #                               Grouped_theta2_mle = c(Bootstrap_Results$theta_mean[2]),
  #                               Grouped_nuget_mle = c(0),
  #                               V_QMLE = Bootstrap_Results$theta_cov);
  #
  #ABC_Result_2 <- ABCMHCHEstimation(DatR = list(DatR), VarR = "X", 
  #                                  DatY = list(DatY), VarY = "Y", 
  #                                  epsilon = 70, W_mat = solve(V),
  #                                  cutoffR = 150, cutoffYR = 70,
  #                                  TSBS_point_estimation = c(Bootstrap_Results$point[2]), 
  #                                  TSBS_sd = sd(Bootstrap_Results$bootstrap_results[2,]),
  #                                  Grouped_vgm_model = c(1),
  #                                  Grouped_theta1_mle = c(Bootstrap_Results$theta_mean[1]),
  #                                  Grouped_theta2_mle = c(Bootstrap_Results$theta_mean[2]),
  #                                  Grouped_nuget_mle = c(0),
  #                                  V_QMLE = Bootstrap_Results$theta_cov);
  #
  #MD_Starting_Value <- c(psill = Bootstrap_Results$theta_mean[1],
  #                     range = Bootstrap_Results$theta_mean[2],
  #                     beta = Bootstrap_Results$point_estimate[2]
  #                     #beta_std = Bootstrap_Results$bootstrap_results$bootstrap_std[2]
  #                     );

  
  
  
  #Robust_Results <- MinimumDistance(DatR, "X", DatY, "Y", "Exp", MD_Starting_Value,
  #                                  cutoffR = 150, cutoffYR = 70);
  #Min_Dist[[SI]] <- Robust_Results;
  
  DatY$TrueR <- rivers_test[, c("FOR_NLCD")];
  DatY$TrueR <- log(DatY$TrueR/(100-DatY$TrueR));
  True_beta <- lm(DatY$Y ~ DatY$TrueR);
  
  #output <- list(Bootstrap_Results, ABC_Results_1, ABC_Result_2, True_beta);
  #return(output);
  
#}

#for(SI in 1:8){
  
#  TSBS_Isserlis[[SI]] <- Output[[4*SI-3]];
#  Min_ABC_MH[[SI]] <- Output[[4*SI-2]];
#  Min_ABC_MH_CH[[SI]] <- Output[[4*SI-1]];
#  True_Reg[[SI]] <- Output[[4*SI]];
#}
