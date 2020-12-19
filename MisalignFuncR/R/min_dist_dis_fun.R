### distance function
Dis <- function(par, gs){

  theta1 <- par[1];
  theta2 <- par[2];
  nuget <- par[3];

  vgm_model <- gs$model; beta <- gs$beta; num_of_bin <- gs$Num_of_bins;
  Dist = gs$Dist; Gamma = gs$Gamma; B = gs$B;

  num_of_beta <- length(beta);
  Par_Gamma <- NULL; index = 0;
  for(v in 1:num_of_beta){
    for(i in 1:num_of_bin[v]){

      index <- index + 1;
      semi_vgm <- Semivariogram_function(vgm_model, Dist[index],
                                         theta1, theta2, nuget, cov = TRUE);
      Par_Gamma <- c(Par_Gamma, 0 + beta[v]*semi_vgm);
    }
  }

  for(i in 1:num_of_bin[num_of_beta+1]){

    index <- index + 1;
    semi_vgm <- Semivariogram_function(vgm_model, Dist[index],
                                       theta1, theta2, nuget, cov = TRUE);
    Par_Gamma <- c(Par_Gamma, 0 + semi_vgm);
  }

  dis_vector <- Par_Gamma - Gamma;

  value <- t(dis_vector)%*%B%*%dis_vector;
  return(value);

}
