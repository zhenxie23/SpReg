SampleCovarianceMatrix <- function(Cov_sample, Dist){

  # DatR: a spatial object, see function coordinates()
  # VarR: name of variable R
  # Dist: distance matrix
  # Return value: sample covariance matrix

  num_Y <- nrow(Dist);
  num_R <- ncol(Dist);
  Distvec <- matrixcalc::vec(Dist);
  y <- rep(1:num_R, each = num_Y);
  x <- rep(seq(1:num_Y), num_R);
  Covvec <- data.frame(Distvec, x, y);
  Covvec <- Covvec[order(Covvec$Distvec), ];

  num_of_bin <- nrow(Cov_sample);
  Cov_sample_mat <- matrix(nrow = num_Y, ncol = num_R, 0);

  for(j in 1:Cov_sample$np[1]){
    Cov_sample_mat[Covvec$x[j], Covvec$y[j]] <- Cov_sample$gamma[1];
  }

  INDEX <- Cov_sample$np[1];

  #print(INDEX);

  for(i in 2:num_of_bin){

    while(Covvec$Distvec[INDEX+1] <= Cov_sample$dist[i]){

      INDEX <- INDEX +1;
      Cov_sample_mat[Covvec$x[INDEX], Covvec$y[INDEX]] <- Cov_sample$gamma[i-1]+
      (Cov_sample$gamma[i]-Cov_sample$gamma[i-1])/(Cov_sample$dist[i]-Cov_sample$dist[i-1])*(Covvec$Distvec[INDEX]-Cov_sample$dist[i-1]);

    }

    #print(INDEX);
  }

  return(Cov_sample_mat);

}
