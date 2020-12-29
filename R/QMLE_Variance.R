QMLE_Variance <- function(W, V, FDV, SDV, Vsample, is_cov_misspecified, is_den_misspecified){

  #library(georob);
  #library(matrixcalc);

  num_of_group <- length(V);
  num_of_par <- length(FDV[[1]]);
  num_of_obs <- nrow(W[[1]]);
  V_QMLE <- list();

  if((is_cov_misspecified == FALSE) && (is_den_misspecified == FALSE)){

    for(k in 1:(num_of_group)){

      Hessian <- sapply(1:num_of_par, function(i){sapply(1:num_of_par, function(j){

        0.5*matrixcalc::matrix.trace(solve(V[[k]],FDV[[k]][[i]])%*%solve(V[[k]],FDV[[k]][[j]]))

      })});

      V_QMLE[[k]] <- solve(Hessian);
      #V_QMLE[[k]] <- round(V_QMLE[[k]],4);
      #V_QMLE[[k]] <- (V_QMLE[[k]]+t(V_QMLE[[k]]))/2;
      #V_QMLE[[k]] <- V_QMLE[[k]]/num_of_obs;
    }
  }

  if((is_cov_misspecified == TRUE) && (is_den_misspecified == FALSE)){

    for(k in 1:(num_of_group)){

      V2 <- solve(V[[k]],Vsample[[k]]);

      Hessian <- sapply(1:num_of_par, function(i){sapply(1:num_of_par, function(j){

        V0 <- solve(V[[k]],FDV[[k]][[i]]);
        V1 <- solve(V[[k]],FDV[[k]][[j]]);

        ans <- -0.5*(matrixcalc::matrix.trace(V0%*%V1))+
        0.5*(matrixcalc::matrix.trace(solve(V[[k]],SDV[[k]][[i]][[j]])))+
        (matrixcalc::matrix.trace(V2%*%V0%*%V1))-
        0.5*(matrixcalc::matrix.trace(V2%*%solve(V[[k]],SDV[[k]][[i]][[j]])));

        return(ans);

      })});


      Info <- sapply(1:num_of_par, function(i){sapply(1:num_of_par, function(j){

        V0 <- solve(V[[k]],FDV[[k]][[i]]);
        V1 <- solve(V[[k]],FDV[[k]][[j]]);

        #-0.25*matrix.trace(V0)*matrix.trace(V1)+
        #0.25*matrix.trace(Vsample[[k]]%*%V0%*%solve(V[[k]]))*matrix.trace(Vsample[[k]]%*%V1%*%solve(V[[k]]))+
        ans <- 0.5*(matrixcalc::matrix.trace(V0%*%V2%*%V1%*%V2));
        return(ans);

      })});

      #print(matrix.trace(V2));

      V_QMLE[[k]] <- solve(Hessian,Info)%*%solve(Hessian);
      #V_QMLE[[k]] <- round(V_QMLE[[k]],4);
      #V_QMLE[[k]] <- (V_QMLE[[k]]+t(V_QMLE[[k]]))/2;
      #V_QMLE[[k]] <- V_QMLE[[k]]/num_of_obs;
    }
  }

  return(V_QMLE);

}
