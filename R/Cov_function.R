
Cov_Def <- function(model, dist, theta1, theta2, nuget, order){
  
  if(order == 1){
    value <- matrix(nrow = 2, ncol = 1, 0);
    
    if(model == 1){
      value[1] <- exp(-dist/theta2);
      value[2] <- (dist/(theta2^2))*theta1*exp(-dist/theta2);
      return(value);
    }
    
    if(model == 2){
      value[1] <- exp(-(dist/theta2)^2);
      value[2] <- (2*dist^2/(theta2^3))*theta1*exp(-(dist/theta2)^2);
    }
    
    if(model == 3){
      if(dist <= theta2){
        value[1] <- 1-dist/theta2;
        value[2] <- theta1*dist/(theta2^2);
      }
    }
    
    if(model == 4){
      if(dist <= theta2){
        value[1] <- 1-(3*dist)/(2*theta2)+(dist^3)/(2*theta2^3);
        value[2] <- theta1*((3*dist)/(2*theta2^2)-1.5*dist^3/(theta2^4));
      }
    }
    
    return(value);
    
  }
  
  if(order == 2){
    
    value <- matrix(nrow = 2, ncol = 2, 0);
    if(model == 1){
      value[1,1] <- 0;
      value[1,2] <- (dist/theta2^2)*exp(-dist/(theta2^2));
      value[2,1] <- value[1,2];
      value[2,2] <- theta1*exp(-dist/theta2)*(((dist^2)/(theta2^4))-2*dist/(theta2^3));
    }
    
    if(model == 2){
      value[1,1] <- 0;
      value[1,2] <- (2*dist^2/(theta2^3))*exp(-(dist/theta2)^2);
      value[2,1] <- value[1,2];
      value[2,2] <- theta1*exp(-(dist/theta2)^2)*(4*dist^4/(theta2^6)-6*dist^2/(theta2^4));
    }
    
    if(model == 3){
      if(dist <= theta2){
        value[1,1] <- 0;
        value[1,2] <- theta1*dist/(theta2^2);
        value[2,1] <- value[1,2];
        value[2,2] <- -2*theta1*dist/(theta2^3);
      }
    }
    
    if(model == 4){
      if(dist <= theta2){
        value[1,1] <- 0;
        value[1,2] <- ((3*dist)/(2*theta2^2)-1.5*dist^3/(theta2^4));
        value[2,1] <- value[1,2];
        value[2,2] <- theta1*((-3*dist)/(theta2^3)+(6*dist^3)/(theta2^5));
      }
    }
    
    return(value);
    
  }
  
}


Semivariogram_function <- function(model, dist, theta1, theta2, nuget, cov){
  
  if(cov == FALSE){
    
    ## "Exp" model
    if(model == 1){
      value <- theta1*(1-exp(-dist/theta2));
      return(value);
    }
    
    ## "Gau" model
    if(model == 2){
      value <- theta1*(1-exp(-(dist/theta2)^2));
      return(value);
    }
    
    ## "Lin" model
    if(model == 3){
      if(dist <= theta2){
        value <- theta1*(dist/theta2);
      }
      else{ 
        value <- theta1*1;
      }
      return(value);
    }
    
    # "Sph" model
    if(model == 4){
      if(dist <= theta2){
        value <- theta1*((3*dist)/(2*theta2)-0.5*(dist/theta2)^3);
      }
      else{ 
        value <- theta1;
      }
      return(value);
    }
    
    # "Log" model
    #if(model == 5){
    #  value <- theta1*log(theta2+dist);
    #  return(value);
    #}
  } 
  
  if(cov == TRUE){
    
    if(dist > 0){
      
      ## "Exp" model
      if(model == 1){
        value <- theta1*(exp(-dist/theta2));
        return(value);
      }
      
      ## "Gau" model
      if(model == 2){
        value <- theta1*(exp(-(dist/theta2)^2));
        return(value);
      }
      
      ## "Lin" model
      if(model == 3){
        if(dist <= theta2){
          value <- theta1*(1-dist/theta2);
        }else{ 
          value <- theta1*(1-1);
        }
        return(value);
      }
      
      # "Sph" model
      if(model == 4){
        if(dist <= theta2){
          value <- theta1*(1-(3*dist)/(2*theta2)+0.5*(dist/theta2)^3);
        }else{ 
          value <- theta1*(1-1);
        }
        return(value);
      }
      
      # "Log" model
      #if(model == 5){
      #  if(theta2)
      #  value <- theta1*(1-log(theta2+dist));
      #  return(value);
      #}
    }
    if(dist == 0){
      
      ## "Exp" model
      if(model == 1){
        
        value <- nuget + theta1*(exp(-dist/theta2));
        
        return(value);
      }
      
      ## "Gau" model
      if(model == 2){
        value <- nuget + theta1*(exp(-(dist/theta2)^2));
        return(value);
      }
      
      ## "Lin" model
      if(model == 3){
        if(dist <= theta2){
          value <- nuget + theta1*(1-dist/theta2);
        }
        else{ 
          value <- nuget + theta1*(1-1);
        }
        return(value);
      }
      
      # "Sph" model
      if(model == 4){
        if(dist <= theta2){
          value <- nuget + theta1*(1-(3*dist)/(2*theta2)+0.5*(dist/theta2)^3);
        }else{ 
          value <- nuget + theta1*(1-1);
        }
        return(value);
      }
      
      # "Log" model
      #if(model == 5){
      #  if(theta2)
      #  value <- theta1*(1-log(theta2+dist));
      #  return(value);
      #}
      
    }
  }
  
}
