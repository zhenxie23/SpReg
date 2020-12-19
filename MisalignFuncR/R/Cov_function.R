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