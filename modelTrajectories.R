coral.model <- function(t, y, parms){
#this code reproduces the figure from paper Thresholds and the resilience of Caribbean coral reefs Mumby et al. 2007. 
  # Args:
  #   t: the current time (of the model)
  #   y: a state vector consisting of:
  #      M: Cover of  microalgae 
  #      C: Cover of coral
  #   parms: a named list of parameters. 
  #     a: algal overgrowth rate
  #     r: coral overgrowth rate
  #     g: grazing rate
  #     d: coral mortality
  #     gamma: macroalgae colonization rate 
  #
  # Returns:
  #    a list containing dM / dt and dC / dt
  
  M <- y["M"]
  C <- y["C"]
  
  
  a <- parms$a
  g <- parms$g
  gamma <- parms$gamma
  r <- parms$r
  d <- parms$d
  
  dM.dt <- a*M*C - g*M/(1-C) + gamma*M*(1-C-M)
  dC.dt <- r*(1-M-C)*C -d*C - a*M*C
  
  
  return(list(c(dM.dt,dC.dt))) 
}


require(deSolve) # load the deSolve package into memory if it isn't loaded already




times <- seq(from = 0, to = 30, length = 100) # sequence of times for which we want a solution

mvec = seq(from= 0.05, to=0.9,by = 0.05) 

par(mfrow=c(1,2))

for (g in c(0.3,0.1)){ #loop over two different values of g
  
  test.parms <- list(a=0.01, d=0.5,  g=g, r=1, gamma=0.75)
  
  for(m0 in mvec ){ #loop for different m0
    
    mvec.lim=0.95-m0
    mvec.lim = round(mvec.lim, digits = 2)
    cvec =seq(from=0.05,to=mvec.lim, by = 0.05)
    for(c0 in cvec){ #loop for different c0
      
      y.initial <- c(M =m0, C=c0) # initial density updated at each iteration of the loop
      
      #simulate the trajectory for a given initial condition 
      output <- ode(y     = y.initial,
                    times = times,
                    func  = coral.model,
                    parms = test.parms)
     
      if(m0==0.05 && c0==0.05){ #inital plot
        plot(output[,2],output[,3],type='l',col='red',xlim=c(0,1), xaxs="i",ylim=c(0,1),yaxs="i",xlab="Microalgae",ylab="Coral",main=sprintf("g = %.1f",g),yaxs="r")
        points(m0,c0,col="red",pch=".",cex=3)
      }
      
      else{ #overlaying trajectories from different starting points 
        lines(output[,2],output[,3],col='red')
        points(m0,c0,col="red",pch=".",cex=3)
        
      }
      
      if (g==0.1){ #plotting fixed points
        points(0.86,0.0068107112,col="black",pch=20,cex=2)
        legend(0.86, 0.0068107112+0.1, bquote(M[s]),cex = .8,text.col="black",box.col="white",bg="white")
        
      }
      
      else{ #plotting fixed points for g=0.3
        points(0.006920419,0.4950222,col="black",pch=20,cex=2)
        legend(0.001920419+0.02, 0.4990222+0.08, bquote(C[s]),cex = .8,text.col="black",box.col="white",bg="white")
        points(0.5911754,0.006904276,col="black",pch=20,cex=2)
        legend(0.5911754+0.02, 0.006904276+0.08, bquote(M[s]),cex = .8,text.col="black",box.col="white",bg="white")
      }
        
      
    }
  }
  
  
}








