# Estimate population size using camera trapping data


# plot camera trapping results
#############################################################################################################################
plot.camtrap = function(x, circle.size=.2, point.scatter=5){
  sum = aggregate(x$Group_size, by = list(x$Lat, x$Lon), sum)
  colnames(sum) = c('Lat','Lon','Count')
  plot(x$Lon, x$Lat, pch=1, cex = 1, xlab='Longitude', ylab='Latitude')
  points(sum$Lon, sum$Lat, pch = 16,
         col = colorRampPalette(c("grey90", "grey50"))(length(sum$Count))[round(rank(sum$Count))], cex=sum$Count*circle.size)
  X = x[x$Group_size > 0, ]
  dates = as.numeric(as.Date(X$Date, origin = "1900-01-01"))
  points(jitter(X$Lon, factor=point.scatter), jitter(X$Lat, factor=point.scatter),
         pch = 1, cex=1, col = colorRampPalette(c("red", "yellow", "green"))(length(dates))[round(rank(dates))])
}
#############################################################################################################################

trap.out = read.csv('D:/ioz/camera_trap/Wild-boar.csv', header=T)
sim.out = read.csv('D:/ioz/camera_trap/sim.out.csv', header=T)


# plot.camtrap(trap.out)
plot.camtrap(trap.out, circle.size = .5, point.scatter = 10)





# simulate camera trapping processes
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
simu.camtrap = function(x, detect = 0.05 / (40000/360), bearing = runif(camera.N, 0, 2*pi), # Detection range = 50m
                        ind = 5, step.N = 5000, step.L = 0.01 / (40000/360),
                        step.V = 0.2, bias = 15/360*2*pi, range = .5,
                        iteration = 1){

  camera = unique(x[,c('Lon','Lat')]) # x=trap.out
  camera = cbind(camera, Count=0)
  camera.N <- nrow(camera)
  out = as.data.frame(matrix(data = 0, nrow = ind*iteration, ncol = camera.N))


  for (ite in 1:iteration){#  iterations
    plot.camtrap(x)
    for (k in 1:ind){
      footchain = data.frame(ID=1:step.N, X=NA, Y=NA) # for plotting footchain
      loc.x = runif(1, min(x$Lon) + 0.25*(max(x$Lon) - min(x$Lon)), max(x$Lon)-0.25*(max(x$Lon) - min(x$Lon))) # random initial location
      loc.y = runif(1, min(x$Lat) + 0.25*(max(x$Lat) - min(x$Lat)), max(x$Lat)-0.25*(max(x$Lat) - min(x$Lat)))
      loc.x.0 = loc.x; loc.y.0 = loc.y

      ### MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  simulating movement  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      theta  <- runif(1, 0, 2*pi) # moving bearing of the first step
      for (i in 1:step.N){   #walking steps
        for (j in 1:camera.N){ # number of cameras
          d = ((loc.x - camera$Lon[j])^2 + (loc.y - camera$Lat[j])^2)^.5
          bear = atan((loc.y - camera$Lat[j])/(loc.x - camera$Lon[j]))
          if(d < detect & abs(bear-bearing[j]) < 50/360*2*pi)   camera$Count[j] <- camera$Count[j] + 1 # 40 degree detetion region
        }
        move = step.L * rnorm(1, 1, step.V) #
        loc.x = loc.x + move*cos(theta) #
        loc.y = loc.y + move*sin(theta)
        theta.f = theta + rnorm(1,0,bias) # move forward
        if (loc.x > loc.x.0)    theta.b <- pi + atan((loc.y - loc.y.0)/(loc.x - loc.x.0))+ rnorm(1,0,bias)#return to origin
        if (loc.x < loc.x.0)    theta.b <-      atan((loc.y - loc.y.0)/(loc.x - loc.x.0))+ rnorm(1,0,bias)#return to origin
        dist = ((loc.x - loc.x.0)^2 + (loc.y - loc.y.0)^2)^.5
        theta = sample(c(theta.f, theta.b), 1, prob=c((1-dist/range), dist/range)) #home range
        footchain$X[i] = loc.x;   footchain$Y[i] = loc.y # for plotting footchain
      }
      # plot footchain
      lines(footchain$X, footchain$Y, col=sample(2:11,1))
      points(loc.x.0,loc.y.0, col='blue',pch=17,cex=.8)
      points(footchain$X[step.N], footchain$Y[step.N], col='red',pch=15,cex=.8)
      # points(camera$Lon, camera$Lat, col='black',pch=16)
      ### MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  simulating movement  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

      camera$Count = sort(camera$Count) # NOT spatial explicit
      out[k+(ite-1)*ind,] = camera$Count            # generate out
      print(paste('Iteration ', ite, ';  Ind. ', k, sep=''))
    }
    camera$Count = camera$Count * 0
  }

  Ind = rep(1:ind, iteration)
  out = cbind(Ind = Ind, out)
  assign('sim.out', out, envir = .GlobalEnv)
  write.csv(out, 'D:/sim.out.csv', row.names = F)
  return(out)
}
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# simu.camtrap(trap.out, ind = 5, iteration = 1)





#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
predict.camtrap = function(simu, x, plot=F){
  library(randomForest)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest needed for this function to work. Please install it.", call. = FALSE)
  }

  RF = randomForest(simu[,c(2:ncol(simu))], simu[,1], prox=TRUE, importance=TRUE, ntree=1000)

  sum = aggregate(x$Group_size, by = list(x$Lat, x$Lon), sum)
  colnames(sum) = c('Lat','Lon','Count')
  obs = sort(sum$Count)
  pred <- predict(RF, obs, type="response", predict.all=TRUE)
  pred.rf.int <- apply( pred$individual, 1, function(x) { quantile(x, c(0.025, 0.5, 0.975) )})
  if(plot){
    plot(density(pred$individual),xlab='Number of individuals', ylab='Frequency',col='darkgrey',xlim=c(0,max(simu$Ind)),lwd=2,
         main=paste('Predicted population size:',round(pred.rf.int[2,1],1), sep=' '))
    abline(v = pred$aggregate, lwd=2)
    abline(v = c(pred.rf.int[1,1], pred.rf.int[3,1]), lwd=1, col='black',lty=2)
  }
  return(round(pred.rf.int, 1))
}
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

# predict.camtrap(sim.out, trap.out, plot=T)
