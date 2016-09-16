calcID <- function(formula, data) {
  y <- attr(terms(formula),"variables")[[2]]
  x <- attr(terms(formula),"variables")[[3]]
  y <- which(names(data) == y)
  x <- which(names(data) == x)
  if(length(x) == 0L | length(x) == 0L) stop("variables in the formula not found in the data set")
  y <- data[,y]
  x <- data[,x]
  id <- 0.5*sum(abs(y/sum(y) - x/sum(x)))
  return(round(id, 3))
}

IDsim <- function(formula, data, total=NULL, weights=NULL, seed=NULL, times=1000) {
  if(!is.null(seed)) set.seed(seed)
  
}


dvals <- function(mlm) {
  
  resids <- residuals(mlm)
  rf <- ranef(mlm)
  results <- matrix(nrow=length(resids), ncol=length(rf)+1)
  rownames(results) <- names(resids)
  colnames(results) <- c("",names(rf))
  
  results[,1] <- resids
  mod.data <- slot(mlm, "frame")
  for(i in 1:length(rf)) {
    
    k <- which(names(mod.data) == names(rf)[i])
    mch <- match(mod.data[,k], rownames(rf[[i]]))
    results[,(i+1)] <- rf[[i]][mch,1]
    
  }
  
  return(results)
  
}


condVar <- function(model) {
  cat("\nCalculating variances, please wait")  
  u0 <- ranef(model, condVar=T)
  return(u0)
}

catplot <- function(model, var=NULL, level=2, method=c("quick","goldstein"), scale=T, labels=F, z=NULL, sigma=NULL, cex=0.7, add=F, ymin=NULL, ymax=NULL) {
  
  alpha <- 0.05
  
  search <- function(zb, se, alpha) {
    gammas <- NULL
    for(i in 1:(length(se)-1)) {
      for(j in (i+1):length(se)) {
        se_i <- se[i]
        se_j <- se[j]
        se_ij <- sqrt(se_i^2 + se_j^2)
        gamma_ij <- 2 * (1 - pnorm(zb * (se_i + se_j) / se_ij))
        gammas <- c(gammas, gamma_ij)
      }
    }
    return(abs(mean(gammas) - alpha))
  }
  
  if(class(model) != "lmerMod") stop("Model object is not of class lmerMod (from the lme4 package)")
  
  nlevels <- length(slot(model, "flist"))
  level <- level - 1
  if(level < 1 | level > nlevels) stop(paste("Please select a level between 2 and",nlevels+1))

  if(is.null(var)) {
    cat("\nCalculating variances, please wait")  
    u0 <- ranef(model, condVar=T)   
  } else {
    u0 <- var
  }
  u0se <- sqrt(attr(u0[[level]], "postVar")[1, , ]) 
  u0 <- u0[[level]][,1]
  
  if(method[1]=="goldstein") {
    cat("\nCalculating average confidence interval, please wait")
    opt <- optimise(search, interval=c(1.39, 1.96), se=u0se, alpha=0.05, tol=0.001)
    z <- opt$minimum
    cat("\nThe z value is",z)
  } else {
    ifelse(is.null(z), z <- 1.39, z <- z)
  }
  
  ci <- data.frame(u0, lwr=u0-z*u0se,  upr=u0+z*u0se)
  ord <- order(ci[,1], decreasing=T)
  ci <- ci[ord,]
  n <- nrow(ci)
  
  names(ci) <- c("u0","lwr","upr")
  if(scale) {
    ifelse(is.null(sigma), ci <- ci/sigma(model), ci <- ci/sigma)
  }
  
  if(!add) {quartz(height=4.3228346, width=4.4173228)
    par(mai=c(0.85,0.90,0.25,0.25))
  }
  
  if(is.null(ymin) & is.null(ymax)) {
    ymin <- min(ci)
    ymax <- max(ci) + 0.5
  }
  
  if(!scale) plot(x=1:n, y=ci$u0, xlim=c(1,n), ylim=c(ymin, ymax), col="white", xlab="Rank",ylab="Residual difference", las=1)
  if(scale) plot(x=1:n, y=ci$u0, xlim=c(1,n), ylim=c(ymin, ymax), col="white", xlab="Rank",ylab="Scaled residual difference", las=1)
  for(i in 1: n) {
    arrows(i, ci$lwr[i], i, ci$upr[i], code=3, angle=90, length=0.025)
  }
  points(x=1:n, y=ci$u0, pch=21, bg="white", cex=0.8)
  abline(h=0, lty="dotted")
  i <- as.integer(rownames(ci))
  if(labels) text(1:n, ci$upr+0.25, labels=rownames(ranef(model)[[level]])[i], cex=cex)
  rownames(ci) <- rownames(ranef(model)[[level]])[i]
  invisible(ci)
}