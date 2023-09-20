myboyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100, 
          PEplot = TRUE, rm.duplicate = TRUE, method = "spearman") 
{
  boycei <- function(interval, obs, fit) {
    # pi <- sum(as.numeric(obs >= interval[nrow(interval),1] & obs <= interval[nrow(interval),2]))/length(obs)
    # ei <- sum(as.numeric(fit >= interval[nrow(interval),1] & fit <= interval[nrow(interval),2]))/length(fit)
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2]))/length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2]))/length(fit)
    return(round(pi/ei, 10)) #rapport de obs et fit dans chaque box
  }
  if (inherits(fit, "RasterLayer")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- extract(fit, obs)
    }
    fit <- getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  mini <- min(fit, obs)
  maxi <- max(fit, obs)
  if (length(nclass) == 1) {
    if (nclass == 0) {
      if (window.w == "default") {
        window.w <- (max(fit) - min(fit))/10
      }
      vec.mov <- seq(from = mini, to = maxi - window.w, #vector of pred values create the boyce bins
                     by = (maxi - mini - window.w)/res)
      # vec.mov[res + 1] <- vec.mov[res + 1] + 1 #add +1 to last value (to be sure to be above 1? )
      interval <- cbind(vec.mov, vec.mov + window.w)
    }
    else {
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - 
                                                     mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else {
    vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < 
                                     mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA
  }
  else {
    r <- 1:length(f)
    if (rm.duplicate == TRUE) {
      r <- c(1:length(f))[f != c(f[-1], TRUE)]
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method) #boyce index
  }
  HS <- apply(interval, 1, sum)/2 #mean of each bin
  if (length(nclass) == 1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1
  }
  HS <- HS[to.keep]
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", 
         col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}