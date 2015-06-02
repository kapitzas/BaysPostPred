#' Plots of jags posterior distributions
#' 
#' This function plots posterior distributions of fits against observations.
#' 
#' @param jags.object Valid jags output (e.g. as produced by \code{jags} in packages \code{jags} or \code{R2jags})
#' @param obs Observations to be plotted
#' @param fits Fits to be plotted against observations. Also model simulations
#' @param store.output optional. Should plots be stored?
#' @param path target path if plots are stored. If \code{store.output = T} and no path is specified, files will be stored in current working directory.
#' @param filename intended output file name as character string without file extension
#' @return Function also returns a list containing the specified posterior distributions for further analysis.
#' @details The function can be used to plot posterior distributions produced by \code{jags} against observations. This helps to visualise systematic misfits of the modelled values compared to actual observations and supports decisions regarding further analysis. If \code{obs} are distributions, their means will be used as values against which \code{fits} are plotted. Plotted are 95% confidence intervals and the means of the specified posterior distributions, as well as a 1:1 line for easy comparison of posterior distributions and observed variables.
#'
#'
#' @examples
#' Both examples use an adoption of the housemartin state-space-model by Kéry and Schaub (2012).

#' fitplot(housemartin.ssm2, obs = "y", fits = "logN.est")
#' fitplot(housemartin.ssm2, obs = "y", fits = "S.logN.est", store.plot = T, path = getwd(), filename = "obsvssims")

#' @references
#' Kéry, M. and Schaub, M. (2012). Bayesian population analysis using WinBUGS: a hierarchical perspective. Academic Press.
#' @author Simon Kapitza \email{simon.kapitza@@uranus.uni-freiburg.de}
#' @export
fitplot <- function(jags.object,obs, fits = NULL, store.plot = F, path = NULL, filename = NULL){
  #store initial options and parameters and restore on exit
  old_options <- options()
  old_par <- par(no.readonly = TRUE) 
  on.exit(options(old_options))
  on.exit(par(old_par), add = T) 
  #suppress browser
  options(error = stop)
  #extract Markov chains
  sims.list <- jags.object[[2]][[8]]
  #check if obs exist
  if (length(which(names(sims.list) == obs)) == 0){
    stop("observations not found", call. = F)
  }
  #check if fits exist
  if (length(which(names(sims.list) == fits)) == 0){
    stop("fits not found", call. = F)
  }
  #extract obs and fits
  x <-as.vector(unique(obs))
  a <- which(names(sims.list) == obs)
  o <- sims.list[a][[1]]
  z <- which(names(sims.list) == fits)
  n <- filename
  f <- sims.list[z][[1]]
  y <- f
  #if plots are stored, is filename specified?
  if (store.plot == T){
    if(is.null(n)){
      stop("specify filename")
    }
    #determine which path to use
    if (is.null(path)){
      output.path <- getwd()
    }else{
      if(!file.exists(path)){
        stop("Target folder does not exist")
      }
      output.path <- path
    }
    #open png device, if plots are stored
    png(paste0(output.path,"/",n,".png"))
  }
  #compute mean, and 95% CI across time series
  quants <-t(apply(y, 2, quantile, c(0.025, 0.5, 0.975)))
  #compute the means of obs at each node (necessary, if there are multiple obs at each node)
  x <- colMeans(o)
  #bring data into correct order for plotting
  data <- data.frame(x, quants)
  data <- data[order(data$x),]
  x <- data[,1]
  data <- data[,-1]
  #plot fits against observations
  plot(x, data[,1], 
       type = "n",
       xlim = c(min(x), max(x)),
       ylim = c(min(data[,1]), max(data[,3])),
       xlab = obs,
       ylab = fits)
  matlines(x, data, lwd = 1, col="black", lty=c(2,1,2))
  abline(0,1, lty = 1, col = "red")
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 1, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend("top", legend = c("95% CI", "Mean", "1:1-line"), 
         lty = c(2,1,1), 
         col = c("black", "black", "red"), 
         horiz = T)
  if(store.plot==T){
    graphics.off()
    #let user know, where their stuff was stored
    message(paste("Fit plot of", obs, "vs", fits, "stored in: ", output.path, "File name:", paste0(filename, ".png")))
  }
  data <- list(o, f)
  names(data) <- c(obs, fits)
  return(data)
}

