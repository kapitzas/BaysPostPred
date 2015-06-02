#' Compute Bayesian p-value using the response variables' posterior distribution, or sums of squares.
#' 
#' @param jags.object Valid jags output (e.g. as produced by jags in packages jags or R2jags)
#' @param obs Observations
#' @param sim.obs Simulated observations
#' @param fits optional. Model fits
#' @param sim.fits optional. Simulated fits
#' @param plots optional. Should plots be produced?
#' @param store.output optional. Should plots be stored?
#' @param method \code{method = "ssq"} or \code{method = "post"}. Method used to compute Bayesian p-value. See details
#' @return Both methods return Bayesian p-values of model fit and optional plots. See details.
##' @details The function calculates Bayesian p-values in two alternative ways that can be specified via the parameter \code{method}.
##' \itemize{
##'  \item{\code{method = "post"}} {produces a p-value as measure of model fit from observations and simulated observations, both for the entire model and at node level. The method checks whether the simulated observations produced by the model are distributed around the actually observed value in such a way that the observed value approximates the median of the simulated observations. The if this is the case, the returned p-value would be at around 0.5., indicating good model fit.}
##'  \item{\code{method = "ssq"}} {looks at the sums of squares (SSQ) of the residuals between observations and fits and the SSQ of the residuals between their respective simulations. The p-value is calculated as the proportion of simulation SSQ that are larger than the SSQ of actual observations and fits (see Kéry and Schaub 2012).}
##' }
##'
#' @examples
#' Both examples use an adoption of the housemartin state-space-model by Kéry and Schaub (2012).
#'  
#' bpvalue(housemartin.ssm2, obs = "y", sim.obs = "S.y", store.output = F, method = "post")
#' 
#' bpvalue(housemartin.ssm2, obs = "y", sim.obs = "S.y", fits = "logN.est", sim.fits = "S
#' .logN.est", store.output = F, method = "ssq") 
#' @references
#' Kéry, M. and Schaub, M. (2012). Bayesian population analysis using WinBUGS: a hierarchical perspective. Academic Press.
#' @author Simon Kapitza \email{simon.kapitza@@uranus.uni-freiburg.de}
#' @export

bpvalue <- function(jags.object, obs, sim.obs, fits = NULL, sim.fits = NULL, plots = T, store.output = F, path = NULL, method = c("post", "ssq")){
  #back up initial options and parameters and restore on exit
  old_options <- options() 
  on.exit(options(old_options))
  options(error = stop)
  #method post
  if(method == "post"){
    #extract markov chains
    sims.list <- jags.object[[2]][[8]]
    #check if fits and sim fits were specified accidentially and give message
    if(!is.null(fits) == T || !is.null(sim.fits) == T){
      message("fits and sim.fits not required for this method. Will be ommitted")
    }
    #check if specified observations exist in jags output
    if (length(which(names(sims.list) == obs)) == 0){
      stop("observations not found", call. = F)
    }
    #check if specified simulated observations exist in jags output
    if (length(which(names(sims.list) == sim.obs)) == 0){
      stop("simulated observations not found", call. = F)
    }
    #extract specified sims and obs chains from sims list
    a <- which(names(sims.list) == obs)
    b <- which(names(sims.list) == sim.obs)
    observed <- sims.list[a][[1]]
    sim.observed <- sims.list[b][[1]]
    #compute Bayesian p-value
    quantiles <- vector()
    bays.pvalues <- vector()
    for (i in 1:ncol(sim.observed)){
      quantiles[i] <- ecdf(sim.observed[,i])(mean(observed[,i]))
      bays.pvalues[i] <- 1-quantiles[i]
    }
    #plots?
    if (plots == T){
      #path?
      if (store.output == T){
        if (is.null(path)){
          output.path <- getwd()
        }else{
          if(!file.exists(path)){
            stop("Target folder does not exist")
          }
          output.path <- path
        }
        #if plots are stored, open png device
        png(paste0(output.path,"/","post_bpv.png"), width = 2000, height = (ncol(sim.observed)/4)*500)
        par(mfrow=c(ceiling(ncol(sim.observed)/4),4))
      }
      #plot node-level posterior distributions with p-values
      for (i in 1:ncol(sim.observed)){
        plot(ecdf(sim.observed[,i]), xlim=c(min(sim.observed[,i]),max(sim.observed[,i])),
             main= paste0(paste("Node ", i)," (p-value: ",round(bays.pvalues[i],3),")"), 
             xlab="observed", cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
        abline(h=quantiles[i], col="red", lty=2)
      }
      if (store.output == T){
        graphics.off()
      }
    }
    #compute average p-value
    average.pvalue <- mean(bays.pvalues)
    bpoutput <- list(average.pvalue, bays.pvalues)
    names(bpoutput) <- c("average p-value", "p-values at node level")
    return(bpoutput)
  }else{
    #method ssq
    if(method == "ssq"){
      #check if fits and sim.fits are specified
      if(is.null(fits) == T || is.null(sim.fits) == T){
        stop("This method requires a vector containing the fitted internal states (fits) and a vector containing simulated fts (sim.fits)")
      }
      #exctract obs, sim.obs, fits and sim.fits
      sims.list <- jags.object[[2]][[8]]
      a <- which(names(sims.list) == obs)
      b <- which(names(sims.list) == sim.obs)
      c <- which(names(sims.list) == fits)
      d <- which(names(sims.list) == sim.fits)
      observed <- sims.list[a][[1]]
      sim.observed <- sims.list[b][[1]]
      fits <- sims.list[c][[1]]
      sim.fits <- sims.list[d][[1]]
      #compute ssq
      residuals <- observed - fits
      sq <- residuals^2
      sim.residuals <- sim.observed - sim.fits
      S.sq <- sim.residuals^2
      ssq.fit <- rowSums(sq)
      ssq.sim <- rowSums(S.sq)
      #plots
      if (plots == T){
        if (store.output == T){
          if (is.null(path)){
            output.path <- getwd()
          }else{
            if(!file.exists(path)){
              stop("Target folder does not exist")
            }
            output.path <- path
          }
          #open png device if plots are stored
          png(paste0(output.path,"/","ssq_bpv.png"), width = 500, height = 500)
        }
        #plot ssq against each other
        plot(ssq.fit, ssq.sim, main="post. pred. check via ssq")
        abline(0,1, col="red")
      }
      if (store.output == T){
        graphics.off()
        message(paste("Output stored in", output.path))
      }
      #compute bayesian p-value from ssq
      ssq.pvalue <- length(which(ssq.sim-ssq.fit >= 0))/length(ssq.fit)
      return(ssq.pvalue)
    }
  }
}
