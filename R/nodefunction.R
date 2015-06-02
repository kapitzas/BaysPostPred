#' Extraction of Markov chains at node level
#' 
#' This function extracts Markov chains of specified parameters at specified nodes and creates optional trace and density plots.
#' 
#' @param jags.object Valid jags output (e.g. as produced by \code{jags} in packages \code{jags} or \code{R2jags})
#' @param parameter Character string containing name of a monitored parameter 
#' @param node optional. If \code{parameter} was monitored at \code{node}-level, \code{node} must be specified.
#' @param plots optional. Should density- and traceplots of \code{parameter} at \code{node} be created?
#' @param store.output optional. Should plots be stored?
#' @param path target path if plots are stored. If \code{store.output = T} and \code{path} is not specified, files will be stored in current working directory.
#' @return Function returns a list containing the Markov chains of specified \code{parameter} at specified \code{node} (if \code{parameter} was monitored at \code{node} level) for further analysis.
#' @details The function can be used to extract Markov chains (e.g. posterior distributions) of all \code{parameters} monitored by \code{jags}. Parameters that were monitored at \code{node}-level (e.g. time series of population), \code{node} must be specified. Other parameters, such as \code{standard deviations}, that are monitored globally, do not require a value for \code{node}. Trace and density plots show the posterior distribution (e.g. all Markov chains) of \code{parameter} at \code{node}.
#'
#'
#' @examples
#' Both examples apply the housemartin state-space-model by Kéry and Schaub (2012).

#' parnode(housemartin.ssm2,"r", 15, plots = T, store.output = T)
#' parnode(housemartin.ssm2,"logN.est", 7, plots = T, store.output = F)

#' @references
#' Kéry, M. and Schaub, M. (2012). Bayesian population analysis using WinBUGS: a hierarchical perspective. Academic Press.
#' @author Simon Kapitza \email{simon.kapitza@@uranus.uni-freiburg.de}
#' @export


parnode <- function(jags.object, parameter, node = NULL, plots = T, store.output = F, path = NULL){
  old_options <- options()
  old_par <- par(no.readonly = TRUE) 
  on.exit(options(old_options), add = T)
  on.exit(par(old_par), add = T) 
  #suppress browser
  options(error = stop)
  #test, if jags.object exists
  if (!exists(as.character(substitute(jags.object)))){
    stop("Specified jags.object could not be found", call. = F)
  }
  j <-jags.object
  #test, if jags object is rjags or bugs object, (to convert to mcmc-class)
  if (!is.element(class(j), c("rjags", "bugs"))){
    stop("jags.object must be of class 'rjags' or 'bugs")
  }
  #test if parameter is character, to avoid confusion with numeric node
  if (!is.character(parameter)){
    stop("parameter has to be class 'character'", call. = F)
  }
  #assign internal parameter names
  p <- parameter
  x <- as.mcmc(j)
  a <- dimnames(x[[2]])[[2]]
  l <- path
  n <-node
  #check if parameter was monitored
  if (length(grep(p,a, fixed = T)) == 0){
    stop("Specified parameter not monitored", call. = F)
  }
  #check if parameter was monitored globally and node unnecessary
  if ((length(which(a == p)) > 0)){
    #assign parameter name
    b <- p
    if (!is.null(n)){
      warning("Parameter monitored across all nodes, specified node omitted")
    }
  }else{
    if (is.null(n)){
      stop("Specified parameter requires node")
    }else{
      ni <- paste0("[", n, "]")
      if (length(grep(n, a, fixed = T)) == 0){
        stop("Specified node not monitored", call. = F)
      }
      #assign parameter name
      b <- paste0(p,ni)
    }
  }
  #extract markov chains of specified parameter (and node)
  x <- x[,which(dimnames(x[[2]])[[2]] == b)]
  #plot?
  if(plots == T){
    #store?
    if(store.output == T){
      #path?
      if (is.null(l)){
        s <- getwd()
      }else{
        if(!file.exists(l)){
          stop("Target folder does not exist")
        }
        #assign target path
        s <- l
      }
      #output file names
      t <- as.character(paste0("traceplot_",p,n))
      h <- as.character(paste0("densplot_",p,n))
      #open png device
      png(paste0(s, "/", t, ".png"), 800, 800)
    }
    #plot traceplot
    par(mar = c(5.2, 4.5, 4.2, 2.2))
    graphics::plot(as.matrix(x[[1]]), type="n",
                   main=paste("Traceplot of", b), 
                   xlab="Iteration", 
                   ylab=b,
                   cex.lab=1.7, cex.axis=1.7, cex.main=1.7)
    
    for (i in 1:length(x)){  
      lines(as.matrix(x[[i]]), col=i+1, lwd=1)
    }
    if(store.output == T){
      png(paste0(s, "/", h, ".png"), 800, 800)
    }
    #plot density plot
    d <- as.matrix(x)
    par(mar = c(5.2, 4.5, 4.2, 2.2))
    graphics::plot(density(d),
                   xlab=b, 
                   col="grey", 
                   main=paste("Density plot of", b),
                   cex.lab=1.7, cex.axis=1.7, cex.main=1.7)
    polygon(density(d), col="grey")
    if(store.output==T){
      graphics.off()
      #let them know where dens and trace plot were saved
      message(paste("Trace plot and Density plot of", b, "stored as png files in:", s, "File names:", t, "and", h))
    }
  }
  if(store.output==T){
    #create output with file names
    out <- list(p, n, b, t, h, x)
    names(out) <- c("parameter", "node", "internal parameter name", "traceplot filename", "histogram filename", "markov chains")
  }else{
    #create output without file names
    out <- list(p, n, b, x)
    names(out) <- c("parameter", "node", "internal parameter name", "markov chains")
  }
  return(out)
}


