

# AKSmooth 
##########################################################################

#' AKSmooth, locally kernel smoothing bisulfite sequencing data
#'
#' Modified nadaraya-watson kernel smoother with different kernel functions for the bisulfite sequencing data, 
#' which performs smoothing on low-coverage methylation data chromosome-wise.
#' 
#' @param bsdata  Bisulfite sequencing data containing genomic positions, read coverage and beta methylation values. 
#' @param h The smoothing window size.
#' @param kernel The kernel function used for smoothing.
#' @details \code{bsdata} is a data.frame contains \code{1.)} genomic location, \code{2.)} the chromosomal position, 
#' \code{3.)} M values, describing the number of methylated read covering a single CpG,
#' \code{4.)} Cov values (read coverage), describing the total number of reads covering a single CpG, 
#' \code{5.)} meth values, describing the beta estimated methylation values computed by M/Cov. 
#' 
#' @return The estimated DNA methylation value for a single sample.
#'   
#' @author Junfang Chen, Pavlo Lutsik and Ruslan Akulenko.
#' @export   
#' @examples
#' ## load the bisulfite sequencing data (chromosome 21 for one single sample)
#' data(bn1chr21)  
#' ## define a bandwidth of 30 CpGs
#' fitChr21gau <- AKSmooth(bn1chr21, 30, "Gaussian")  
#' 
AKSmooth <- function(bsdata, h, kernel=c("Gaussian", "Epanechnikov", "Tricube")){
     
  kernel=match.arg(kernel)
  
  if(kernel=="Gaussian")
  {
    
    ksGau <- function(x, y, Cov, h){
        
        n=length(x)   
        fit = vector(,n)
        
        wa <- c()
        wb <- c()  
        
        
        for(j in 1:n){
          
          if (j <= h){  # the meth values close the boundary should be calculted differently
            ind <- j + (-(j-1):h)
            N <- Cov[ind]
            N[-j] <- 0  # set the value of the target point as # of its coverage and all the neighboring points as 0
            wa <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N) * y[ind]
            wb <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N)  # compute the weights of each point within the window
            
            fit[j] <- sum(wa) / sum(wb)
          }
          
          if (j >= h+1 && j <= n-h){
            ind <- j + (-h:h)
            N <- Cov[ind]
            N[-(h+1)] <- 0  # set the value of the target point as # of its coverage and all the neighboring points as 0
            wa <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N) * y[ind]
            wb <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N) # compute the weights of each point within the window
            
            fit[j] <- sum(wa) / sum(wb)   
          }
          
          if (j > n-h){
            ind <- j + (-h:(n-j))
            N <- Cov[ind]
            N[-(h+1)] <- 0 # set the value of the target point as # of its coverage and all the neighboring points as 0
            wa <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N) * y[ind]
            wb <- exp(-(x[j] - x[ind])^2 / (2*h^2) + N) # compute the weights of each point within the window
            
            fit[j] <- sum(wa) / sum(wb)
          }
        }
        return(fit)
      }
    
    chromo <- as.character(unique(data$chr))
    ksMeth <- matrix(nrow=0, ncol=1)
    
    for (i in 1:length(chromo)){
      
      wh <- which(data$chr == chromo[i])
      meth <- data$meth[wh]
      meth <- replace(meth, is.nan(meth), 0.00)  
      Cov <- data$Cov[wh]
          
      chrwiseMeth <- ksGau(1:length(wh), meth, Cov, h)
      chrwiseMeth <- round(chrwiseMeth, 2)
      chrwiseMeth <- as.matrix(chrwiseMeth)
      ksMeth <- rbind(ksMeth, chrwiseMeth)
      
    } 
    return(ksMeth)  
  }
  
  
  if(kernel=="Epanechnikov")
  {
    
    ksEpan <- function(x, y, Cov, h){
      # the total number of CpGs in a whole chr 
      n=length(x) 
      fit = vector(,n) 
      wa <- c()
      wb <- c()   
      for(j in 1:n){
        # CpGs on the boundary with an unsymmetric window mask 
        if (j <= h){  # window mask indices
          ind <- j + (-(j-1):h)
          N <- Cov[ind]
          # define the coverage of the neighboring loci as 1  
          N[-j] <- 1   
          # bi-weights (kernel and coverage weight) are assigned to each CpG within the window 
          wa <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N) * y[ind]
          wb <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N)  
          # return AKSmoothed result 
          fit[j] <- sum(wa) / sum(wb)
        } 
        
        if (j >= h+1 && j <= n-h){
          ind <- j + (-h:h) # window mask indices
          N <- Cov[ind]
          N[-(h+1)] <- 1   
          # assign the bi-weights to each CpG within the window  
          wa <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N) * y[ind]
          wb <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N)   
          
          fit[j] <- sum(wa) / sum(wb)   
        }
        # CpGs on the boundary with an unsymmetric window mask 
        if (j > n-h){
          ind <- j + (-h:(n-j))
          N <- Cov[ind]
          # define the coverage of the neighboring loci as 1 
          N[-(h+1)] <- 1  
          # assign the bi-weights to each CpG within the window   
          wa <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N) * y[ind]
          wb <- 3/4 * (( 1-((x[j] - x[ind])/h)^2 ) * N)   
          
          fit[j] <- sum(wa) / sum(wb)
        }
      }
      return(fit)
    }
  
    chromo <- as.character(unique(data$chr))
    ksMeth <- matrix(nrow=0, ncol=1)
    
    for (i in 1:length(chromo)){
      
      wh <- which(data$chr == chromo[i])
      meth <- data$meth[wh]
      meth <- replace(meth, is.nan(meth), 0.00)  
      Cov <- data$Cov[wh]
      
      chrwiseMeth <- ksEpan(1:length(wh), meth, Cov, h)
      chrwiseMeth <- round(chrwiseMeth, 2)
      chrwiseMeth <- as.matrix(chrwiseMeth)
      ksMeth <- rbind(ksMeth, chrwiseMeth)
      
    } 
    return(ksMeth)  
  }
  
  
  if(kernel=="Tricube")
  {
    ksTricube <- function(x, y, Cov, h){
      # the total number of CpGs in a whole chr  
      n=length(x) 
      fit = vector(,n) 
      wa <- c()
      wb <- c()  
      for(j in 1:n){
        # CpGs on the boundary with an unsymmetric window mask  
        if (j <= h){   
          ind <- j + (-(j-1):h) # window mask indices
          N <- Cov[ind]
          # define the coverage of the neighboring loci as 1   
          N[-j] <- 1   
          # bi-weights (kernel and coverage weight) are assigned to each CpG within the window  
          wa <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N) * y[ind]
          wb <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N)   
          # return AKSmoothed result  
          fit[j] <- sum(wa) / sum(wb)
        }
        
        if (j >= h+1 && j <= n-h){
          ind <- j + (-h:h) # window mask indices
          N <- Cov[ind]
          # assign the bi-weights to each CpG within the window   
          N[-(h+1)] <- 1   
          wa <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N) * y[ind]
          wb <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N)  
          
          fit[j] <- sum(wa) / sum(wb)   
        }
        # CpGs on the boundary with an unsymmetric window mask  
        if (j > n-h){
          ind <- j + (-h:(n-j))
          N <- Cov[ind]
          # define the coverage of the neighboring loci as 1  
          N[-(h+1)] <- 1  
          # assign the bi-weights to each CpG within the window    
          wa <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N) * y[ind]
          wb <- 70/81 * (( 1-abs(((x[j] - x[ind])/h))^3 )^3 * N)   
          
          fit[j] <- sum(wa) / sum(wb)
        }
      }
      return(fit)
    }

    chromo <- as.character(unique(data$chr))
    ksMeth <- matrix(nrow=0, ncol=1)
    
    for (i in 1:length(chromo)){
      
      wh <- which(data$chr == chromo[i])
      meth <- data$meth[wh]
      meth <- replace(meth, is.nan(meth), 0.00)  
      Cov <- data$Cov[wh]
      
      chrwiseMeth <- ksTricube(1:length(wh), meth, Cov, h)
      chrwiseMeth <- round(chrwiseMeth, 2)
      chrwiseMeth <- as.matrix(chrwiseMeth)
      ksMeth <- rbind(ksMeth, chrwiseMeth)
    
    } 
    return(ksMeth)  
  }
  
  
}



 