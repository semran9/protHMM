#' @title hmm_ac
#'
#' @description Calculated auto-covariance for profile hidden markov model outputs from HHblits
#'
#' @param hmm path of an HMM file
#' @param lg lag number
#'
#' @return A vector of length 20 * the lag value; by default this is a vector of length 80.
#'
#' @importFrom utils read.table

hmm_ac<- function(hmm, lg = 4){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[which_emmission])[,3:22])
  x[x == "*"]<- 0
  x[]<- 2^-((0.001)*as.numeric(x))
  x[x == 1]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  out_m<- matrix(0, nrow = lg, ncol = 20)
  for(m in 1:lg){
    for(n in 1:20){
      sum = 0
      for(k in 1:(nrow(x)-lg)){
        sum <- sum + (x[k,n] - mean(x[,n])) * (x[k,n] - mean(x[,n]))
      }
      out_m[m,n]<- sum/(nrow(x)-m)
    }
  }
  return(as.vector(out_m))
}
