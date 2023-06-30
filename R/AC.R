#' @title hmm_ac
#' @description This feature calculates the covariance between two residues separated by a lag value within
#' the same amino acid emission frequency column along the protein sequence.
#' @param hmm The name of a profile hidden markov model file.
#' @param lg The lag value, which indicates the distance between residues.
#' @note
#' The lag value must be less than the length of the protein sequence
#' @return A vector of length 20 \eqn{\times} the lag value; by default this is a vector of length 80.
#' @references Dong, Q., Zhou, S., & Guan, J. (2009).
#' A new taxonomy-based protein fold recognition approach based on autocross-covariance transformation.
#' Bioinformatics, 25(20), 2655–2662.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_ac(system.file("extdata", "1DLHA2-7", package="protHMM"))

hmm_ac<- function(hmm, lg = 4){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emmission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emmission])[,3:22])
  x[x == "*"]<- 0
  #x[]<- 2^-((0.001)*as.numeric(x))
  #x[x == 1]<- 0
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
