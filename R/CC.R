#' @title hmm_cc
#' @description The feature calculates the covariance between different residues separated along the protein sequences
#' by a lag value across different amino acid emission frequency columns.
#' @param hmm The name of a profile hidden markov model file.
#' @param lg The lag value, which indicates the distance between residues.
#' @note
#' The lag value must less than the length of the amino acid sequence.
#' @return A vector of length 20 x 19 x the lag value; by default this is a vector of length 1520.
#' @references Dong, Q., Zhou, S., & Guan, J. (2009).
#' A new taxonomy-based protein fold recognition approach based on autocross-covariance transformation.
#' Bioinformatics, 25(20), 2655–2662.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_cc(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_cc<- function(hmm, lg = 4){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  x[x == "*"]<- 0
  #x[]<- 2^-((0.001)*as.numeric(x))
  #x[x == 1]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  out_v<- vector(length = 20*19*lg)
  count<- 1
  for(z in 1:lg){
    for(m in 1:20){
      for(n in 1:20){
        sum <- 0
        if(m != n){
          for(i in 1:(nrow(x)-z)){
            sum = sum + (x[i, m] - mean(x[, m])) * (x[i + z, n] - mean(x[, n]))
          }
          out_v[count]<- sum/(nrow(x)-z)
          count = count + 1
          sum = 0
        }
      }
    }
  }
  return(out_v)
}
