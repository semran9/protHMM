#' @title pse_hmm
#' @description The first twenty numbers of this feature correspond to the means of each column of the HMM matrix
#' \eqn{H}. The rest of the features in the feature vector are given by correlation of the \eqn{ith} most
#' contiguous values along the chain per each amino acid column, where \eqn{0<i<g+1}. This creates a vector of 20 \eqn{\times}
#' g, and this combines with the first 20 features to create the final feature vector.
#' @param hmm The name of a profile hidden markov model file.
#' @param g The contiguous distance between residues.
#' @note
#' g must be less than the length of the protein sequence
#' @return A vector of length \eqn{20 + g \times 20}, by default this is 320.
#' @references Chou, K., & Shen, H. (2007). MemType-2L: A Web server for predicting membrane proteins and their types
#' by incorporating evolution information through Pse-PSSM. Biochemical and Biophysical Research Communications,
#' 360(2), 339–345.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- pse_hmm(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
pse_hmm<- function(hmm, g = 15){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  x[x == "*"]<- 0
  x[]<- 2^-((0.001)*as.numeric(x))
  x[x == 1]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  out_v <- vector(length = 320)
  for(m in 1:20){
    out_v[m]<- mean(x[, m])
  }
  for(m in 1:20){
    for(l in 1:g){
      for(n in 1:(nrow(x)-l)){
        out_v[20 + m + 20 * (l-1)]<- (1/(nrow(x)-l)) * (x[n, m] - x[n + l, m])^2
      }
    }
  }
  return(out_v)
}
