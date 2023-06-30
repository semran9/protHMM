#' @title hmm_svd
#' @description This feature uses singular value decomposition (SVD) to reduce the dimensionality of the inputted hidden
#' markov model matrix. SVD factorizes a matrix C of dimensions \eqn{i, j} to \eqn{U[i, r] \times \Sigma[r, r] \times V[r, j]}.
#' The diagonal values of \eqn{\Sigma} are known as the singular values of matrix C, and are what are returned with this function.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 20.
#' @references Song, X., Chen, Z., Sun, X., You, Z., Li, L., & Zhao, Y. (2018).
#' An Ensemble Classifier with Random Projection for Predicting Protein–Protein Interactions Using Sequence and Evolutionary Information.
#' Applied Sciences, 8(1), 89.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_svd(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'

hmm_svd<- function(hmm){
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
  return(svd(x)[[1]])
}
