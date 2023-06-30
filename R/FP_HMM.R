#' @title fp_hmm
#' @description This feature consists of two vectors, \eqn{d, s}. Vector \eqn{d} corresponds to the sums across the sequence
#' for each of the 20 amino acid columns. Vector \eqn{s} corresponds to a flattened matrix \eqn{S[i, j] = \sum_{k = 1}^{L}
#' H[k, j] \times \delta[k, i]} in which \eqn{\delta[k, i] = 1} when \eqn{A_i = H[k, j]}. \eqn{A} refers to a list of
#' all possible amino acids, \eqn{i, j} span from \eqn{1:20}.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 20.
#' @return A vector of length 400.
#' @references Zahiri, J., Yaghoubi, O., Mohammad-Noori, M., Ebrahimpour, R., & Masoudi-Nejad, A. (2013).
#' PPIevo: Protein–protein interaction prediction from PSSM based evolutionary information.
#' Genomics, 102(4), 237–242.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- fp_hmm(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
fp_hmm<- function(hmm){
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
  f<-x
  f[is.na(f)]<- 0
  f[f<0]<- 0
  seq<- as.matrix(read.table(text = text[emission])[,1])
  d<-vector(length = 20)
  v<-c("A" ,"C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  for(i in 1:20){
    d[i] = sum(f[, i])
  }
  s<- matrix(0, 20, 20)
  for(m in 1:20){
    for(n in 1:20){
      s[m,n]<- sum(f[,n][which(seq==v[m])])
    }
  }
  return(list(as.vector(d), as.vector(s)))
}
