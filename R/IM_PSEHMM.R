#' @title IM_psehmm
#' @description The first twenty numbers of this feature correspond to the means of each column of the HMM matrix
#' \eqn{H}. The rest of the features in the feature vector are found in matrix \eqn{T[i,j]}, where \eqn{T[i,j] =
#' \frac{1}{L-i}\sum_{n = 1}^{20-i} [H_{m,n}-H_{m, n+i}]^2, m = 1:L,\space i = 1:d\space and\space j = 1:20}.
#' @param hmm The name of a profile hidden markov model file.
#' @param d The maximum distance between residues column-wise.
#' @note
#' d must be less than 20.
#' @return A vector of length \eqn{20+20\times d-d\times\frac{d+1}{2}}
#' @references Ruan, X., Zhou, D., Nie, R., & Guo, Y. (2020).
#' Predictions of Apoptosis Proteins by Integrating Different Features Based on Improving Pseudo-Position-Specific Scoring Matrix.
#' BioMed Research International, 2020, 1–13.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- IM_psehmm(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
IM_psehmm<- function(hmm, d = 13){
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
  out_mean<- vector('integer',20)
  for(i in 1:20){
    out_mean[i]<- mean(x[, i])
  }
  out_v<- vector('integer', 20*d-d*(d+1)/2)
  count = 1
  for(s in 1:(d)){
    for(n in 1:(20-s)){
      diff<- sum((x[, n] - x[, n+s])^2)
      out_v[count]<- diff
      count<- count + 1
    }
  }
  out_v<- c(out_mean, out_v)
  return(out_v)
}
