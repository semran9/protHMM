#' @title hmm_distance
#' @description This feature calculates the cosine distance matrix between two HMMs \eqn{A} and \eqn{B} before dynamic time warp is applied to
#' the distance matrix calculate the cumulative distance between the HMMs, which acts as a measure of similarity,
#' The cosine distance matrix \eqn{D} is found to be \eqn{D[a_i, b_j] = 1 - \frac{a_ib_j^{T}}{a_ia_i^Tb_jb_j^T}},
#' in which \eqn{a_i} and \eqn{a_i} refer to row vectors of \eqn{A} and \eqn{B} respectively.
#' This in turn means that \eqn{D} is of dimensions \eqn{nrow(A), nrow(b)}. Dynamic time warp then calculates the
#' cumulative distance by calculating matrix \eqn{C[i, j] = min(C[i-1, j], C[i, j-1], C[i-1, j-1]) + D[i, j]},
#' where \eqn{C_{i,j}} is 0 when \eqn{i} or \eqn{j} are less than 1. The lower rightmost point of the matrix \eqn{C}
#' is then returned as the cumulative distance between proteins.
#' @param hmm_1 The name of a profile hidden markov model file.
#' @param hmm_2 The name of another profile hidden markov model file.
#' @return A double that indicates distance between the two proteins.
#' @references Lyons, J., Paliwal, K. K., Dehzangi, A., Heffernan, R., Tsunoda, T., & Sharma, A. (2016).
#' Protein fold recognition using HMM–HMM alignment and dynamic programming.
#' Journal of Theoretical Biology, 393, 67–74.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_distance(system.file("extdata", "1DLHA2-7", package="protHMM"),
#' system.file("extdata", "1TEN-7", package="protHMM"))
hmm_distance<- function(hmm_1, hmm_2){
  text= readLines(hmm_1)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  x[x == "*"]<- 0.0001
  x<- matrix(as.numeric(x), ncol = ncol(x))
  text= readLines(hmm_2)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  y = as.matrix(read.table(text = text[emission])[,3:22])
  y[y == "*"]<- 0.0001
  y<- matrix(as.numeric(y), ncol = ncol(y))
  I<- matrix(0, nrow(x), nrow(y))
  #for(m in 1:nrow(x)){
   # for(n in 1:nrow(y)){
    #  number<- 1 - sum(x[m,] * t(y[n,]))/(sum(x[m,] * t(x[m,])) * sum(y[n,] * t(y[n,])))
    #  I[m,n]<- number
  #  }
  #}
  a_den <- sum(x %*% t(x))
  b_den <- sum(y %*% t(y))
  I<- 1 - (x %*% t(y))/(a_den*b_den)
  D<- I
  for(m in 1:nrow(I)){
    for(n in 1:ncol(I)){
      if(m == 1 && n > 1){
        D[m,n] = I[m,n] + D[m, n-1]
      }
      else if(n == 1 && m > 1){
        D[m,n]= I[m,n] + D[m-1, n]
      }
      else if(n == 1 && m == 1){
        D[m,n] = D[m,n]
      }
      else{
        D[m,n]= I[m,n] + min(c(D[m, n-1],D[m-1, n-1], D[m-1, n]))
      }
    }
  }
  return(D[nrow(D), ncol(D)])
}
