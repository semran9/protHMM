#' @title chmm
#' @description This feature begins by creating a CHMM, which is created by constructing 4 matrices, \eqn{A, B, C, D} from
#' the original HMM \eqn{H}. \eqn{A} contains the first 75 percent of the original matrix \eqn{H} row-wise, \eqn{B} the
#' last 75 percent, \eqn{C} the middle 75 percent and \eqn{D} the entire original matrix. These are then merged to create the new
#' CHMM \eqn{Z}. From there, the Bigrams feature is calculated with a flattened 20 x 20 matrix \eqn{B}, in which \eqn{B[i, j] = \sum_{a = 1}^{L-1} Z_{a, i} \times Z_{a+1, j}}.
#' \eqn{H} corresponds to the original HMM matrix, and \eqn{L} is the number of rows in \eqn{Z}. Local Average Group,
#' or LAG is then calculated by splitting up the CHMM into 20 groups along the length of the protein sequence and calculating
#' the sums of each of the columns of each group, making a 1 x 20 vector per group, and a length 20 x 20 vector for all groups. These features are then fused.
#' @param hmm The name of a profile hidden markov model file.
#' @return A fusion vector of length 800.
#' @return A LAG vector of length 400.
#' @return A Bigrams vector of length 400.
#' @references An, J., Zhou, Y., Zhao, Y., & Yan, Z. (2019).
#' An Efficient Feature Extraction Technique Based on Local Coding PSSM and Multifeatures Fusion for Predicting Protein-Protein Interactions.
#' Evolutionary Bioinformatics, 15, 117693431987992.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- chmm(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
chmm<- function(hmm){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  x[x == "*"]<- 0
  x[]<- 2^-((0.001)*as.numeric(x))
  x[x == 1]<- 10^(-10)
  x<- matrix(as.numeric(x), ncol = ncol(x))
  l = nrow(x)
  l_1 = round(l*.25, digits = 0)
  l_2 = round(l*.75, digits = 0)
  A<- x[1:l_2,]
  B<- x[l_2:l,]
  C<- x[l_1:l_2,]
  D<- x
  hmm_mat<- rbind(A, B, C, D)
  ##lag
  n = nrow(hmm_mat)
  groups<- vector('integer', 21)
  groups[1]<- 1
  lag_v<- vector('integer', 400)
  for(i in 2:21){

    groups[i]<- round(0.05*(i-1)*(n), 0)
  }
  count = 1
  for(i in 1:20){
    group_m<- hmm_mat[groups[i]:groups[i+1],]
    for(j in 1:20){
      sum<- sum(group_m[, j])
      sum<- sum*(20/n)
      lag_v[count]<- sum
      count = count +1
    }
  }
  ##bigrams
  bi_m<- matrix(0, 20, 20)
  for(i in 1:20){
    for(j in 1:20){
      sum_bi = 0
      for(k in 1:(n-1)){
        sum_bi = sum_bi + hmm_mat[k, i]*hmm_mat[k+1, j]
      }
      bi_m[i,j]<- sum_bi
    }
  }
  out_v<- c(lag_v, as.vector(bi_m))
  return(list(out_v, lag_v, as.vector(bi_m)))
}
