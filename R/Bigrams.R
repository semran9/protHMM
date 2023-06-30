#' @title hmm_bigrams
#' @description This feature is calculated with a 20 x 20 matrix \eqn{B}, in which \eqn{B[i, j] = \sum_{a = 1}^{L-1} H_{a, i}H_{a+1, j}}.
#' \eqn{H} corresponds to the original HMM matrix, and \eqn{L} is the number of rows in \eqn{H}. Matrix \eqn{B} is then flattened to
#' a feature vector of length 400, and returned.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 400
#' @references Lyons, J., Dehzangi, A., Heffernan, R., Yang, Y., Zhou, Y., Sharma, A., & Paliwal, K. K. (2015).
#' Advancing the Accuracy of Protein Fold Recognition by Utilizing Profiles From Hidden Markov Models.
#' IEEE Transactions on Nanobioscience, 14(7), 761–772.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_bigrams(system.file("extdata", "1DLHA2-7", package="protHMM"))
hmm_bigrams<- function(hmm){
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
  out_df<- matrix(nrow = 20, ncol = 20)
  for(n in 1:20){
    for(m in 1:20){
      sum<-0
      for(k in 1:(nrow(x)-1)){
        sum<- sum + x[k, m]*x[k+1, n]
      }
      out_df[m,n]<- sum
    }
  }
  out_v<- as.vector(out_df)
  return(out_v)
}
