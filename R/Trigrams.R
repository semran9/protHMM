#' @title hmm_trigrams
#' @description This feature is calculated with a 20 x 20 x 20 block \eqn{B}, in which \eqn{B[i, j, k] = \sum_{a = 1}^{L-2} H_{a, i}H_{a+1, j}H_{a+2, k}}.
#' \eqn{H} corresponds to the original HMM matrix, and \eqn{L} is the number of rows in \eqn{H}. Matrix \eqn{B} is then flattened to
#' a feature vector of length 8000, and returned.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 8000
#' @references Lyons, J., Dehzangi, A., Heffernan, R., Yang, Y., Zhou, Y., Sharma, A., & Paliwal, K. K. (2015).
#' Advancing the Accuracy of Protein Fold Recognition by Utilizing Profiles From Hidden Markov Models.
#' IEEE Transactions on Nanobioscience, 14(7), 761–772.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_trigrams(system.file("extdata", "1DLHA2-7", package="protHMM"))

hmm_trigrams<- function(hmm){ ##HMM trigrams from Lyons et al. (dev)
  text = readLines(hmm)
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
  out_v<- vector(length = 8000, mode = 'integer')
  count = 0
  for(n in 1:20){
    for(m in 1:20){
      for(r in 1:20){
        sum<-0
        for(k in 1:(nrow(x)-2)){
          sum<- sum + x[k, m] * x[k+1, n] * x[k+2, r]
        }
        out_v[count]<- sum
        count = count+1
      }
    }
  }
  return(out_v)
}
