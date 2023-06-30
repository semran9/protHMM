#' @title hmm_smooth
#' @description This feature smooths the HMM matrix \eqn{H} by using sliding window of length \eqn{sw} to incorporate information
#' from up and downstream residues into each row of the HMM matrix. Each HMM row \eqn{r_i} is made into the summation
#' of \eqn{r_{i-(sw/2)}+... r_i...+r_{i+(sw/2)}}, for \eqn{i = 1:L}, where \eqn{L} is the number of rows in \eqn{H}.
#' For rows such as the beginning and ending rows, \eqn{0} matrices of dimensions \eqn{sw/2, 20} are appended to the
#' original matrix \eqn{H}.
#' @param hmm The name of a profile hidden markov model file.
#' @param sw The size of the sliding window.
#' @return A matrix of dimensions L \eqn{\times} 20.
#' @references Fang, C., Noguchi, T., & Yamana, H. (2013).
#' SCPSSMpred: A General Sequence-based Method for Ligand-binding Site Prediction.
#' IPSJ Transactions on Bioinformatics, 6(0), 35–42.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_smooth(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_smooth<- function(hmm, sw = 7){
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
  L<- nrow(x)
  x<- matrix(as.numeric(x), ncol = ncol(x))
  smooth<- x
  zero_m<- matrix(0, (sw-1)/2, 20)
  start<- (sw-1)/2
  x<- rbind(zero_m, x, zero_m)
  for(i in start:(L+start)){
    sum<-x[i]
    for(j in 1:start){
      sum<- rbind(x[i+j, ], sum,x[i-j, ])
    }
    smooth[i-start,]<- colSums(sum)
  }

  return (smooth)
}
