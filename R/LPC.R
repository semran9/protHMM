#' @title hmm_LPC
#' @description This feature uses linear predictive coding (LPC) to map each HMM to a \eqn{20 \times 14 = 280} dimensional vector,
#' where for each of the 20 columns of the HMM, LPC is used to extract a 14 dimensional vector \eqn{D_n}
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 280.
#' @references Qin, Y., Zheng, X., Wang, J., Chen, M., & Zhou, C. (2015).
#' Prediction of protein structural class based on Linear Predictive Coding of PSI-BLAST profiles.
#' Central European Journal of Biology, 10(1).
#' @importFrom utils read.table
#' @importFrom phonTools lpc
#' @export
#' @examples
#' h<- hmm_LPC(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_LPC<- function(hmm){
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
  out_m<- matrix(0, 14, 20)
  for(j in 1:20){
    out_m[, j]<- phonTools::lpc(x[, j])
  }
  return(as.vector(out_m))
}
