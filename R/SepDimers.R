#' @title hmm_SepDim
#'
#' @description This feature calculates the probabilistic expression of amino acid dimers that are spatially separated by a distance \eqn{l}.
#' Mathematically, this is done with a 20 x 20 matrix \eqn{F}, in which \eqn{F[m, n] = \sum_{i = 1}^{L-l} H_{i, m}H_{i+k, n}}.
#' \eqn{H} corresponds to the original HMM matrix, and \eqn{L} is the number of rows in \eqn{H}. Matrix \eqn{F} is then flattened to
#' a feature vector of length 400, and returned.
#'
#' @param hmm The name of a profile hidden markov model file.
#' @param l Spatial distance between dimer residues.
#' @return A vector of length 400
#' @references Saini, H., Raicar, G., Sharma, A., Lal, S. K., Dehzangi, A., Lyons, J., Paliwal, K. K., Imoto, S., & Miyano, S. (2015).
#' Probabilistic expression of spatially varied amino acid dimers into general form of Chou's pseudo amino acid composition for protein fold recognition.
#' Journal of Theoretical Biology, 380, 291–298.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_SepDim(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'

hmm_SepDim<- function(hmm, l = 7){ ##Separated Dimers for HMMs
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
      for(k in 1:(nrow(x)-l)){
        sum<- sum + x[k, m]*x[k+l, n]
      }
      out_df[m,n]<- sum
    }
  }
  out_v<- as.vector(out_df)
  return(out_v)
}
