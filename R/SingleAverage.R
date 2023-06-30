#' @title hmm_Single_Average
#' @description This feature groups together rows that are related to the same amino acid. This is done using a vector
#' \eqn{SA(k)}, in which \eqn{k} spans \eqn{1:400} and \eqn{SA(k) = avg_{i = 1, 2... L}H[i, j] \times \delta(P(i), A(z))},
#' in which \eqn{H} is the HMM matrix, \eqn{P} in the protein sequence, \eqn{A} is an ordered set of amino acids,
#'  the variables \eqn{j, z = 1:20}, the variable \eqn{k = j + 20 \times (z-1)} when creating the vector,
#' and \eqn{\delta()} represents Kronecker's delta.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 400.
#' @references Nanni, L., Lumini, A., & Brahnam, S. (2014).
#' An Empirical Study of Different Approaches for Protein Classification. The Scientific World Journal, 2014, 1–17.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_Single_Average(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_Single_Average<- function(hmm){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  seq<- as.matrix(read.table(text = text[emission])[,1])
  x[x == "*"]<- 0
  x[]<- 2^-((0.001)*as.numeric(x))
  x[x == 1]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  v<-c("A" ,"C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  out_v<- vector(length = 400)
  for(m in 1:20){
    index <- which(seq == v[m])
    if(length(index) == 0){
      index = 1:length(seq)
    }
    for(n in 1:20){
      out_v[(m-1)*20+n]<- sum(x[,n][index])
    }
  }
  return(out_v)
}
