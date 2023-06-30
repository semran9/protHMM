#' @title hmm_read
#' @description Reads in the amino acid emission frequency columns of a profile hidden markov model matrix
#' and converts each position to frequencies.
#' @param hmm The name of a profile hidden markov model file.
#' @return A 20 x L matrix, in which L is the sequence length.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_read(system.file("extdata", "1DLHA2-7", package="protHMM"))

hmm_read<- function(hmm){
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
  return(x)
}
