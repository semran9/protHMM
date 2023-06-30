#' @title hmm_GA
#' @description This feature calculates the Geary autocorrelation of each amino acid type for each distance
#' d less than or equal to the lag value and greater than or equal to 1.
#' @param hmm The name of a profile hidden markov model file.
#' @param lg The lag value, which indicates the distance between residues.
#' @note
#' The lag value must be less than the length of the protein sequence
#' @return A vector of length lg \eqn{\times} 20, by default this is 180.
#' @references Liang, Y., Liu, S., & Zhang. (2015).
#' Prediction of Protein Structural Class Based on Different Autocorrelation Descriptors of Position–Specific Scoring Matrix.
#' MATCH: Communications in Mathematical and in Computer Chemistry, 73(3), 765–784.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_GA(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_GA<- function(hmm, lg = 9){
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
  out_m<- matrix(0, lg, 20)
  for(d in 1:lg){
    for(j in 1:20){
      s<- 0
      for(i in 1:(nrow(x)-d)){
        s<- s + (x[i, j] - x[i+d, j])^2
      }
      n<- 0
      for(i in 1:nrow(x)){
        n<- n+(x[i,j] - mean(x[,j]))^2
      }
      s<-s/(2*(nrow(x)-d))
      n<-n/(nrow(x)-1)
      out_m[d,j]<- s/n
    }
  }
  out_m[is.na(out_m)]<-0
  return(as.vector(out_m))
}
