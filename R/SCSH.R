#' @title hmm_SCSH
#' @description This feature returns the 2 and 3-mer compositions of the protein sequence. This is done by first
#' finding all possible 2 and 3-mers for any protein (\eqn{20^2} and \eqn{20^3} permutations for 2 and 3-mers respectively).
#' With those permutations, vectors of length 400 and 8000 are created, each point corresponding to one 2 or 3-mer.
#' Then, the protein sequence that corresponds to the HMM scores is extracted, and put into a bipartite graph with the protein sequence.
#' Each possible path of length 1 or 2 is found, and the corresponding vertices on the graph are noted as 2 and 3-mers.
#' For each 2 or 3-mer found from these paths, 1 is added to the position that responds to that 2/3-mer in the
#' 2-mer and 3-mer vectors , which are the length 400 and 8000 vectors created previously. The vectors are then returned.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 400.
#' @return A vector of length 8000.
#' @references Mohammadi, A. M., Zahiri, J., Mohammadi, S., Khodarahmi, M., & Arab, S. S. (2022).
#' PSSMCOOL: a comprehensive R package for generating evolutionary-based descriptors of protein sequences from PSSM profiles.
#' Biology Methods and Protocols, 7(1).
#' @importFrom utils read.table
#' @importFrom gtools permutations
#' @export
#' @examples
#' h_400<- hmm_SCSH(system.file("extdata", "1DLHA2-7", package="protHMM"))[[1]]
#' h_8000<- hmm_SCSH(system.file("extdata", "1DLHA2-7", package="protHMM"))[[2]]
hmm_SCSH<- function(hmm){ ## gtools
  text= readLines(hmm)
  seq<- unlist(strsplit(text[18], ""))
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
  seq2<- c()
  v<-c("A" ,"C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  for(i in 1:nrow(x)){
    index<- which.max(x[i, ])
    seq2<- c(seq2, v[index])
  }
  k2<- permutations(20, 2, v, repeats.allowed = TRUE)
  k_2<-c()
  for(f in 1:400){
    k_2<-c(k_2, paste0(k2[f, 1],k2[f, 2]))
  }
  k3<- permutations(20, 3, v, repeats.allowed = TRUE)
  k_3<-c()
  for(f in 1:8000){
    k_3<-c(k_3, paste0(k3[f, 1],k3[f, 2], k3[f, 3]))
  }
  out_2<- vector(mode = 'integer', length = 400)
  two_mer<- c()
  three_mer<- c()
  for(i in 1:(length(seq)-1)){
    two_mer<- c(two_mer, paste0(seq[i], seq[i+1]), paste0(seq[i], seq2[i+1]), paste0(seq2[i], seq[i+1]), paste0(seq2[i], seq2[i+1]))
  }
  for(n in 1:length(two_mer)){
    ind<- which(k_2[]==two_mer[n])
    out_2[ind]<- out_2[ind] + 1
  }
  out_3<- vector(length = 8000, mode = 'integer')
  for(i in 1:(length(seq)-2)){
    three_mer<- c(three_mer,
                  paste0(seq[i], seq[i+1], seq[i+2]),
                  paste0(seq2[i], seq[i+1], seq[i+2]),
                  paste0(seq[i], seq2[i+1], seq[i+2]),
                  paste0(seq[i], seq[i+1], seq2[i+2]),
                  paste0(seq2[i], seq2[i+1], seq[i+2]),
                  paste0(seq2[i], seq[i+1], seq2[i+2]),
                  paste0(seq[i], seq2[i+1], seq2[i+2]),
                  paste0(seq2[i], seq2[i+1], seq2[i+2])
    )
  }
  for(m in 1:length(three_mer)){
    ind<- which(k_3[]==three_mer[m])
    out_3[ind]<- out_3[ind] + 1
  }
  return(list(out_2, out_3))
}
