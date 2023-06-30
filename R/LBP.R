#' @title hmm_LBP
#' @description This feature uses local binary pattern with a neighborhood of radius 1 and 8 sample points to extract
#' features from the HMM. A 256 bin histogram is extracted as a 256 length feature vector.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 256.
#' @references Li, Y., Li, L., Wang, L., Yu, C., Wang, Z., & You, Z. (2019).
#' An Ensemble Classifier to Predict Protein–Protein Interactions by Combining PSSM-based Evolutionary Information with Local Binary Pattern Model.
#' International Journal of Molecular Sciences, 20(14), 3511.
#' @importFrom utils read.table
#' @export
#' @examples
#' h<- hmm_LBP(system.file("extdata", "1DLHA2-7", package="protHMM"))
#'
hmm_LBP<- function(hmm){
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
  col_num<- ncol(x)
  row_num<- nrow(x)
  x<- matrix(as.numeric(x), ncol = ncol(x))
  r<- 1
  points<- list(c(1, 0),
                c(1, 1),
                c(0, 1),
                c(-1, 1),
                c(-1, 0),
                c(-1, -1),
                c(0, -1),
                c(1, -1))
  bin_points<- list()
  centers<-x[2:(row_num-1),2:(col_num-1)]
  for(i in 1:length(points)){
    dc<- points[[i]][1]
    dr<- points[[i]][2]
    bin_points[[i]]<- x[(2+dr):(row_num-1+dr), (2+dc):(col_num-1+dc)] - centers
  }
  for(i in 1:length(bin_points)){
    bin_points[[i]]<- ifelse(bin_points[[i]]>=0, 1, 0)
  }
  pixel_matrix<- matrix(0, nrow = row_num - 2, ncol = col_num - 2)
  for(j in 1:ncol(pixel_matrix)){
    for(i in 1:nrow(pixel_matrix)){
      binary_num<- c()
      for(z in 1:length(bin_points)){
        binary_num<- c(binary_num, bin_points[[z]][i, j])
      }
      pixel_matrix[i, j]<- sum(2^(which(binary_num == 1)-1))
    }
  }
  out_v<- vector('integer', 256)
  for(i in 0:255){
    out_v[i+1]<- length(which(pixel_matrix == i))
  }
  return(out_v)
}
