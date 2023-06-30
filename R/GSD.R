#' @title hmm_GSD
#' @description This feature initially creates a grouping matrix \eqn{G} by assigning each position a number \eqn{1:3} based on
#' the value at each position of HMM matrix \eqn{H}; \eqn{1} represents the low probability group, \eqn{2} the medium and \eqn{3} the high probability group.
#' The number of total points in each group for each column is then calculated, and the sequence is then split
#' based upon the the positions of the 1st, 25th, 50th, 75th and 100th percentile (last) points for each of the three groups,
#' in each of the 20 columns of the grouping matrix. Thus for column \eqn{j}, \eqn{S(k, j, z) = \sum_{i = 1}^{(z)*.25*N} |G[i, j] = k|},
#' where \eqn{k} is the group number, \eqn{z = 1:4} and \eqn{N} corresponds to number of rows in matrix \eqn{G}.
#' @param hmm The name of a profile hidden markov model file.
#' @return A vector of length 300.
#' @references Jin, D., & Zhu, P. (2021).
#' Protein Subcellular Localization Based on Evolutionary Information and Segmented Distribution.
#' Mathematical Problems in Engineering, 2021, 1–14.
#' @importFrom utils read.table
#' @importFrom stats quantile
#' @importFrom stats sd
#' @export
#' @examples
#' h<- hmm_GSD(system.file("extdata", "1DLHA2-7", package="protHMM"))

hmm_GSD<- function(hmm){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  emission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[emission])[,3:22])
  seq<- as.matrix(read.table(text = text[emission])[,1])
  x[x == "*"]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  x[]<- 2^-((0.001)*as.numeric(x))
  x[x == 1]<- 0
  g<- x
  out_v<- c()
  for(j in 1:20){
    for(i in 1:nrow(x)){
      if(x[i,j] >= 0 && x[i,j] <= quantile(x, 0.33)){
        g[i,j]<-  1
      }
      else if(x[i,j] > quantile(x, 0.33) && x[i,j] <= quantile(x, 0.66)){
        g[i,j]<-  2
      }
      else{
        g[i, j]<- 3
      }
    }
  }
  t<-matrix(0, 3, 20)
  for(j in 1:20){
    t[1, j]<- sum(g[, j][g[, j] == 1])
    t[2, j]<- sum(g[, j][g[, j] == 2])/2
    t[3, j]<- sum(g[, j][g[, j] == 3])/3
  }
  out_v<- vector(300, mode = 'integer')
  count<- 1
  for(j in 1:20){
    for(i in 1:3){
      if(t[i, j] > 0){
        indc<- which(g[, j] == i)
        out_v[count]<- indc[1]
        count<- count +1
        sum<-g[indc[1], j]
        c<-1
        while(sum < 0.25*t[i, j]){
          sum<- sum + g[indc[c], j]
          c<- c+1
        }
        out_v[count]<- indc[c]
        count<- count+1
        while(sum < 0.5*t[i, j]){
          sum<- sum + g[indc[c], j]
          c<- c+1
        }
        out_v[count]<- indc[c]
        count<- count+1
        while(sum < 0.75*t[i, j]){
          sum<- sum + g[indc[c], j]
          c<- c+1
        }
        out_v[count]<- indc[c]
        count<- count+1
        c<- length(indc)
        out_v[count]<- indc[c]
        count<- count+1
      }
      else{
        k<- 1
        while(k < 6){
          out_v[count]<- 0
          count<- count+1
          k<- k + 1
        }
      }
    }
  }
  s_m<- mean(out_v)
  s_sd<- sd(out_v)
  out_v<- (s_m-out_v)/s_sd
  return(out_v)
}
