## ex: /Users/shayaanemran/Downloads/DD/DDtrain.hhm/1LPBB1-23
## ex2: /Users/shayaanemran/Downloads/DD/DDtrain.hhm/2SNV-14
## goals for functions: folding = 6 (check), structural class = 4, classification = 2, PPI = 2
## progress: SVD, bigrams, trigrams, SD
## check normalization for hmm ac/cc
## protHMM name

hmm_ac<- function(hmm, lg = 4){
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
  out_m<- matrix(0, nrow = lg, ncol = 20)
  for(m in 1:lg){
    for(n in 1:20){
      sum = 0
      for(k in 1:(nrow(x)-lg)){
        sum <- sum + (x[k,n] - mean(x[,n])) * (x[k,n] - mean(x[,n]))
      }
      out_m[m,n]<- sum/(nrow(x)-m)
    }
  }
  return(as.vector(out_m))
} # protein folding 1 (not sure if this works)

hmm_cc<- function(hmm, lg = 4){
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  x = as.matrix(read.table(text = text[which_emmission])[,3:22])
  x[x == "*"]<- 0
  #x[]<- 2^-((0.001)*as.numeric(x))
  #x[x == 1]<- 0
  x<- matrix(as.numeric(x), ncol = ncol(x))
  out_v<- vector(length = 20*19*lg)
  count<- 1
  for(m in 1:20){
    for(n in 1:20){
      sum <- 0
      if(m != n){
        for(i in 1:(nrow(x)-lg)){
          sum = sum + (x[i, m] - mean(x[, m])) * (x[i + lg, n] - mean(x[, n]))
          out_v[count]<- sum/(nrow(x)-lg)
          count = count + 1
          sum = 0
        }
      }
    }
  }
  return(out_v)
} # protein folding 2 (unfinished)

hmm_bigrams<- function(hmm){ ##HMM bigrams from Xia et al. (2017)
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  df_emmission = as.matrix(read.table(text = text[which_emmission])[,3:22])
  df_emmission[df_emmission == "*"]<- 0
  df_emmission[]<- 2^-((0.001)*as.numeric(df_emmission))
  df_emmission[df_emmission == 1]<- 0
  out_df<- matrix(nrow = 20, ncol = 20)
  for(n in 1:20){
    for(m in 1:20){
      sum<-0
      for(k in 1:(nrow(df_emmission)-1)){
        sum<- sum + as.numeric(df_emmission[k, m])*as.numeric(df_emmission[k+1, n])
      }
      out_df[m,n]<- sum
    }
  }
  out_v<- as.vector(out_df)
  return(out_v)
} # protein folding 3

hmm_SD<- function(hmm){ ##Separated Dimers for HMMs
  text= readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  df_emmission = as.matrix(read.table(text = text[which_emmission])[,3:22])
  df_emmission[df_emmission == "*"]<- 0
  df_emmission[]<- 2^-((0.001)*as.numeric(df_emmission))
  df_emmission[df_emmission == 1]<- 0
  out_df<- matrix(nrow = 20, ncol = 20)
  for(n in 1:20){
    for(m in 1:20){
      sum<-0
      for(k in 1:(nrow(df_emmission)-7)){
        sum<- sum + as.numeric(df_emmission[k, m])*as.numeric(df_emmission[k+7, n])
      }
      out_df[m,n]<- sum
    }
  }
  out_v<- as.vector(out_df)
  return(out_v)
} # protein folding 4

hmm_trigrams<- function(hmm){ ##HMM trigrams from Lyons et al. (dev)
  text = readLines(hmm)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  df_emmission = as.matrix(read.table(text = text[which_emmission])[,3:22])
  df_emmission[df_emmission == "*"]<- 0
  df_emmission[]<- 2^-((0.001)*as.numeric(df_emmission))
  df_emmission[df_emmission == 1]<- 0
  out_v<- vector(length = 8000, mode = 'integer')
  count = 0
  for(n in 1:20){
    for(m in 1:20){
      for(r in 1:20){
        sum<-0
        for(k in 1:(nrow(df_emmission)-2)){
          sum<- sum + as.numeric(df_emmission[k, m])*as.numeric(df_emmission[k+1, n])*as.numeric(df_emmission[k+2, r])
        }
        out_v[count]<- sum
        count = count+1
      }
    }
  }
  return(out_v)
} # protein folding 5

hmm_align<- function(hmm_1, hmm_2){ ## this should work... hopefully
  ##read hmms
  text= readLines(hmm_1)
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
  text= readLines(hmm_2)
  start = grep("HMM", (text))
  start = start[length(start)]
  end = grep("//", text)
  text = text[start:end]
  which_emmission = grep(" [0-9]{1,9} ", text)
  y = as.matrix(read.table(text = text[which_emmission])[,3:22])
  y[y == "*"]<- 0
  y[]<- 2^-((0.001)*as.numeric(y))
  y[y == 1]<- 0
  y<- matrix(as.numeric(y), ncol = ncol(y))
  I<- matrix(0, nrow(x), nrow(y))
  ## dissimilarity matrix
  for(m in 1:nrow(x)){
    for(n in 1:nrow(y)){
      number<- 1 - sum(x[m,] * t(y[n,]))/sqrt(sum(x[m,] * t(x[m,])) * sum(y[n,] * t(y[n,])))
      I[m,n]<- number
    }
  }
  D<- I
  for(m in 1:nrow(I)){
    for(n in 1:ncol(I)){
      if(m == 1 && n > 1){
        D[m,n]= I[m,n] + D[m, n-1]
      }
      else if(n == 1){
        D[m,n]= I[m,n] + D[m-1, n]
      }
      else{
        D[m,n]= I[m,n] + min(c(D[m, n-1],D[m-1, n-1], D[m-1, n]))
      }
    }
  }
  return(D[nrow(D), ncol(D)])
} # protein folding 6 (unfinished)

hmm_SVD<- function(hmm){
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
  return(svd(x)[[1]])
} # classification 1

hmm_Single_Average<- function(hmm){
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
  return(svd(x)[[1]])
} # classification 2

hmm_discrete_cosine<- function(hmm){
  
} # PPI 1

hmm_SCSH2<- function(hmm){
  
} # PPI 2

hmm_read<- function(hmm){
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
  return(x)
} # generic read function