ReadCoGAPSResults <- function(path=getwd(), output.list=TRUE) {
  origDir <- getwd()

  message(origDir)
  message('\n')
  
  setwd(path)

  # identify experiments in the folder
  PMeanFiles <- list.files(pattern='Pmean',full.names=F)
  fileIDS <- sapply(strsplit(PMeanFiles,'\\.'),function(x){x[3]})
  names(PMeanFiles) <- fileIDS
  
  # find other files
  PSDFiles <- paste('Psd','0',fileIDS,'txt',sep=".")
  names(PSDFiles) <- fileIDS

  AMeanFiles <- paste('Amean','0',fileIDS,'txt',sep=".")
  names(AMeanFiles) <- fileIDS
  
  ASDFiles <- paste('Asd','0',fileIDS,'txt',sep=".")
  names(ASDFiles) <- fileIDS
  
  # check if files are in the folder
  files <-   c(PSDFiles,AMeanFiles,ASDFiles)
  if (!all(file.exists(files))) {
    missingFiles <- files[which(!file.exists(files))]
    stop(paste('Cannot read CoGAPS results: missing files:',
         paste(missingFiles,collapse=",")))
  }

  # read in data from the files
  A.mean <- list()
  P.mean <- list()
  M <- list()
  A.sd <- list()
  P.sd <- list()
  for (ID in fileIDS) {
    A.mean[[ID]] <- as.matrix(read.table(AMeanFiles[ID], 
                              header=T,row.names=1,sep="\t"))
    A.sd[[ID]] <- as.matrix(read.table(ASDFiles[ID], 
                            header=T,row.names=1,sep="\t"))

    P.mean[[ID]] <- as.matrix(read.table(PMeanFiles[ID], 
                              header=T,row.names=1,sep="\t"))
    P.sd[[ID]] <- as.matrix(read.table(PSDFiles[ID], 
                            header=T,row.names=1,sep="\t"))

    M[[ID]] <- A.mean[[ID]]%*%P.mean[[ID]]

  }
  
  # return to original directory
  setwd(origDir)

  #return files
  if (output.list) { 
    results <- list(A.mean=A.mean, A.sd=A.sd, P.mean=P.mean, P.sd=P.sd, M=M)
  } else {
    A.mean.matrix <- A.mean[[fileIDS[1]]]
    colnames(A.mean.matrix) <- paste(colnames(A.mean[[fileIDS[1]]]),
       rep(fileIDS[1],ncol(A.mean[[fileIDS[1]]])), sep=".")
    A.sd.matrix <- A.sd[[fileIDS[1]]]
    colnames(A.sd.matrix) <- paste(colnames(A.sd[[fileIDS[1]]]),
       rep(fileIDS[1],ncol(A.sd[[fileIDS[1]]])), sep=".")
    P.mean.matrix <- P.mean[[fileIDS[1]]]
    row.names(P.mean.matrix) <- paste(row.names(P.mean[[fileIDS[1]]]),
       rep(fileIDS[1],nrow(P.mean[[fileIDS[1]]])), sep=".")
    P.sd.matrix <- P.sd[[fileIDS[1]]]
    row.names(P.sd.matrix) <- paste(row.names(P.sd[[fileIDS[1]]]),
       rep(fileIDS[1],nrow(P.sd[[fileIDS[1]]])), sep=".")
    M.matrix <- M[[fileIDS[1]]]
    colnames(M.matrix) <- paste(colnames(M[[fileIDS[1]]]),
       rep(fileIDS[1],ncol(M[[fileIDS[1]]])), sep=".")  
    }
    
    if (length(fileIDS)>1) {
      for (ID in fileIDS[2:length(fileIDS)]) {

        A.mean.matrix <- cbind(A.mean.matrix, A.mean[[ID]])
        colnames(A.mean.matrix)[(ncol(A.mean.matrix)-ncol(A.mean[[ID]])+1):ncol(A.mean.matrix)] <- paste(colnames(A.mean[[ID]]),rep(ID,ncol(A.mean[[ID]])), sep=".")
        A.sd.matrix <- cbind(A.sd.matrix, A.sd[[ID]])
        colnames(A.sd.matrix)[(ncol(A.sd.matrix)-ncol(A.sd[[ID]])+1):ncol(A.sd.matrix)] <- paste(colnames(A.sd[[ID]]),rep(ID,ncol(A.sd[[ID]])), sep=".")
        
        M.matrix <- cbind(M.matrix, M[[ID]])
        colnames(M.matrix)[(ncol(M.matrix)-ncol(M[[ID]])+1):ncol(M.matrix)] <- paste(colnames(M[[ID]]),rep(ID,ncol(M[[ID]])), sep=".")
        
        P.mean.matrix <- rbind(P.mean.matrix, P.mean[[ID]])
        row.names(P.mean.matrix)[(nrow(P.mean.matrix)-nrow(P.mean[[ID]])+1):nrow(P.mean.matrix)] <- paste(row.names(P.mean[[ID]]),rep(ID,nrow(P.mean[[ID]])), sep=".")
        P.sd.matrix <- rbind(P.sd.matrix, P.sd[[ID]])
        row.names(P.sd.matrix)[(nrow(P.sd.matrix)-nrow(P.sd[[ID]])+1):nrow(P.sd.matrix)] <- paste(row.names(P.sd[[ID]]),rep(ID,nrow(P.sd[[ID]])), sep=".")
        
      }
    }

    results <- list(A.mean=A.mean.matrix, A.sd=A.sd.matrix,
                    P.mean=P.mean.matrix, P.sd=P.sd.matrix, M=M.matrix)
  }

  return(results)
  
}
