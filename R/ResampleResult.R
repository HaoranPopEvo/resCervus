#ResampleResult---------------------------
#' Summarize Parentage Result
#'
#' \code{ResampleResult} reports parentage result from summarizing parent-offspring
#' pairs that occur in each replicate.
#' @param filepath path of directory that contains parentage result files. It is the
#' same directory assigned by parameter \code{Data_directory} in function "\code{ResampleCervus}".
#' @param Add_sex whether sex data is used in parentage analysis.
#' @return a \code{data.frame} of four columns. The first two columns are parents. The
#' third column shows information about offspring. The last column is per–triad support value,
#' which is defined as the number of replicates that support a triad divided by the
#' number of samples that include the individuals of the triad.
#' @importFrom utils read.csv
#' @importFrom stats na.omit
#' @examples
#' \dontrun{
#' Data_directory <- "D:/xxx"
#' result <- ResampleResult(Data_directory, Add_sex = TRUE)
#' print(result)
#' }
#' @author Hao-Ran Wu
#' @export
ResampleResult <- function(filepath, Add_sex = FALSE){
  current_path <- getwd()

  setwd(filepath)
  files <- list.files(filepath)
  files.res <- files[grepl("Res.csv",files)]
  if(Add_sex){
    files.dad <- files[grepl("Dad.csv",files)]
    files.mom <- files[grepl("Mom.csv",files)]
  } else{
    files.par <- files[grepl("Par.csv",files)]
  }
  files.off <- files[grepl("Off.csv",files)]

  trioALL <- do.call(rbind
                     ,sapply(strsplit(files.res," Res.csv"),function(xx){
                       dat.res <- read.csv(paste(xx,"Res.csv"),row.names = NULL)
                       if(all(is.na(dat.res[,3:ncol(dat.res)]))){ #nothing assigned
                         return(NULL)
                       }

                       if(all(is.na(dat.res[,ncol(dat.res)]))){ #column misplacement
                         colNames <- colnames(dat.res)[-1]
                         dat.res <- dat.res[,-ncol(dat.res)]
                         colnames(dat.res) <- colNames
                       }

                       if(is.numeric(na.omit(dat.res[,ncol(dat.res)])) || all(is.na(dat.res[,ncol(dat.res)]))){
                         return(NULL) #nothing assigned
                       }

                       is.signif <- (dat.res[,ncol(dat.res)]=="*")
                       if(sum(is.signif)==0){
                         NULL #nothing assigned
                       } else{
                         dat.res <- dat.res[is.signif,]
                         if(Add_sex){
                           data.frame(
                             P1=dat.res$Candidate.father.ID,
                             P2=dat.res$Candidate.mother.ID,
                             Off=dat.res$Offspring.ID
                           )
                         } else{
                           data.frame(
                             P1=dat.res$First.candidate.ID,
                             P2=dat.res$Second.candidate.ID,
                             Off=dat.res$Offspring.ID
                           )
                         }
                       }
                     },simplify = FALSE))
  colnames(trioALL)<-c("P1.ID","P2.ID","Off.ID")
  family_table <- table(paste(trioALL$P1.ID, ":", trioALL$P2.ID,
                              ":", trioALL$Off.ID, sep = ""))
  IndList <- lapply(1:length(files.off),function(xx){
    if(Add_sex){
      c(read.csv(files.dad[xx])[,1],read.csv(files.mom[xx])[,1],read.csv(files.off[xx])[,1])
    } else{
      c(read.csv(files.par[xx])[,1],read.csv(files.off[xx])[,1])
    }
  })
  famComplete_count <- c()
  for (i in 1:length(family_table)) {
    fnID <- strsplit(names(family_table[i]), ":")[[1]]
    fncount <- sum(sapply(IndList, function(x) all(fnID %in%
                                                     x)))
    famComplete_count <- c(famComplete_count, fncount)
  }
  family_freq <- family_table/famComplete_count
  triodata <- data.frame(do.call(rbind, strsplit(names(family_freq), ":")), support = as.numeric(family_freq))
  triodata <- triodata[order(triodata$support, decreasing = TRUE),]
  dat <- triodata

  #combine same families
  large <- ifelse(dat[,1]>=dat[,2],dat[,1],dat[,2])
  small <- ifelse(dat[,1]>=dat[,2],dat[,2],dat[,1])
  dat[,1] <- large
  dat[,2] <- small

  dup <- duplicated(paste(dat[,1],dat[,2],dat[,3]))
  while(any(dup)){
    index <- which(dup)[1]
    index.other <- which(dat[index,1]==dat[,1]&dat[index,2]==dat[,2]&dat[index,3]==dat[,3])
    resample.new <- sum(dat[index.other,]$support)
    new.record <- dat[index,]
    new.record$support <- resample.new

    dat <- rbind(dat[-index.other,],new.record)
    dup <- duplicated(paste(dat[,1],dat[,2],dat[,3]))
  }

  dat[order(dat$support,decreasing = TRUE),]
}
