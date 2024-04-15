#ResampleCervus---------------------------
#' Parentage from Resampling
#'
#' \code{ResampleCervus} is used to perform parentage analysis by CervusCL using
#' resampling approach.
#' @param CervusCL_directory path of directory that contains "\code{CervusCL.exe}"
#' @param Data_directory path of directory that contains data files, which include
#' genotypic data (as well as individual ID), parent ID, offspring ID, and gender
#' information (optional). All files should be presented as CSV format.
#' @param Genotype_file name of genotypic data file. Genotype information should be
#' organized as a format that can be read by CervusCL. Each row is an individual (except
#' the header row), and every two columns is a loci (except the ID column). The ID
#' column should be previous to loci information.
#' @param Genotype_file.ID_in_column index of ID column for genotypic data file.
#' @param Genotype_file.First_allele_in_column index of the start column of loci information
#' for genotypic data file.
#' @param Parent_file name of parent ID file. Only the first column is used unless
#' \code{Parent_file.ID_in_column} is modified. It contains all IDs for parental
#' individuals.
#' @param Offspring_file name of offspring ID file. Only the first column is used
#' unless \code{Offspring_file.ID_in_column} is modified. It contains all IDs for
#' offspring individuals.
#' @param Simulate.offspring number of simulated offspring. Defaulted is 10000, as
#' consistent with Cervus. Change of this value is not suggested.
#' @param Simulate.Prop.sampled proportion of sampled individuals in the population.
#' Default is 0.8 but its actual value may vary among different populations although
#' estimation of this value is challenging.
#' @param Simulate.Prop.mistyped genotyping error rate which may due to alleleic dropout,
#' stochastic typing error, mutation, etc. Estimation of this value is also challenging.
#' @param Simulate.Minimum_typed_loci minimum typed loci. Individuals with typed loci
#' lower than this value will be discarded. Defaulted is the half of available loci, as
#' consistent with Cervus.
#' @param Genotype_file.HeaderRow whether to contain header row in genotypic data file.
#' @param Parent_file.HeaderRow whether to contain header to in parental ID file.
#' @param Parent_file.ID_in_column index of ID column in parental ID file.
#' @param Offspring_file.HeaderRow whether to contain header to in offspring ID file.
#' @param Offspring_file.ID_in_column index of ID column in offspring ID file.
#' @param Add_sex whether to add gender information.
#' @param Add_sex.Gender_file name of gender data file. The file should contain only
#' one column with each row showing gender information for an individual. The order of
#' row should be consistent with that in genotypic file. Gender can be donated by any
#' numbers and/or characters. For example, use 1 to donate males and 2 donating females.
#' Note that \code{1 < 2}, and the smaller level is considered as males while higher level
#' for females.
#' @param Add_sex.Simulate.Prop.father.sampled when using gender information,
#' \code{Simulate.Prop.sampled} is used to donate sampling rate for mothers while
#' this parameter is for fathers.
#' @param resample.prop proportion of resampled individuals for each replicate.
#' @param resample.nsim number of replicates.
#' @param ignore.stdout inhibit outputs when executing CervusCL.
#' @param ignore.hasexecuted ignore files that has been executed. Sometimes, you may
#' finish CervusCL algorithm for part of the replicates while the remaining is not
#' executed yet. In the next time, you can continue your work by setting this parameter to
#' \code{TRUE} to perform parentage for the remaining files.
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @return no returned values.
#' @examples
#' \dontrun{
#' CervusCL_directory <- "D:/xxx/Cervus/Cervus CL"
#' Data_directory <- "D:/xxx"
#' Genotype_file <- "Genotype_data.csv"
#' Genotype_file.ID_in_column <- 1
#' Genotype_file.First_allele_in_column <- 2
#' Parent_file <- "Parent_ID_data.csv"
#' Offspring_file <- "Offspring_ID_data.csv"
#' Gender_file <- "Gender_data.csv"
#'
#' ResampleCervus(CervusCL_directory, Data_directory,
#'                Genotype_file, Genotype_file.ID_in_column, Genotype_file.First_allele_in_column,
#'                Parent_file, Offspring_file, ignore.hasexecuted = TRUE,
#'                Add_sex = TRUE, Add_sex.Gender_file = Gender_file, resample.prop = 0.8)
#' }
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998). Statistical confidence for likelihood-based paternity inference in natural populations. Molecular Ecology, 7, 639-655.
#' @author
#' Hao-Ran Wu
#' @export
#'
ResampleCervus <- function(CervusCL_directory, Data_directory,

                           Genotype_file, Genotype_file.ID_in_column, Genotype_file.First_allele_in_column,
                           Parent_file, Offspring_file,

                           Simulate.offspring = 10000, Simulate.Prop.sampled = 0.8, Simulate.Prop.mistyped = 0.001,
                           Simulate.Minimum_typed_loci = "defaulted",

                           Genotype_file.HeaderRow = 1,
                           Parent_file.HeaderRow = 1,
                           Parent_file.ID_in_column = 1,
                           Offspring_file.HeaderRow = 1,
                           Offspring_file.ID_in_column = 1,

                           Add_sex = FALSE,
                           Add_sex.Gender_file = NULL,
                           Add_sex.Simulate.Prop.father.sampled = 0.8,

                           resample.prop = 0.8,
                           resample.nsim = 100,
                           ignore.stdout = TRUE,
                           ignore.hasexecuted = FALSE){

  if(Add_sex && is.null(Add_sex.Gender_file)) stop("no gender file specified.")

  current_path <- getwd() #get current path
  if(grepl("/",CervusCL_directory)){
    CervusCL_directory <- paste(strsplit(CervusCL_directory,"/")[[1]],collapse="\\")
  } #change format of filepath if "/" occurs
  if(grepl("/",Data_directory)){
    Data_directory <- paste(strsplit(Data_directory,"/")[[1]],collapse="\\")
  } #change format of filepath if "/" occurs

  #generate `.crv` and related files
  setwd(Data_directory)
  if(ignore.hasexecuted==TRUE){
    #ignore creating files
    local_files <- list.files()
    local_files <- local_files[grepl(".crv",local_files)]
    max_resample_value_in_local_files <- max(as.numeric(sapply(strsplit(local_files," "),function(xx) xx[2])))
    if(max_resample_value_in_local_files!=resample.nsim){
      stop("Some files are not generated successfully. Please use `ignore.hasexecuted = FALSE`.")
    }
  } else{

    Genotype_file.Data <- read.csv(Genotype_file)
    Parent_file.Data <- read.csv(Parent_file)[,1]
    Offspring_file.Data <- read.csv(Offspring_file)[,1]
    if(Add_sex) Sex_file.Data <- read.csv(Add_sex.Gender_file)[,1]
    Genotype_file.N_loci <- length(Genotype_file.First_allele_in_column:ncol(Genotype_file.Data))/2
    if(Simulate.Minimum_typed_loci=="defaulted"){
      Simulate.Minimum_typed_loci <- round(Genotype_file.N_loci/2)
    }

    if(round(Genotype_file.N_loci)!=Genotype_file.N_loci){
      stop("How many loci in the Genotype file? Odd number of columns detected. You should use two columns to express an locus.")
    }
    cat("create `.crv` and relevant data files... ")
    for(resample.isim in 1:resample.nsim){
      #basic information
      FileNamePrefix <- paste("Resample",resample.isim) #file prefix
      CurrentTime_TEXT <- paste(strsplit(strsplit(as.character(Sys.time())," CST")[[1]],"-")[[1]],collapse = "/")
      if(length(strsplit(CurrentTime_TEXT,"/0")[[1]])>1){
        CurrentTime_TEXT <- paste(strsplit(CurrentTime_TEXT,"/0")[[1]],collapse = "/")
      } #get current time

      #resampling
      Sampling_Index <- sample(1:nrow(Genotype_file.Data),nrow(Genotype_file.Data)*resample.prop)
      Sampling_Data <- Genotype_file.Data[Sampling_Index,]
      write.csv(Sampling_Data,paste(FileNamePrefix,"Gene.csv"),row.names = FALSE)
      Sampling_IndID <- Genotype_file.Data[Sampling_Index,Genotype_file.ID_in_column]
      if(!Add_sex){
        Simulate.CandParents <- length(Parent_file.Data[Parent_file.Data%in%Sampling_IndID])
        write.csv(Parent_file.Data[Parent_file.Data%in%Sampling_IndID],paste(FileNamePrefix,"Par.csv"),row.names = FALSE)
      } else{
        Sampled_ParentVector <- Parent_file.Data[Parent_file.Data%in%Sampling_IndID]
        Sampled_ParentGender <- Sex_file.Data[as.numeric(lapply(Sampled_ParentVector,function(xx) which(xx==Genotype_file.Data[,Genotype_file.ID_in_column])))]
        Gender_Levels <- names(table(Sampled_ParentGender))
        Is_Male <- Sampled_ParentGender==Gender_Levels[1]
        cat("Read gender information:",Gender_Levels[1],"as male while",Gender_Levels[2],"as female")
        Is_Male[is.na(Is_Male)] <- TRUE
        Sampled_MaleVector <- Sampled_ParentVector[Is_Male]
        Is_Female <- Sampled_ParentGender==Gender_Levels[2]
        Is_Female[is.na(Is_Female)] <- TRUE
        Sampled_FemaleVector <- Sampled_ParentVector[Is_Female]
        Simulate.CandParents <- length(Sampled_FemaleVector)
        Simulate.CandFathers <- length(Sampled_MaleVector)
        write.csv(Sampled_FemaleVector,paste(FileNamePrefix,"Mom.csv"),row.names = FALSE)
        write.csv(Sampled_MaleVector,paste(FileNamePrefix,"Dad.csv"),row.names = FALSE)
      }
      write.csv(Offspring_file.Data[Offspring_file.Data%in%Sampling_IndID],paste(FileNamePrefix,"Off.csv"),row.names = FALSE)

      #calculate prop% of loci typed (null alleles)
      IsNullAllele <- Sampling_Data[,Genotype_file.First_allele_in_column:ncol(Sampling_Data)]==0
      Simulate.Prop_Of_Loci_Typed <- round(1-sum(IsNullAllele)/(nrow(IsNullAllele)*ncol(IsNullAllele)),3)


      #write .crv file
      proj_content <- paste(
        "[ProgramInfo]
ProgramName=Cervus
ProgramVersion=3.0
FileVersion=3.0.7.0

[Registration]
UserName=
UserCompany=
Code=

[FileInfo]
FileName=",Data_directory,"\\",FileNamePrefix," Proj.crv
FileType=.crv
CreationDate=",CurrentTime_TEXT,"

[GenotypeFile]
FileName=",Data_directory,"\\",FileNamePrefix," Gene.csv
HeaderRow=",Genotype_file.HeaderRow,"
ReadLocusNames=1
FirstAlleleColumnNumber=",Genotype_file.First_allele_in_column,"
IDColumnNumber=",Genotype_file.ID_in_column,"
NLoci=",Genotype_file.N_loci,"
PropLociTyped=",Simulate.Prop_Of_Loci_Typed,"
ColumnsPerLocus=2
SexColumn=0
UnknownSexLabel=

[CodecFile]
FileName=
HeaderRow=1
UseSameCodingForAllLoci=1
GenotypeFileName=",Data_directory,"\\",FileNamePrefix," Gene.csv

[AlleleFrequencySummaryFile]
FileName=",Data_directory,"\\",FileNamePrefix," Freq.txt
DoHardyWeinberg=1
HWMinExpectedFrequency=5
UseYatesCorrection=1
UseBonferroniCorrection=1
DoNullAllele=1

[AlleleFrequencyDataFile]
FileName=",Data_directory,"\\",FileNamePrefix," Freq.alf
HeaderRow=1

[SimulationParameters]
AnalysisType=Parent pair (sexes ",ifelse(Add_sex,"known","unknown"),")
NOffspring=",Simulate.offspring,"
NCandidateMales=",ifelse(Add_sex,Simulate.CandFathers,0),"
PropCandidateMalesSampled=",ifelse(Add_sex,Add_sex.Simulate.Prop.father.sampled,0),"
NCandidateFemales=",Simulate.CandParents,"
PropCandidateFemalesSampled=",Simulate.Prop.sampled,"
PropLociTyped=",Simulate.Prop_Of_Loci_Typed,"
PropLociMistyped=",Simulate.Prop.mistyped,"
MinTypedLoci=",Simulate.Minimum_typed_loci,"
CriticalStatisticName=LOD
TruncateAtZero=0
RelaxedConfidence=80
StrictConfidence=95
SimulateInbreeding=0
ParentRelatedness=0
InbreedingRate=0
AlwaysTestSelfing=0
SimulateFemaleRelatives=0
FemalePropRelatives=0
FemaleRelatedTo=Offspring
FemaleRelatedness=0
SimulateMaleRelatives=0
MalePropRelatives=0
MaleRelatedTo=Offspring
MaleRelatedness=0
UseCorrectedLikelihoods=1
UseMistypingRateAsLikelihoodErrorRate=1
LikelihoodErrorRate=0

[SimulationSummaryFile]
FileName=",Data_directory,"\\",FileNamePrefix," Simu.txt

[SimulationDataFile]
FileName=",Data_directory,"\\",FileNamePrefix," Simu.sim

[SimulationOutput]
RepeatSimulation=0
NRepeats=1
ApplyPreviousSimulationData=0
GenerateTables=0
SaveRawStatisticScores=0
GenerateHistograms=0
NCategories=0
MinStatistic=0
MaxStatistic=0

[PreviousSimulationDataFile]
FileName=

[ParentageParameters]
AnalysisType=Parent pair (sexes ",ifelse(Add_sex,"known","unknown"),")
UseSimulationParameters=1
CalculateConfidenceLevels=1
AlwaysTestSelfing=0
MinTypedLoci=",Simulate.Minimum_typed_loci,"
UseCorrectedLikelihoods=1
UseMistypingRateAsLikelihoodErrorRate=0
LikelihoodErrorRate=",Simulate.Prop.mistyped,"
CriticalStatisticName=LOD
TruncateAtZero=0

[OffspringFile]
FileName=",Data_directory,"\\",FileNamePrefix," Off.csv
HeaderRow=",Offspring_file.HeaderRow,"
OffspringIDColumnNumber=",Offspring_file.ID_in_column,"
IncludesKnownParents=0
KnownParentIDColumnNumber=0
IncludesCandidateParents=0
CandidateParentIDColumnNumber=0

[CandidateFemaleFile]
FileName=",Data_directory,"\\",FileNamePrefix," ",ifelse(Add_sex,"Mom","Par"),".csv
HeaderRow=",Parent_file.HeaderRow,"
CandidateParentFormat=One column for all offspring
OffspringIDColumnNumber=0
CandidateParentIDColumnNumber=",Parent_file.ID_in_column,"

[CandidateMaleFile]
FileName=",ifelse(Add_sex,paste(Data_directory,"\\",FileNamePrefix," Dad.csv",sep=""),""),"
HeaderRow=",Parent_file.HeaderRow,"
CandidateParentFormat=One column for all offspring
OffspringIDColumnNumber=0
CandidateParentIDColumnNumber=",Parent_file.ID_in_column,"

[ParentageSummaryFile]
FileName=",Data_directory,"\\",FileNamePrefix," Res.txt

[ParentageDataFile]
FileName=",Data_directory,"\\",FileNamePrefix," Res.csv
OutputType=All parents with positive LOD scores
SortedBy=Joint LOD score
IncludeNonExclusionProbabilities=0

[IdentityParameters]
TestSexesSeparately=0
MinMatchingLoci=0
AllowFuzzyMatching=0
MaxMismatchingLoci=0
ShowAllComparisons=0

[IdentitySummaryFile]
FileName=

[IdentityDataFile]
FileName=

[ConversionSourceFile]
FileName=
HeaderRow=0
AlleleIDLength=2
UseFirstIDAsPopulationName=0
PopulationColumnNumber=0
IDColumnNumber=0
MumIDColumnNumber=0
DadIDColumnNumber=0
FirstAlleleColumnNumber=0
NLoci=0
ContainsAlleleFrequencyData=0

[ConversionSourceAlleleFrequencyDataFile]
FileName=
HeaderRow=0

[ConversionDestinationFile]
FileName=
WriteAlleleFrequencyData=0
AlleleIDLength=2
WritePopulationNames=0
",sep = "")

      writeLines(proj_content,paste(FileNamePrefix,"Proj.crv"))
    }
    cat("done.\n")
  }


  #execute parentage
  setwd(CervusCL_directory)

  if(ignore.hasexecuted){
    local_files <- list.files(path = Data_directory)
    local_files <- local_files[grepl(" Res.csv",local_files)]
    local_files <- local_files[grepl("Resample",local_files)]
    if(length(local_files)>0){
      resample_vector <- (1:resample.nsim)[-as.numeric(sapply(strsplit(local_files," "),function(xx) xx[2]))]

      cat(paste("Start resample -- Remains ",length(resample_vector)," files.\n",sep=""))
    } else{
      resample_vector <- 1:resample.nsim
      cat(paste("Start resample -- Finds ",length(resample_vector)," files.\n",sep=""))
    }
  } else{
    resample_vector <- 1:resample.nsim
    cat(paste("Start resample -- Finds ",length(resample_vector)," files.\n",sep=""))
  }

  cout_progress <- 0
  for(resample.isim in resample_vector){
    cout_progress <- cout_progress + 1
    cat(paste(cout_progress,". ",sep=""))
    FileNamePrefix <- paste("Resample",resample.isim)
    exe_path <- paste(Data_directory,"\\",FileNamePrefix," Proj.crv",sep="")

    is.error <- system(paste('R CMD CervusCL "',exe_path,'" /ALF /O',sep=""),ignore.stdout = ignore.stdout)
    if(is.error) stop(paste('R CMD CervusCL "',exe_path,'" /ALF /O',sep=""))
    is.error <- system(paste('R CMD CervusCL "',exe_path,'" /SIM /O',sep=""),ignore.stdout = ignore.stdout)
    if(is.error) stop(paste('R CMD CervusCL "',exe_path,'" /SIM /O',sep=""))
    is.error <- system(paste('R CMD CervusCL "',exe_path,'" /PAR /O',sep=""),ignore.stdout = ignore.stdout)
    if(is.error) stop(paste('R CMD CervusCL "',exe_path,'" /PAR /O',sep=""))
  }

  cat("done.\n")
  return(NULL)
}
