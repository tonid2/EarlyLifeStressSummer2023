#Preparing our limma results for meta-analysis
#Toni Duan
#July 27 2023

###############

#Start with an empty workspace

###############

#Cleaning up the annotation for our results

setwd("/Users/toniduan/Desktop/R early life stress")

TempResultsJoined<-read.delim("116416Limma_results.txt", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

#################

#Reading in the function

FilteringDEResults_GoodAnnotation<-function(TempResultsJoined){
  
  print("# of rows in results")
  print(nrow(TempResultsJoined))
  
  print("# of rows with missing NCBI annotation:")
  print(sum(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="null"))
  
  print("# of rows with NA NCBI annotation:")
  print(sum(is.na(TempResultsJoined$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(TempResultsJoined$GeneSymbol==""|TempResultsJoined$GeneSymbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', TempResultsJoined$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', TempResultsJoined$GeneSymbol)))
  
  #I only want the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  TempResultsJoined_NoNA<-TempResultsJoined[(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="null")==FALSE & is.na(TempResultsJoined$NCBIid)==FALSE,]
  
  if(length(grep('\\|', TempResultsJoined_NoNA$NCBIid))==0){
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$NCBIid)),]
  }
  
  print("# of rows with good annotation")
  print(nrow(TempResultsJoined_NoNA_NoMultimapped))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")
  
  rm(TempResultsJoined_NoNA)
  
  print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}

#################

#Example of using the function for a dataset

FilteringDEResults_GoodAnnotation(TempResultsJoined)
#[1] "# of rows in results"
#[1] 13140
#[1] "# of rows with missing NCBI annotation:"
#[1] NA
#[1] "# of rows with NA NCBI annotation:"
#[1] 5
#[1] "# of rows with missing Gene Symbol annotation:"
#[1] 5
#[1] "# of rows mapped to multiple NCBI_IDs:"
#[1] 0
#[1] "# of rows mapped to multiple Gene Symbols:"
#[1] 0
#[1] "# of rows with good annotation"
#[1] 13135
#[1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"


#################

#Adapting the Log2FC and SE extraction from last year:

#We need to make a matrix of the Log2FoldChange values for the comparisons of interest. 
#These columns start with "Coef" and are specific to your comparisons of interest


colnames(TempResultsJoined_NoNA_NoMultimapped)

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.treatment_factorSeparation.insecurity.lipopolysaccharide)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.treatment_factorSeparation.insecurity.lipopolysaccharide)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid


#Let's rename our columns to something nicer

ComparisonsOfInterest<-c("GSE116416SeparationInsecurity")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest


###################

#Reading in the Extracting Results function:

ExtractingDEResults<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats){
  
  #Let's rename our columns to something nicer
  
  #We calculate the standard error by dividing the log2FC by the tstat
  TempResultsJoined_NoNA_NoMultimapped_SE<-TempResultsJoined_NoNA_NoMultimapped_FoldChanges/TempResultsJoined_NoNA_NoMultimapped_Tstats
  str(TempResultsJoined_NoNA_NoMultimapped_SE)
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE)^2
  str(TempResultsJoined_NoNA_NoMultimapped_SV)
  
  TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChanges, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstats, SE=TempResultsJoined_NoNA_NoMultimapped_SE, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
  
}


###################


#Example usage for the extracting results function:

ExtractingDEResults(GSE_ID="GSE116416", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
# num [1:13135, 1] 0.0491 0.0367 0.0334 0.0508 0.0571 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# ..$ : chr "GSE116416SeparationInsecurity"
# num [1:13135, 1] 0.00241 0.00135 0.00112 0.00258 0.00326 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# ..$ : chr "GSE116416SeparationInsecurity"
# [1] "Output: Named DEResults_GSE116416"

str(DEResults_GSE116416)
# List of 4
# $ Log2FC: num [1:13135, 1] -0.03009 0.01078 -0.00299 0.07725 -0.05633 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# .. ..$ : chr "GSE116416SeparationInsecurity"
# $ Tstat : num [1:13135, 1] -0.613 0.2934 -0.0894 1.5197 -0.9872 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# .. ..$ : chr "GSE116416SeparationInsecurity"
# $ SE    : num [1:13135, 1] 0.0491 0.0367 0.0334 0.0508 0.0571 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# .. ..$ : chr "GSE116416SeparationInsecurity"
# $ SV    : num [1:13135, 1] 0.00241 0.00135 0.00112 0.00258 0.00326 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:13135] "17427" "20807" "11757" "68870" ...
# .. ..$ : chr "GSE116416SeparationInsecurity"" "68870" ...
# # .. ..$ : chr "GSE116416SeparationInsecurity"

############SKIPPED##############

#Sanity Check:
# 
# #Calculating the SE from the Log2FC and Tstat for the first gene(row)
# 1.386/8.37
# #[1] 0.1655914
# 
# #Calculating the SV for the first gene
# 0.1655914^2
# #[1] 0.02742051
# #Looks good (the decimal differences are probably just due to rounding)
# 
# #Calculating the SE from the Log2FC and Tstat for the second gene(row)
# 0.668/7.09
# #[1] 0.09421721
# 
# #Calculating the SV for the second gene
# 0.09421721^2
# #[1] 0.008876883
# #Looks good (the decimal differences are probably just due to rounding)


#################

#Cleaning up our environment before moving on to the next dataset:

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

#################

#Running the next dataset:

##

#This part of the code will be dataset specific:

#Setting the working directory to where the limma results for the dataset are located:

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("89692_Limma_results.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

###

#This part is not dataset specific:
FilteringDEResults_GoodAnnotation(TempResultsJoined)
#[1] "# of rows in results"
#[1] 28703
#[1] "# of rows with missing NCBI annotation:"
#[1] NA
#[1] "# of rows with NA NCBI annotation:"
#[1] 1
#[1] "# of rows with missing Gene Symbol annotation:"
#[1] 1
#[1] "# of rows mapped to multiple NCBI_IDs:"
#[1] 0
#[1] "# of rows mapped to multiple Gene Symbols:"
#[1] 0
#[1] "# of rows with good annotation"
#[1] 28702
#[1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

colnames(TempResultsJoined_NoNA_NoMultimapped)

##

#This part is (somewhat) dataset specific
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.")
#With the appropriate names for this dataset:

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.treatment_factorELS.reference.subject.role)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.treatment_factorELS.reference.subject.role)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:

ComparisonsOfInterest<-c("GSE89692EarlyLifeStressVS_Ctrl")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE#:

ExtractingDEResults(GSE_ID="GSE89692", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)

##

#This part is not dataset specific:

str(DEResults_GSE89692)

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

##################


#Move on to the next dataset

##

#This part of the code will be dataset specific:

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("153043Limma_results.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

###

#This part is not dataset specific:
FilteringDEResults_GoodAnnotation(TempResultsJoined)
#[1] "# of rows in results"
#[1] 18726
#[1] "# of rows with missing NCBI annotation:"
#[1] 1
#[1] "# of rows with NA NCBI annotation:"
#[1] 0
#[1] "# of rows with missing Gene Symbol annotation:"
#[1] 1
#[1] "# of rows mapped to multiple NCBI_IDs:"
#[1] 1
#[1] "# of rows mapped to multiple Gene Symbols:"
#[1] 0
#[1] "# of rows with good annotation"
#[1] 18724
#[1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

colnames(TempResultsJoined_NoNA_NoMultimapped)

##

#This part is (somewhat) dataset specific
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.")
#With the appropriate names for this dataset:

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.treatment_factorearly.life.stress)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.treatment_factorearly.life.stress)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:

ComparisonsOfInterest<-c("GSE153043EarlyLifeStressVS_Ctrl")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE#:

ExtractingDEResults(GSE_ID="GSE153043", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)

##

#This part is not dataset specific:

str(DEResults_GSE153043)

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

##################

#Next dataset

#Running the next dataset:

##

#This part of the code will be dataset specific:

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("14720Limma_results.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(TempResultsJoined)

###

#This part is not dataset specific:
FilteringDEResults_GoodAnnotation(TempResultsJoined)
#[1] "# of rows in results"
#[1] 8193
#[1] "# of rows with missing NCBI annotation:"
#[1] NA
#[1] "# of rows with NA NCBI annotation:"
#[1] 1
#[1] "# of rows with missing Gene Symbol annotation:"
#[1] 1
#[1] "# of rows mapped to multiple NCBI_IDs:"
#[1] 0
#[1] "# of rows mapped to multiple Gene Symbols:"
#[1] 0
#[1] "# of rows with good annotation"
#[1] 8192
#[1] "Outputted object: TempResultsJoined_NoNA_NoMultimapped"

colnames(TempResultsJoined_NoNA_NoMultimapped)

##

#This part is (somewhat) dataset specific
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.")
#With the appropriate names for this dataset:

TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.disease_factormaternal.separation)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest

TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.disease_factormaternal.separation)

str(TempResultsJoined_NoNA_NoMultimapped_Tstats)

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:

ComparisonsOfInterest<-c("GSE14720SeparationVS_Ctrl")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE#:

ExtractingDEResults(GSE_ID="GSE14720", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)

##

#This part is not dataset specific:

str(DEResults_GSE14720)

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)
