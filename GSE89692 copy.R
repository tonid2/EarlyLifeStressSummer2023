#Example code: Exploring a Dataset (with more detail)
#Toni Duan
#July 13, 2023

##################################

#Example Dataset:
#GSE89692

SummarizedExperiment<-gemma.R::get_dataset_object("GSE89692", type = 'se', consolidate="average")
str(SummarizedExperiment)

library(SummarizedExperiment)

#Getting a matrix of expression data
ExpressionData<-assay(SummarizedExperiment[[1]])
str(ExpressionData)
# num [1:25638, 1:163] -2.03 -5.8 -5.8 -5.7 -5.8 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:25638] "100009600" "100009609" "100009614" "100009664" ...
# ..$ : chr [1:163] "FPSS2" "FNSS1" "RNAseq_StdDefeat1" "RNAseq_ELSDefeat2" ...

#Making a histogram of the log2 cpm values:
hist(ExpressionData)

#The minimum log2 cpm value:
min(ExpressionData)
#[1] -9.104648

#Let's look at the mean vs. variance curve:
ExpressionData_MeanPerGene<-apply(ExpressionData, 1, mean)
ExpressionData_SDPerGene<-apply(ExpressionData, 1, sd)

plot(ExpressionData_SDPerGene~ExpressionData_MeanPerGene)
#This curve is definitely not flat - the log2 transformation did not solve heteroskedasticity issues, especially at the low end

#How many genes have zero variance?
sum(ExpressionData_SDPerGene==0)
#[1] 3371

########################
str(rowData(SummarizedExperiment_Filtered[[1]]))

write.csv(rowData(SummarizedExperiment_Filtered[[1]]), "GSE89692_Annotation.csv")

str(colData(SummarizedExperiment_Filtered[[1]]))

write.csv(colData(SummarizedExperiment_Filtered[[1]])[,-1], "GSE89692_MetaData.csv")

write.csv(ExpressionData_Filtered, "GSE89692ExpressionData_Filtered.csv")

########################

#Reading in a summarized experiment object for the filtered data:

SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE89692", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment_Filtered
# $`14636`
# class: SummarizedExperiment 
# dim: 29208 163 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(29208): 100009600 100009664 ... Averaged from 100039542 100993 Averaged from 11641
# 677884
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(163): RNAseq_ELSDefeat2 MPSC1 ... FVSS1 CJP_NAc_18
# colData names(5): factorValues organism part developmental stage treatment block

#Library Size:
#Within the Gemma Diagnostics tab, download the MultiQC report
#Within the MultiQC report, choose "copy table"
#Paste the table within an excel file and save as a .csv file in your working directory
GSE89692_LibrarySize<-read.csv("GSE89692_LibrarySize.csv", header=TRUE, stringsAsFactors = FALSE)

#The library size information is annotated with sample information using GSM number. 
head(GSE89692_LibrarySize)
#Sample.Name X..Alignable X..Aligned M.Aligned X..Dups X..GC Length M.Seqs
#1  GSM2386791       98.40%     85.80%      37.0                          NA
#2  GSM2386792       98.40%     85.70%      32.6                          NA
#3  GSM2386793       98.40%     85.80%      30.2                          NA
#4  GSM2386794       98.40%     85.50%      33.5                          NA
#5  GSM2386795       98.30%     85.70%      39.2                          NA
#6  GSM2386796       98.30%     84.90%      35.6                          NA

 
colData(SummarizedExperiment_Filtered[[1]])

row.names(colData(SummarizedExperiment_Filtered[[1]]))
# [1] "RNAseq_ELSDefeat2" "MPSC1"             "FVES4"             "CJP_NAc_5"         "FVEC3"            
# [6] "MPEC2"             "FPSC3"             "FVEC5"             "FVSC2"             "FVSC5"            
# [11] "CJP_NAc_11"        "RNAseq_ELSDefeat5" "FVSS3"             "FNSC5"             "MHSD4"    

#So we need to read in the experimental design info from the Gemma website, which contains GSM#
#To do this, go to the Experimental Design tab -> "show details"->"download design file"
#It is a tab-delimited text file
GSE89692_expdesign<-read.delim("14636_GSE89692_expdesign.data.txt", sep="\t", comment.char="#", header=TRUE, stringsAsFactors = FALSE)
str(GSE89692_expdesign)
#'data.frame':	163 obs. of  6 variables:

#This column has the GSM #
GSE89692_expdesign$ExternalID
# [1] "GSM3736136" "GSM3736087" "GSM3736037"

#This column includes the Name... and a bunch of other stuff:
GSE89692_expdesign$Bioassay
# [1] "GSE89692_Biomat_164___BioAssayId=465749Name=CJP_NAc_20"      
# [2] "GSE89692_Biomat_115___BioAssayId=465899Name=FNSS4"  
temp<-do.call(rbind.data.frame, strsplit(GSE89692_expdesign$Bioassay, "Name="))
#This adds the name information to the expdesign matrix:
GSE89692_expdesign$SampleName<-temp[,2]
GSE89692_expdesign$SampleName

#For this dataset, The formatting is still a little different from the Summarized Experiment object, e.g.
GSE89692_expdesign$SampleName[1]
# [1] "CJP_NAc_20"
row.names(colData(SummarizedExperiment_Filtered[[1]]))[1]
#[1] "RNAseq_ELSDefeat2"

temp<-gsub("[.]", ": ", GSE89692_expdesign$SampleName)

GSE89692_expdesign$SampleName_toJoin<-gsub(" ", "", temp)
SampleName_toJoin<-gsub(" ", "", row.names(colData(SummarizedExperiment_Filtered[[1]])))

#We can't just add library size as a column to experimental design because the samples aren't in the same order
#We have to use a "join" or "merge" function to align them instead
#To join the design matrix with the library size, we need to have two data.frames that have columns that have the same name:
str(GSE89692_LibrarySize)
#In this data.frame, the GSM# is called Sample.Name, whereas in the exp. design data.frame it is called ExternalID
colnames(GSE89692_LibrarySize)[1]<-"ExternalID"

#Now we can join by these columns (the function "merge" also works for this):
GSE89692_expdesign_wLibrarySize<-join(GSE89692_expdesign, GSE89692_LibrarySize, by="ExternalID", type="left")
str(GSE89692_expdesign_wLibrarySize)


#Now we have to add this information into the Summarized Experiment object
#Which also has a different sample order
#I'm going to grab that sample order information and name it after the column with Sample Names in the experimental design object

SamplesOrderedLikeSummarizedExperiment<-data.frame(SampleName_toJoin=SampleName_toJoin)

GSE89692_expdesign_wLibrarySize_Ordered<-join(SamplesOrderedLikeSummarizedExperiment, GSE89692_expdesign_wLibrarySize, by="SampleName_toJoin", type="left")                   
str(GSE89692_expdesign_wLibrarySize_Ordered)
#Double checking that things are actually in the same order
cbind(row.names(colData(SummarizedExperiment_Filtered[[1]])), GSE89692_expdesign_wLibrarySize_Ordered$SampleName_toJoin)

#Adding library size to the Summarized Experiment object
colData(SummarizedExperiment_Filtered[[1]])$LibrarySize<-GSE89692_expdesign_wLibrarySize_Ordered$M.Aligned

hist(colData(SummarizedExperiment_Filtered[[1]])$LibrarySize)

boxplot(colData(SummarizedExperiment_Filtered[[1]])$LibrarySize~colData(SummarizedExperiment_Filtered[[1]])$treatment)

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
#  num [1:29208, 1:163] -1.32 -5.84 2.54 6.76 -2.02 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:29208] "100009600" "100009664" "100017" "100019" ...
# ..$ : chr [1:163] "RNAseq_ELSDefeat2" "MPSC1" "FVES4" "CJP_NAc_5" ...

hist(ExpressionData_Filtered)


min(ExpressionData_Filtered)
#[1] -9.104648

#Let's look at the mean vs. variance curve:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)

plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)
#no visible difference

#How many genes have zero variance?
sum(ExpressionData_Filtered_SDPerGene==0)
#[1] 0


#Unfortunately, if we subset the data (e.g., to focus on a region or specific groups) we may still have genes with zero variance 
#So we'll have to come back and filter some more later

###########################

#Dataset Subsetting:

#Before we do much more with the dataset, let's subset down to the samples that we actually plan to use:

#First, we need to know what we have:

#How to access different parts of the Summarized Experiment object:

colData(SummarizedExperiment_Filtered[[1]])
#Sample data - organism part, developmental stage, treatment, block, etc

rowData(SummarizedExperiment_Filtered[[1]])
#All of the annotation - Probe, GeneSymbol, GeneName, NCBIid

#The distribution of samples from different brain regions
table(SummarizedExperiment_Filtered[[1]]$`organism part`)
# medial prefrontal cortex        nucleus accumbens      ventral hippocampus   ventral tegmental area 
#                      46                       51                       24                       42                                                                   24
table(SummarizedExperiment_Filtered[[1]]$treatment, SummarizedExperiment_Filtered[[1]]$`developmental stage`, SummarizedExperiment_Filtered[[1]]$`organism part`)

SampleFilter<-
  SummarizedExperiment_Filtered[[1]]$`organism part`=="medial prefrontal cortex" &
  (SummarizedExperiment_Filtered[[1]]$treatment=="ELS,reference subject role"|
     SummarizedExperiment_Filtered[[1]]$treatment=="reference subject role,reference subject role")

#Subsetting the data to have only PFC:
SummarizedExperiment_Subset<-SummarizedExperiment_Filtered[[1]][,SampleFilter]

SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 29208 23
#This dataset is much smaller now

#Sanity Check: Double-checking that the subsetting worked properly:

table(SummarizedExperiment_Subset$`organism part`)
#PFC 
#23

table(SummarizedExperiment_Subset$treatment)
# ELS,reference subject role reference subject role,reference subject role 
#                         11                                            12 

#Pulling out a matrix of gene expression data for this subset (to use in functions that require matrices)
ExpressionData_Subset<-assay(SummarizedExperiment_Subset)
str(ExpressionData_Subset)
#num [1:29208, 1:23] -0.731 -5.705 2.118 6.549 -3.381 ...

###########################

#Outlier Removal:

#This would be a good time to check for outliers and remove them if they are present.

#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
###########ERROR###################
plot(ExpressionData_Subset[,1]~ExpressionData_Subset[,2])

#Creating a matrix showing how each sample correlates with every other sample:
CorMatrix<-cor(ExpressionData_Subset)

#Writing that matrix out to your working directory to save it:
write.csv(CorMatrix, "CorMatrix.csv")

#Creating a hierarchically clustered heatmap illustrating the sample-sample correlation matrix:
heatmap(CorMatrix)
#Cool - there is actually some clustering amongst the samples that are from the antidepressant vs. saline groups

#Creating a boxplot illustrating the sample-sample correlations for each sample. Outliers should be obvious at this point.
boxplot(CorMatrix)
#None of these samples look like clear outliers

# SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset[,OutlierFilter]
# SummarizedExperiment_Subset_noBad
# # class: SummarizedExperiment 
# # dim: 29208 23 

#If we don't have any outliers to remove, we can just rename our object so it works with the downstream code:
SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset

#And then we will need to recreate the ExpressionData Matrix as well:
ExpressionData_Subset_noBad<-assay(SummarizedExperiment_Subset_noBad)
str(ExpressionData_Subset_noBad)
#num [1:29208, 1:23] -0.731 -5.705 2.118 6.549 -3.381 ...
colnames(ExpressionData_Subset_noBad)

FemaleOrNot<-grep("[F]", colnames(ExpressionData_Subset_noBad))

Sex<-colnames(ExpressionData_Subset_noBad)

Sex[FemaleOrNot]<-"F"
Sex[-(FemaleOrNot)]<-"M"

SummarizedExperiment_Subset_noBad$Sex<-Sex
###########################


#Filtering Genes...Again:

#Now that we have subsetted our data, we have a new problem:
#Gemma filtered out genes that lacked variability in the full dataset
#..but that doesn't mean all of the remaining genes have variability in this particular subset of samples

ExpressionData_Subset_noBad_SDperGene<-apply(ExpressionData_Subset_noBad, 1, sd)
min(ExpressionData_Subset_noBad_SDperGene)
#[1] 0

sum(ExpressionData_Subset_noBad_SDperGene==0)
#[1] 428

#These genes are going to cause problems - you can't run stats on a variable with no variability! - let's get rid of them.

#...But there is more:
#We may still have issues with genes that lack any variability *for any particular subgroup of interest* as well.

#This function calculates the sd for each treatment group for a particular gene (row of data):
tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd)
# ELS,reference subject role reference subject role,reference subject role 
#                  0.5354965                                     0.4970677 

#We want the minimum sd for all of our groups to not be zero:
min(tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd))!=0
#[1] TRUE

#... and we need to know that for all rows.
GenesToFilter<-apply(ExpressionData_Subset_noBad, 1, function(y) (min(tapply(y, SummarizedExperiment_Subset_noBad$treatment, sd))!=0))

#How many genes we'll end up keeping
sum(GenesToFilter)
#[1] 28703
#So we're not losing that many

SummarizedExperiment_Subset_noBad_Filtered<-SummarizedExperiment_Subset_noBad[GenesToFilter,]
SummarizedExperiment_Subset_noBad_Filtered
#class: SummarizedExperiment 
#dim: 28703 23

#And again, remaking the expression set:
ExpressionData_Subset_noBad_Filtered<-assay(SummarizedExperiment_Subset_noBad_Filtered)
str(ExpressionData_Subset_noBad_Filtered)
#num [1:28703, 1:23] -0.731 -5.705 2.118 6.549 -3.381 ...

###########################

#Checking for batch confound:

#Since we are down to the subset of samples that we plan to use, this would be a good time to check for confounds
#Unfortunately, processing batches are often unbalanced in regards to variables of interest

#Currently, Gemma has all of the processing batch information lumped into one variable

table(SummarizedExperiment_Subset_noBad_Filtered$block)
#Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=1 Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=2 
#3                                                   5 
#Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=3 Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=4 
#2                                                   2 
#Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=5 Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=7 
#1                                                   2 
#Device=HWI-D00381:Run=363:Flowcell=C8DHEANXX:Lane=8 Device=HWI-D00381:Run=364:Flowcell=C8DV7ANXX:Lane=1 
#2                                                   2 
#Device=HWI-D00381:Run=364:Flowcell=C8DV7ANXX:Lane=2 Device=HWI-D00381:Run=364:Flowcell=C8DV7ANXX:Lane=7 
#2                                                   2 

#this function breaks apart these "blocks" into specific batch-related variables:
strsplit(SummarizedExperiment_Subset_noBad_Filtered$block, ":")

#To make that into an easier-to-use data.frame
BatchVariables<-do.call(rbind.data.frame, strsplit(SummarizedExperiment_Subset_noBad_Filtered$block, ":")) 
str(BatchVariables)
#I'm going to rename these (note - this will be dataset specific)
colnames(BatchVariables)<-c("Device", "Run", "FlowCell", "Lane")

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Device)
#                                              Device=HWI-D00381
#ELS,reference subject role                                   11
#reference subject role,reference subject role                12


table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Run)
#                                              Run=363 Run=364
#ELS,reference subject role                          8       3
#reference subject role,reference subject role       9       3

table(SummarizedExperiment_Subset_noBad$Sex, BatchVariables$Run)
#Run=363 Run=364
#F      11       0
#M       6       6

#Looks like run and device may be redundant
table(BatchVariables$Device, BatchVariables$Run)
#                  Run=363 Run=364
#Device=HWI-D00381      17       6

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$FlowCell)
#Looks also potentially redundant
#                                              Flowcell=C8DHEANXX Flowcell=C8DV7ANXX
#ELS,reference subject role                                     8                  3
#reference subject role,reference subject role                  9                  3

table(BatchVariables$FlowCell, BatchVariables$Run)
#                   Run=363 Run=364
#Flowcell=C8DHEANXX      17       0
#Flowcell=C8DV7ANXX       0       6

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$Lane)
#                                              Lane=1 Lane=2 Lane=3 Lane=4 Lane=5 Lane=7 Lane=8
#ELS,reference subject role                         2      3      1      1      1      2      1
#reference subject role,reference subject role      3      4      1      1      0      2      1
#Too many Lanes to take into consideration

###########################
table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#ELS,reference subject role reference subject role,reference subject role 
#                        11                                            12 

#I'm going to rename the control condition because it is a mouthful
SummarizedExperiment_Subset_noBad_Filtered$treatment[SummarizedExperiment_Subset_noBad_Filtered$treatment=="reference substance role,saline"]<-"reference subject role,reference subject role "

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#ELS,reference subject role reference subject role,reference subject role 
#                        11                                            12 

SummarizedExperiment_Subset_noBad_Filtered$treatment[SummarizedExperiment_Subset_noBad_Filtered$phenotype=="ELS,reference subject role"]<-"ELS,reference subject role"
SummarizedExperiment_Subset_noBad_Filtered$treatment[SummarizedExperiment_Subset_noBad_Filtered$phenotype=="reference subject role,reference subject role"]<-"reference subject role,reference subject role"

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#ELS,reference subject role reference subject role,reference subject role 
#.                       11                                            12 

###########################

#Principal components analysis:

library(stats)
pca_output<-prcomp(t(ExpressionData_Subset_noBad_Filtered), scale=TRUE)

PCeigenvectors<-pca_output$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, rowData(SummarizedExperiment_Subset_noBad_Filtered))
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]
PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]

#Output a scree plot for the PCA (no outliers):
#This plot illustrates the proportion of variance explained by each principal component (PC):
png("10 PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis")
dev.off()

#You can color these plots in using different variables in the dataset 
#This can help explain the main sources of variation in the data
#Technical variables (brain region, batch) tend to be the main culprits
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(BatchVariables$Run ))

plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))
#PC1 and 2, don't seem to be related to our main batch variable
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$LibrarySize))

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$LibrarySize))

plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad$Sex))

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$block))

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad$Sex))

plot(PC1~SummarizedExperiment_Subset_noBad_Filtered$LibrarySize)

plot(PC2~SummarizedExperiment_Subset_noBad_Filtered$LibrarySize)

plot(PC3~SummarizedExperiment_Subset_noBad_Filtered$LibrarySize)

plot(PC4~SummarizedExperiment_Subset_noBad_Filtered$LibrarySize)

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Subset_noBad_Filtered$treatment, las=2, xlab="")

boxplot(SummarizedExperiment_Subset_noBad_Filtered$LibrarySize~SummarizedExperiment_Subset_noBad_Filtered$treatment, las=2, xlab="")

boxplot(SummarizedExperiment_Subset_noBad_Filtered$LibrarySize~SummarizedExperiment_Subset_noBad_Filtered$Sex, las=2, xlab="")

#############################

#Making a design matrix:

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#ELS,reference subject role reference subject role,reference subject role 
#11                                            12 
table(SummarizedExperiment_Subset_noBad_Filtered$Sex)
#F  M 
#11 12                              

str(SummarizedExperiment_Subset_noBad_Filtered$treatment)
# chr [1:23] "reference subject role,reference subject role" "ELS,reference subject role" ...

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment)
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "ELS,reference subject role"                    "reference subject role,reference subject role"

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor, ref="reference subject role,reference subject role")
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "reference subject role,reference subject role" "ELS,reference subject role"    

str(SummarizedExperiment_Subset_noBad_Filtered$Sex)
#chr [1:23]
#currently a character vector

SummarizedExperiment_Subset_noBad_Filtered$Sex_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$Sex)
levels(SummarizedExperiment_Subset_noBad_Filtered$Sex_factor)



library(limma)

colData(SummarizedExperiment_Subset_noBad_Filtered)$treatment_factor

design <- model.matrix(~treatment_factor+LibrarySize+Sex_factor, data=colData(SummarizedExperiment_Subset_noBad_Filtered))
design

#testing it out with the data for 1 gene 
#we have to tell the lm function to not add an intercept
summary.lm(lm(ExpressionData_Subset_noBad_Filtered[1,]~0+design))


#Applying the model to all genes using limma:

fit <- lmFit(ExpressionData_Subset_noBad_Filtered, design)

#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
efit <- eBayes(fit, trend=TRUE)

dt<-decideTests(efit)
summary(decideTests(efit))
#(Intercept) treatment_factorELS,reference subject role LibrarySize Sex_factorM
#Down          8845                                          0           0          88
#NotSig        6930                                      28703       28703       28389
#Up           12928                                          0           0         226
#str(efit)

#Writing out the results into your working directory:
write.fit(efit, adjust="BH", file="GSE89692Limma_results.txt")
write.csv(rowData(SummarizedExperiment_Subset_noBad_Filtered), "GSE89692Annotation_LimmaResults.csv")


