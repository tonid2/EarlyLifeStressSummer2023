#Example code: Exploring a Dataset (with more detail)
#Toni Duan
#July 13, 2023

##################################

#Example Dataset:
#GSE116416 

SummarizedExperiment<-gemma.R::get_dataset_object("GSE116416", type = 'se')
str(SummarizedExperiment)

library(SummarizedExperiment)

#Getting a matrix of expression data
ExpressionData<-assay(SummarizedExperiment[[1]])
str(ExpressionData)
#num [1:36090, 1:47] 6.49 7.25 6.4 5.67 6.14 ...

#Making a histogram of the log2 cpm values:
hist(ExpressionData)

#The minimum log2 cpm value:
min(ExpressionData)
#[1] NaN

sum(is.na(ExpressionData))
# 36090

head(ExpressionData)

tail(ExpressionData)

#Let's look at the mean vs. variance curve:
ExpressionData_MeanPerGene<-apply(ExpressionData, 1, mean)
ExpressionData_SDPerGene<-apply(ExpressionData, 1, sd)


plot(ExpressionData_SDPerGene~ExpressionData_MeanPerGene)
#This curve is definitely not flat - the log2 transformation did not solve heteroskedasticity issues, especially at the low end

#How many genes have zero variance?
sum(ExpressionData_SDPerGene==0)
#[1] N/A

################

str(rowData(SummarizedExperiment_Filtered[[1]]))

write.csv(rowData(SummarizedExperiment_Filtered[[1]]), "GSE116416Example_Annotation.csv")

str(colData(SummarizedExperiment_Filtered[[1]]))

write.csv(colData(SummarizedExperiment_Filtered[[1]])[,-1], "GSE116416Example_MetaData.csv")

write.csv(ExpressionData_Filtered, "GSE116416ExpressionData_Filtered.csv")


################
#Reading in a summarized experiment object for the filtered data:

SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE116416", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment_Filtered

# #Library Size:
# #Within the Gemma Diagnostics tab, download the MultiQC report
# #Within the MultiQC report, choose "copy table"
# #Paste the table within an excel file and save as a .csv file in your working directory
# GSE116416_LibrarySize<-read.csv("GSE116416_LibrarySize.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# #The library size information is annotated with sample information using GSM number. 
# head(GSE116416_LibrarySize)
# 
# #As far as I can tell, this identifier is not in the Summarized experiment object
# #e.g.,
# colData(SummarizedExperiment_Filtered[[1]])
# #Sample names are given in the format
# #"Sample 18: AMY_resilient_saline"
# row.names(colData(SummarizedExperiment_Filtered[[1]]))
# 
# #So we need to read in the experimental design info from the Gemma website, which contains GSM#
# #To do this, go to the Experimental Design tab -> "show details"->"download design file"
# #It is a tab-delimited text file
# GSE116416_expdesign<-read.delim("15023_GSE116416_expdesign.data.txt", sep="\t", comment.char="#", header=TRUE, stringsAsFactors = FALSE)
# str(GSE116416_expdesign)
# 
# #This column has the GSM #
# GSE116416_expdesign$ExternalID
# 
# #This column includes the Name... and a bunch of other stuff:
# GSE116416_expdesign$Bioassay
# # [1] "GSE116416_Biomat_64___BioAssayId=450411Name=Sample102.NAC_susceptible_ketamine_non_responder"
# # [2] "GSE116416_Biomat_58___BioAssayId=450412Name=Sample108.HIP_susceptible_saline"
# #This splits these long names up at the "Name=" and then combines things back into a matrix again
# temp<-do.call(rbind.data.frame, strsplit(GSE116416_expdesign$Bioassay, "Name="))
# #This adds the name information to the expdesign matrix:
# GSE116416_expdesign$SampleName<-temp[,2]
# GSE116416_expdesign$SampleName
# 
# #For this dataset, The formatting is still a little different from the Summarized Experiment object, e.g.
# GSE116416_expdesign$SampleName[1]
# 
# row.names(colData(SummarizedExperiment_Filtered[[1]]))[1]
# 
# #the ": " version is straight off of GEO, so I'm going to assume that for some reason Gemma doesn't like : in the design data.frame
# #This is likely to be a dataset specific issue - I'm not sure how to generalize this code.
# #For many datasets, things may be fine at this point.
# #So we essentially need a find and replace - in R, these are often done with sub() or gsub()
# #https://www.digitalocean.com/community/tutorials/sub-and-gsub-function-r
# #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/grep
# #this is a little funky, because we're not replacing simple alphanumeric characters, but punctuation that also has a meaning for coding
# temp<-gsub("[.]", ": ", GSE116416_expdesign$SampleName)
# 
# #But it turns out that the names also differ because of an extra space.
# #To take care fo this we just removed all spaces from both vectors because it is easier
# GSE116416_expdesign$SampleName_toJoin<-gsub(" ", "", temp)
# SampleName_toJoin<-gsub(" ", "", row.names(colData(SummarizedExperiment_Filtered[[1]])))
# 
# #We can't just add library size as a column to experimental design because the samples aren't in the same order
# #We have to use a "join" or "merge" function to align them instead
# #To join the design matrix with the library size, we need to have two data.frames that have columns that have the same name:
# str(GSE116416_LibrarySize)
# #In this data.frame, the GSM# is called Sample.Name, whereas in the exp. design data.frame it is called ExternalID
# colnames(GSE116416_LibrarySize)[1]<-"ExternalID"
# 
# #Now we can join by these columns (the function "merge" also works for this):
# GSE116416_expdesign_wLibrarySize<-join(GSE116416_expdesign, GSE116416_LibrarySize, by="ExternalID", type="left")
# str(GSE116416_expdesign_wLibrarySize)
# 
# str(GSE116416_expdesign)
# 
# GSE116416_expdesign$ExternalID%in%GSE116416_LibrarySize$ExternalID
# 
# #Now we have to add this information into the Summarized Experiment object
# #Which also has a different sample order
# #I'm going to grab that sample order information and name it after the column with Sample Names in the experimental design object
# 
# SamplesOrderedLikeSummarizedExperiment<-data.frame(SampleName_toJoin=SampleName_toJoin)
# 
# GSE116416_expdesign_wLibrarySize_Ordered<-join(SamplesOrderedLikeSummarizedExperiment, GSE116416_expdesign_wLibrarySize, by="SampleName_toJoin", type="left")
# str(GSE116416_expdesign_wLibrarySize_Ordered)
# #Double checking that things are actually in the same order
# cbind(row.names(colData(SummarizedExperiment_Filtered[[1]])), GSE116416_expdesign_wLibrarySize_Ordered$SampleName_toJoin)
# 
# #Adding library size to the Summarized Experiment object
# colData(SummarizedExperiment_Filtered[[1]])$LibrarySize<-GSE116416_expdesign_wLibrarySize_Ordered$M.Aligned
# 
# colnames(SummarizedExperiment_Filtered[[1]])
# 
# SummarizedExperiment_Filtered_Temp<-SummarizedExperiment_Filtered[[1]][, colnames(SummarizedExperiment_Filtered[[1]])!="mouse26_f_control_2_rep1"]
# SummarizedExperiment_Filtered_Temp
# SummarizedExperiment_Filtered
# 
# colnames(SummarizedExperiment_Filtered_Temp)
# 
# SummarizedExperiment_Filtered<-SummarizedExperiment_Filtered_Temp

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered)
str(ExpressionData_Filtered)


hist(ExpressionData_Filtered)


min(ExpressionData_Filtered)
#4.53942

#Let's look at the mean vs. variance curve:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)


plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)
#This curve is more flat, but still has heteroskedasticity issues, especially at the low end
#The voom function will help this data still be useable in regression equations

#How many genes have zero variance?
sum(ExpressionData_Filtered_SDPerGene==0)
#[1] 0
# 

#Unfortunately, if we subset the data (e.g., to focus on a region or specific groups) we may still have genes with zero variance 
#So we'll have to come back and filter some more later

###########################

#Dataset Subsetting:

#Before we do much more with the dataset, let's subset down to the samples that we actually plan to use:

#First, we need to know what we have:

#How to access different parts of the Summarized Experiment object:

colData(SummarizedExperiment_Filtered)
#Sample data - treatment

rowData(SummarizedExperiment_Filtered)
#All of the annotation - Probe, GeneSymbol, GeneName, NCBIid

#The distribution of samples from treatment
table(SummarizedExperiment_Filtered$treatment)
#reference substance role Separation insecurity,lipopolysaccharide 
#                     21                                       25                        

#Subsetting the data to have only hippocampus:
SummarizedExperiment_Subset<-SummarizedExperiment_Filtered

SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 13136 46 


table(SummarizedExperiment_Subset$treatment)


table(SummarizedExperiment_Subset$treatment)
#reference substance role Separation insecurity,lipopolysaccharide 
#                      21                                       25 

#Pulling out a matrix of gene expression data for this subset (to use in functions that require matrices)
ExpressionData_Subset<-assay(SummarizedExperiment_Subset)

###########################

#Outlier Removal:

#This would be a good time to check for outliers and remove them if they are present.

#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
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

#If there was an outlier sample, you could remove it using subsetting similar to above by identifying it's column name
#E.g.
OutlierFilter<-colnames(ExpressionData_Subset)!="mouse31_f_early adversity_2_rep1"& colnames(ExpressionData_Subset)!="mouse6_m_control_1_rep1"
OutlierFilter
# [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
# [18]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
# [35]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset[,OutlierFilter]
SummarizedExperiment_Subset_noBad
# class: SummarizedExperiment 
# dim: 13136 44 

#And then we will need to recreate the ExpressionData Matrix as well:
ExpressionData_Subset_noBad<-assay(SummarizedExperiment_Subset_noBad)
str(ExpressionData_Subset_noBad)
# num [1:13136, 1:44] 8.12 10.87 10.05 7.13 10 ...

###########################


#Filtering Genes...Again:

#Now that we have subsetted our data, we have a new problem:
#Gemma filtered out genes that lacked variability in the full dataset
#..but that doesn't mean all of the remaining genes have variability in this particular subset of samples

ExpressionData_Subset_noBad_SDperGene<-apply(ExpressionData_Subset_noBad, 1, sd)
min(ExpressionData_Subset_noBad_SDperGene)
#[1] 0.03789484

sum(ExpressionData_Subset_noBad_SDperGene==0)
#[1] 0

#These genes are going to cause problems - you can't run stats on a variable with no variability! - let's get rid of them.

#...But there is more:
#We may still have issues with genes that lack any variability *for any particular subgroup of interest* as well.

#This function calculates the sd for each treatment group for a particular gene (row of data):
tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd)
# reference substance role Separation insecurity,lipopolysaccharide 
# 0.1790095                                0.1545023 

#We want the minimum sd for all of our groups to not be zero:
min(tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd))!=0
#[1] TRUE

#... and we need to know that for all rows.
GenesToFilter<-apply(ExpressionData_Subset_noBad, 1, function(y) (min(tapply(y, SummarizedExperiment_Subset_noBad$treatment, sd))!=0))

#How many genes we'll end up keeping
sum(GenesToFilter)
#[1] 13136

SummarizedExperiment_Subset_noBad_Filtered<-SummarizedExperiment_Subset_noBad[GenesToFilter,]
SummarizedExperiment_Subset_noBad_Filtered
#class: SummarizedExperiment 
#dim: 13136 44  

#And again, remaking the expression set:
ExpressionData_Subset_noBad_Filtered<-assay(SummarizedExperiment_Subset_noBad_Filtered)
str(ExpressionData_Subset_noBad_Filtered)
#num [1:13136, 1:44] 8.12 10.87 10.05 7.13 10 ...

###########################

#Checking for batch confound:

#Since we are down to the subset of samples that we plan to use, this would be a good time to check for confounds
#Unfortunately, processing batches are often unbalanced in regards to variables of interest

#Currently, Gemma has all of the processing batch information lumped into one variable

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
# reference substance role Separation insecurity,lipopolysaccharide 
# 20                                       24 
#this function breaks apart these "blocks" into specific batch-related variables:
strsplit(SummarizedExperiment_Subset_noBad_Filtered$treatment, ":")

#To make that into an easier-to-use data.frame
BatchVariables<-do.call(rbind.data.frame, strsplit(SummarizedExperiment_Subset_noBad_Filtered$treatment, ":")) 
str(BatchVariables)

table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$c..reference.substance.role....Separation.insecurity.lipopolysaccharide...)
#                                          reference substance role Separation insecurity,lipopolysaccharide
# reference substance role                                       20                                        0
# Separation insecurity,lipopolysaccharide                        0                                       24paration insecurity,lipopolysaccharide                        0                                       24

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

plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment))

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Subset_noBad_Filtered$treatment, las=2, xlab="")

######################

#But, of course, running a differential expression analysis for an entire dataset is a little more complicated than that...
#And for RNA-Seq it is even more complicated...

#Making a design matrix:

table(SummarizedExperiment_Subset_noBad_Filtered$treatment)
# reference substance role Separation insecurity,lipopolysaccharide 
# 20                                       24
str(SummarizedExperiment_Subset_noBad_Filtered$treatment)
#chr [1:44] "reference substance role" "Separation insecurity,lipopolysaccharide" ...

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$treatment)
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "reference substance role"                 "Separation insecurity,lipopolysaccharide"

SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor, ref="reference substance role")
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "reference substance role"                 "Separation insecurity,lipopolysaccharide"


library(limma)

colData(SummarizedExperiment_Subset_noBad_Filtered)$treatment_factor

design <- model.matrix(~treatment_factor, data=colData(SummarizedExperiment_Subset_noBad_Filtered))
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
#(Intercept) treatment_factorSeparation insecurity,lipopolysaccharide
#Down             0                                                        0
#NotSig           0                                                    13136
#Up           13136                                                        0
#str(efit)

#Writing out the results into your working directory:
write.fit(efit, adjust="BH", file="116416Limma_results.txt")
write.csv(rowData(SummarizedExperiment_Subset_noBad_Filtered), "116416Annotation_LimmaResults.csv")