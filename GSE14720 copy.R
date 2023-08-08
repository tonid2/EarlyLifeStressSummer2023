#Example code: Exploring a Dataset (with more detail)
#Toni Duan
#July 13, 2023

##################################

#Example Dataset:
#GSE14720 

SummarizedExperiment<-gemma.R::get_dataset_object("GSE14720", type = 'se', consolidate="average")
str(SummarizedExperiment)

library(SummarizedExperiment)

#Getting a matrix of expression data
ExpressionData<-assay(SummarizedExperiment[[1]])
str(ExpressionData)
# num [1:10517, 1:11] 4.48 6.26 9.11 7 6.88 ...

#Making a histogram of the log2 cpm values:
hist(ExpressionData)


#The minimum log2 cpm value:
min(ExpressionData)
#[1] 0.849513

#Let's look at the mean vs. variance curve:
ExpressionData_MeanPerGene<-apply(ExpressionData, 1, mean)
ExpressionData_SDPerGene<-apply(ExpressionData, 1, sd)


plot(ExpressionData_SDPerGene~ExpressionData_MeanPerGene)
#This curve is definitely not flat - the log2 transformation did not solve heteroskedasticity issues, especially at the low end

#How many genes have zero variance?
sum(ExpressionData_SDPerGene==0)
#[1] 0

###################

str(rowData(SummarizedExperiment_Filtered[[1]]))

write.csv(rowData(SummarizedExperiment_Filtered[[1]]), "GSE14720Example_Annotation.csv")

str(colData(SummarizedExperiment_Filtered[[1]]))

write.csv(colData(SummarizedExperiment_Filtered[[1]])[,-1], "GSE14720Example_MetaData.csv")

write.csv(ExpressionData_Filtered, "GSE14720ExpressionData_Filtered.csv")

##################


#Reading in a summarized experiment object for the filtered data:

SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE14720", type = 'se', filter=TRUE, consolidate="average")


ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
# num [1:8193, 1:11] 5.59 6.36 9.07 6.74 6.31 ...

hist(ExpressionData_Filtered)


min(ExpressionData_Filtered)
# 0.8494513

#Let's look at the mean vs. variance curve:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)


plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)


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
#Sample data - treatment, disease

rowData(SummarizedExperiment_Filtered[[1]])
#All of the annotation - Probe, GeneSymbol, GeneName, NCBIid

#The distribution of samples from different brain regions
table(SummarizedExperiment_Filtered[[1]]$treatment)
#5-hydroxytryptamine 2A receptor inverse agonist                                          saline 
#                                            4                                               7 

#The distribution of the different phenotypes in this experiment
table(SummarizedExperiment_Filtered[[1]]$disease)
#maternal separation reference subject role 
#                  4                      7


#When combining together criteria "&" means "AND", "|" (also called the pipe) means "OR"
#"==" means equals, whereas "!=" means doesn't equal
SampleFilter<-
  SummarizedExperiment_Filtered[[1]]$`treatment`=="saline"


SummarizedExperiment_Subset<-SummarizedExperiment_Filtered[[1]][,SampleFilter]

SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 8193 7 


table(SummarizedExperiment_Subset$`treatment`)
# saline 
#      7

table(SummarizedExperiment_Subset$disease)
#maternal separation reference subject role 
#                  4                      3 

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

#If we don't have any outliers to remove, we can just rename our object so it works with the downstream code:
SummarizedExperiment_Subset_noBad<-SummarizedExperiment_Subset

#And then we will need to recreate the ExpressionData Matrix as well:
ExpressionData_Subset_noBad<-assay(SummarizedExperiment_Subset_noBad)
str(ExpressionData_Subset_noBad)
#num [1:8193, 1:7] 5.59 6.36 9.07 6.74 6.31 ...


###########################


#Filtering Genes...Again:

#Now that we have subsetted our data, we have a new problem:
#Gemma filtered out genes that lacked variability in the full dataset
#..but that doesn't mean all of the remaining genes have variability in this particular subset of samples

ExpressionData_Subset_noBad_SDperGene<-apply(ExpressionData_Subset_noBad, 1, sd)
min(ExpressionData_Subset_noBad_SDperGene)
#[1] 0.04299843

sum(ExpressionData_Subset_noBad_SDperGene==0)
#[1] 0
#121 genes have zero variability in the HC dataset.

#These genes are going to cause problems - you can't run stats on a variable with no variability! - let's get rid of them.

#...But there is more:
#We may still have issues with genes that lack any variability *for any particular subgroup of interest* as well.

#This function calculates the sd for each treatment group for a particular gene (row of data):
tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd)
#saline 
#0.4094036 

#We want the minimum sd for all of our groups to not be zero:
min(tapply(ExpressionData_Subset_noBad[1,], SummarizedExperiment_Subset_noBad$treatment, sd))!=0
#[1] TRUE

#... and we need to know that for all rows.
GenesToFilter<-apply(ExpressionData_Subset_noBad, 1, function(y) (min(tapply(y, SummarizedExperiment_Subset_noBad$treatment, sd))!=0))

#How many genes we'll end up keeping
sum(GenesToFilter)
#[1] 8193




SummarizedExperiment_Subset_noBad_Filtered<-SummarizedExperiment_Subset_noBad[GenesToFilter,]
SummarizedExperiment_Subset_noBad_Filtered
#class: SummarizedExperiment 
#dim: 8193 7

#And again, remaking the expression set:
ExpressionData_Subset_noBad_Filtered<-assay(SummarizedExperiment_Subset_noBad_Filtered)
str(ExpressionData_Subset_noBad_Filtered)
# num [1:8193, 1:7] 5.59 6.36 9.07 6.74 6.31 ...

###########################

#Checking for batch confound:

#Since we are down to the subset of samples that we plan to use, this would be a good time to check for confounds
#Unfortunately, processing batches are often unbalanced in regards to variables of interest

#Currently, Gemma has all of the processing batch information lumped into one variable

table(SummarizedExperiment_Subset_noBad_Filtered$disease)
# maternal separation reference subject role 
#4                      3 

#this function breaks apart these "blocks" into specific batch-related variables:
strsplit(SummarizedExperiment_Subset_noBad_Filtered$disease, ":")

#To make that into an easier-to-use data.frame
BatchVariables<-do.call(rbind.data.frame, strsplit(SummarizedExperiment_Subset_noBad_Filtered$disease, ":")) 
str(BatchVariables)

#I'm going to rename these (note - this will be dataset specific)
#colnames(BatchVariables)<-c("Device", "Run", "FlowCell", "Lane")

table(SummarizedExperiment_Subset_noBad_Filtered$disease, BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)
#                       maternal separation reference subject role
#maternal separation                      4                      0
#reference subject role                   0                      3

table(SummarizedExperiment_Subset_noBad_Filtered$disease, BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)
#maternal separation reference subject role
#maternal separation                      4                      0
#reference subject role                   0                      3

#Looks like run and device may be redundant
table(BatchVariables$c..maternal.separation....maternal.separation....maternal.separation..., BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)

table(SummarizedExperiment_Subset_noBad_Filtered$disease, BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)
#Looks also potentially redundant

table(BatchVariables$c..maternal.separation....maternal.separation....maternal.separation..., BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)



table(SummarizedExperiment_Subset_noBad_Filtered$treatment, BatchVariables$c..maternal.separation....maternal.separation....maternal.separation...)
#       maternal separation reference subject role
#saline                   4                      3


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
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))
#PC1 and 2, don't seem to be related to our main batch variable
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$disease))
#But they do seem related to antidepressant treatment, esp. PC1
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset_noBad_Filtered$disease))
#Less related to phenotype

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Subset_noBad_Filtered$disease, las=2, xlab="")

#############################

#Plotting data for a particular gene:

#Getting the expression data for a particular gene:
Pvalb<-assay(SummarizedExperiment_Subset_noBad_Filtered)[rowData(SummarizedExperiment_Subset_noBad_Filtered)$GeneSymbol=="Pvalb",]

#######ERROR##############
#Plotting a boxplot of the expression of Pvalb Log2 CPM across antidepressant groups:
boxplot(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$disease, xlab="", ylab="Log2 CPM", main="Pvalb", las=2)

#We can add jittered data points to our graph so that we can see the values for the individual samples in each group.
#To do this, we use the function stripchart and the exact same y~x formula that we used for the boxplot
#The parameter pch determines the size of the data points.
#The parameter "add" places the points on top of the boxplot that we already created.
stripchart(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$disease, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)

#It looks like antidepressant treatment might decrease PVALB
#We would need to run inferential statistics to learn whether this effect is significant

#If we were looking at the data from a single gene, and our data was truly numeric (which RNA-Seq is not quite...)
#We could run a simple linear regression

#First we would just set up our variables so that they have intuitive reference groups:
SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$disease)
levels(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)
#[1] "maternal separation"    "reference subject role"


#We want saline to be our control (reference) group:
SummarizedExperiment_Subset_noBad_Filtered$treatment_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$treatment_factor, ref="reference subject role")

###########ERROR##########
#And then we could just run linear regression:
summary.lm(lm(Pvalb~SummarizedExperiment_Subset_noBad_Filtered$treatment_factor))
# Call:
#  lm(formula = Pvalb ~ SummarizedExperiment_Subset_noBad_Filtered$treatment_factor)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -1.27666 -0.73710  0.00939  0.49644  1.35772 

# Coefficients:
#  Estimate
# (Intercept)                                                                                               4.28211
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorreference subject role,reference subject role -0.05324
# Std. Error
# (Intercept)                                                                                                 0.22450
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorreference subject role,reference subject role    0.31080
# t value
# (Intercept)                                                                                               19.074
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorreference subject role,reference subject role  -0.171
# Pr(>|t|)
# (Intercept)                                                                                              9.66e-15
# SummarizedExperiment_Subset_noBad_Filtered$treatment_factorreference subject role,reference subject role    0.866

# (Intercept)                                                                                              ***
#   SummarizedExperiment_Subset_noBad_Filtered$treatment_factorreference subject role,reference subject role    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.7446 on 21 degrees of freedom
# Multiple R-squared:  0.001395,	Adjusted R-squared:  -0.04616 
# F-statistic: 0.02934 on 1 and 21 DF,  p-value: 0.8656
######################

#But, of course, running a differential expression analysis for an entire dataset is a little more complicated than that...
#And for RNA-Seq it is even more complicated...

#Making a design matrix:

table(SummarizedExperiment_Subset_noBad_Filtered$disease)
#maternal separation reference subject role 
#4                      3 

str(SummarizedExperiment_Subset_noBad_Filtered$disease)
#chr [1:16]
#currently a character vector

SummarizedExperiment_Subset_noBad_Filtered$disease_factor<-as.factor(SummarizedExperiment_Subset_noBad_Filtered$disease)
levels(SummarizedExperiment_Subset_noBad_Filtered$disease_factor)

SummarizedExperiment_Subset_noBad_Filtered$disease_factor<-relevel(SummarizedExperiment_Subset_noBad_Filtered$disease_factor, ref="reference subject role")
levels(SummarizedExperiment_Subset_noBad_Filtered$disease_factor)
#[1] "reference subject role" "maternal separation"   

library(limma)

colData(SummarizedExperiment_Subset_noBad_Filtered)$disease_factor

design <- model.matrix(~disease_factor, data=colData(SummarizedExperiment_Subset_noBad_Filtered))
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
#(Intercept) disease_factormaternal separation
#Down             0                                 0
#NotSig           6                              8193
#Up            8187                                 0

#Writing out the results into your working directory:
write.fit(efit, adjust="BH", file="14720Limma_results.txt")
write.csv(rowData(SummarizedExperiment_Subset_noBad_Filtered), "14720Annotation_LimmaResults.csv")





