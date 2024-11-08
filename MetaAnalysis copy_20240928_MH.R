#Code for aligning limma results from different datasets
#Toni Duan
#July 27, 2023

###############################

library(plyr)

#I tweaked this code from last year to align using NCBIid within species

AligningRatDatasets<-function(ListOfRatDEResults){
  
  Rat_MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[1]]),ListOfRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_FoldChange_Dfs:")
  print(str(Rat_MetaAnalysis_FoldChange_Dfs))
  
  Rat_MetaAnalysis_FoldChanges<<-join_all(Rat_MetaAnalysis_FoldChange_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_FoldChanges:")
  print(str(Rat_MetaAnalysis_FoldChanges))
  
  Rat_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[4]]),ListOfRatDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_SV_Dfs:")
  print(str(Rat_MetaAnalysis_SV_Dfs))
  
  Rat_MetaAnalysis_SV<<-join_all(Rat_MetaAnalysis_SV_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_SV:")
  print(str(Rat_MetaAnalysis_SV))
  
  rm(Rat_MetaAnalysis_SV_Dfs, Rat_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfRatDEResults<-list(DEResults_GSE14720, DEResults_GSE153043, DEResults_GSE124387)

AligningRatDatasets(ListOfRatDEResults)


###########

AligningMouseDatasets<-function(ListOfMouseDEResults){
  
  Mouse_MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[1]]),ListOfMouseDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_FoldChange_Dfs:")
  print(str(Mouse_MetaAnalysis_FoldChange_Dfs))
  
  Mouse_MetaAnalysis_FoldChanges<<-join_all(Mouse_MetaAnalysis_FoldChange_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_FoldChanges:")
  print(str(Mouse_MetaAnalysis_FoldChanges))
  
  Mouse_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[4]]),ListOfMouseDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_SV_Dfs:")
  print(str(Mouse_MetaAnalysis_SV_Dfs))
  
  Mouse_MetaAnalysis_SV<<-join_all(Mouse_MetaAnalysis_SV_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_SV:")
  print(str(Mouse_MetaAnalysis_SV))
  
  rm(Mouse_MetaAnalysis_SV_Dfs, Mouse_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

ListOfMouseDEResults<-list(DEResults_GSE89692, DEResults_GSE116416)

AligningMouseDatasets(ListOfMouseDEResults)


############

#Code for aligning the rat and mice results:

#We have the ortholog database that we downloaded from Jackson Labs on July 27, 2023
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20230727.R"
MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20230727.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	29778 obs. of  5 variables:


#If there are rat datasets:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	35532 obs. of  8 variables:


Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_SV_wOrthologs)

#If there are rat datasets:
MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)

#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")
MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")


#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:
colnames(MetaAnalysis_FoldChanges)

plot(MetaAnalysis_FoldChanges$GSE89692EarlyLifeStressVS_Ctrl~MetaAnalysis_FoldChanges$GSE116416SeparationInsecurity)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

write.csv(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"), "CorMatrix_ELS.csv")
#There isn't much similarity across conditions here (even the two chronic stress conditions)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
pdf("Heatmap_CorMatrix_ELS.pdf", height=5, width=5)
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))
dev.off()



###############################################################################

#Adapting the old meta-analysis code for our new objects:



#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 


#We can only run a meta-analysis if there are differential expression results from more than one comparison.
#Since I know that the differential expression results from the same study (dataset) are artificially correlated, I would actually prefer that there are results from more than one dataset.

#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)

#Or, since there are a limited number of integer answers (0-3), I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
#     0     1     2     3     4     5 
  # 5247  6642  2678  4698 15751   516

#Let's try running a meta-analysis using genes that were found in at least 4 sets of differential expression results
#Since there are 5 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#2 NA is too many

library(metafor)

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

#Example Usage:
NumberOfComparisons=5
CutOffForNAs=2
#I have 5 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5 
# 5247  6642  2678  4698 15751   516 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	11889 obs. of  8 variables:
#   $ Rat_EntrezGene.ID              : chr  "498097" "114521" "360576" "24628" ...
# $ Mouse_EntrezGene.ID            : chr  "68980" "21665" "237858" "18591" ...
# $ MouseVsRat_EntrezGene.ID       : chr  "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# $ GSE89692EarlyLifeStressVS_Ctrl : num  -0.01747 0.0231 0.09934 0.01031 0.00505 ...
# $ GSE116416SeparationInsecurity  : num  0.00262 0.05223 NA -0.02839 -0.01151 ...
# $ GSE14720SeparationVS_Ctrl      : num  NA -0.1799 0.0552 0.1398 -0.1989 ...
# $ GSE153043EarlyLifeStressVS_Ctrl: num  0.13014 0.05731 0.39797 -0.02019 -0.00453 ...
# $ GSE124387SeparationVS_Ctrl     : num  -0.1897 0.119 0.6058 -0.0259 -0.0677 ...
# NULL

str(metaOutput)

# num [1:11889, 1:6] -0.00205 0.04131 0.27717 -0.00807 -0.01103 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:11889] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)

# Log2FC_estimate         SE       pval       CI_lb
# 68980_498097     -0.002048013 0.03693402 0.95577950 -0.07443735
# 21665_114521      0.041313641 0.02637182 0.11721221 -0.01037417
# 237858_360576     0.277167680 0.12194105 0.02302832  0.03816761
# 18591_24628      -0.008066286 0.04064188 0.84267559 -0.08772290
# 116870_64520     -0.011027201 0.02278374 0.62838981 -0.05568250
# 57869_312258      0.025139793 0.02089312 0.22887691 -0.01580997
# CI_ub Number_Of_Comparisons
# 68980_498097  0.07034133                     4
# 21665_114521  0.09300146                     5
# 237858_360576 0.51616775                     4
# 18591_24628   0.07159033                     5
# 116870_64520  0.03362810                     5
# 57869_312258  0.06608956                     5

tail(metaOutput)

# Log2FC_estimate         SE       pval
# 240066_362845         0.03139048 0.03492010 0.36869368
# 240041_100362641      0.29495189 0.25734659 0.25174259
# 71640_100233177       0.08239618 0.03940649 0.03653442
# 232853_308320         0.07100005 0.02878765 0.01365027
# 101197_312301        -0.03709950 0.02763068 0.17937172
# 243308_304336         0.02658531 0.04766786 0.57703564
# CI_lb      CI_ub Number_Of_Comparisons
# 240066_362845    -0.037051655 0.09983261                     4
# 240041_100362641 -0.209438160 0.79934193                     4
# 71640_100233177   0.005160877 0.15963148                     4
# 232853_308320     0.014577296 0.12742281                     4
# 101197_312301    -0.091254639 0.01705565                     4
# 243308_304336    -0.066841969 0.12001260                     4

########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

library(multtest)

#Let's functionalize it!
FalseDiscoveryCorrection<-function(metaOutput, HOM_MouseVsRat){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  
  TempDF2<-join(TempDF, HOM_MouseVsRat[,c(4:5,9:11)], by="Mouse_EntrezGene.ID", type="left", match="first")
  
  TempDF3<-join(TempDF2, HOM_MouseVsRat[,c(15:16,20:22)], by="Rat_EntrezGene.ID", type="left", match="first")
  
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR_annotated[,9]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

#Example usage:

HOM_MouseVsRat <- read.csv("HOM_MouseVsRat.csv", header = TRUE, row.names = 1)
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)
HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat)
# [1] "metaOutputFDR:"
# num [1:11889, 1:7] -0.00205 0.04131 0.27717 -0.00807 -0.01103 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:11889] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 17
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID
# 21802_24827               24827               21802
# 18417_84588               84588               18417
# 404710_310621            310621              404710
# 17136_29409               29409               17136
# 12814_25654               25654               12814
# 13482_25253               25253               13482
# Log2FC_estimate         SE         pval      CI_lb
# 21802_24827        -0.3456419 0.06479706 9.595494e-08 -0.4726418
# 18417_84588        -0.2965111 0.06084605 1.098407e-06 -0.4157671
# 404710_310621      -0.5352452 0.11765462 5.382672e-06 -0.7658440
# 17136_29409        -0.2550212 0.05903334 1.560631e-05 -0.3707244
# 12814_25654        -0.4809033 0.11305363 2.102076e-05 -0.7024844
# 13482_25253         0.3346969 0.08051491 3.225179e-05  0.1768905
# CI_ub Number_Of_Comparisons         FDR
# 21802_24827   -0.2186420                     4 0.001140808
# 18417_84588   -0.1772550                     5 0.006529483
# 404710_310621 -0.3046464                     4 0.021331529
# 17136_29409   -0.1393180                     5 0.046385858
# 12814_25654   -0.2593223                     4 0.049983168
# 13482_25253    0.4925032                     4 0.055426276
# MouseVsRat_EntrezGene.ID Mouse_Symbol
# 21802_24827                21802_24827         Tgfa
# 18417_84588                18417_84588       Cldn11
# 404710_310621            404710_310621       Iqgap3
# 17136_29409                17136_29409          Mag
# 12814_25654                12814_25654      Col11a1
# 13482_25253                13482_25253         Dpp4
# Mouse_Genetic.Location
# 21802_24827                 Chr6  cM
# 18417_84588                 Chr3  cM
# 404710_310621               Chr3  cM
# 17136_29409                 Chr7  cM
# 12814_25654                 Chr3  cM
# 13482_25253                 Chr2  cM
# Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 21802_24827                                   Chr6:86172205-86252701(+)
# 18417_84588                                   Chr3:31204069-31218473(+)
# 404710_310621                                 Chr3:87989309-88028355(+)
# 17136_29409                                   Chr7:30598601-30614298(-)
# 12814_25654                                 Chr3:113824189-114014367(+)
# 13482_25253                                   Chr2:62160417-62242575(-)
# Mouse_Name
# 21802_24827                  transforming growth factor alpha
# 18417_84588                                        claudin 11
# 404710_310621 IQ motif containing GTPase activating protein 3
# 17136_29409                    myelin-associated glycoprotein
# 12814_25654                        collagen, type XI, alpha 1
# 13482_25253                             dipeptidylpeptidase 4
# Rat_Symbol Rat_Genetic.Location
# 21802_24827         Tgfa             Chr4 q34
# 18417_84588       Cldn11             Chr2 q24
# 404710_310621     Iqgap3             Chr2 q34
# 17136_29409          Mag             Chr1 q21
# 12814_25654      Col11a1             Chr2 q41
# 13482_25253         Dpp4             Chr3 q21
# Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 21802_24827                                                        NA
# 18417_84588                                                        NA
# 404710_310621                                                      NA
# 17136_29409                                                        NA
# 12814_25654                                                        NA
# 13482_25253                                                        NA
# Rat_Name
# 21802_24827                  transforming growth factor alpha
# 18417_84588                                        claudin 11
# 404710_310621 IQ motif containing GTPase activating protein 3
# 17136_29409                    myelin-associated glycoprotein
# 12814_25654                    collagen type XI alpha 1 chain
# 13482_25253                             dipeptidylpeptidase 4
############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:17)]
# [1] "Tgfa"     "Cldn11"   "Iqgap3"   "Mag"      "Col11a1" 
# [6] "Dpp4"     "Plaat5"   "Henmt1"   "Mal"      "Klk6"    
# [11] "Ugcg"     "Pex7"     "Fa2h"     "Tmem88b"  "Amn1"    
# [16] "Cfap70"   "Ppp1r14a" 

metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:17)]
# [1] "21802"  "18417"  "404710" "17136"  "12814"  "13482" 
# [7] "66727"  "66715"  "17153"  "19144"  "22234"  "18634" 
# [13] "338521" "320587" "232566" "76670"  "68458"  

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -0.5 to 0.5, but there are a few with Log2FC as big as -1.5-1.5

#MH: For the paper, I updated this code to label with GeneSymbols for convenience
MakeForestPlots<-function(EntrezIDAsCharacter, GeneSymbol){
  
  pdf(paste("ForestPlot_", EntrezIDAsCharacter, GeneSymbol, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
  
  mtext(paste(GeneSymbol), line=-1.5, cex=2)
  dev.off()
}


# [1] "21802"  "18417"  "404710" "17136"  "12814"  "13482" 
# [7] "66727"  "66715"  "17153"  "19144"  "22234"  "18634" 
# [13] "338521" "320587" "232566" "76670"  "68458"  


ToPlot<-metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:17)]
SymbolsToPlot<-metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:17)]
#Example Usage:
for(i in c(1:length(ToPlot))){
MakeForestPlots(ToPlot[i], SymbolsToPlot[i])
}
warnings()

#The myelin-related gene results and Tfga are really pretty

#I'm also going to plot FST because it is the most significant gene in the MDD pathway that is enriched with ELS-DE
MakeForestPlots("14313", "Fst")
#Looks like the results for that one are mostly driven by one study