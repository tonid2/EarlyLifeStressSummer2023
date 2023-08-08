#Code for aligning limma results from different datasets
#Toni Duan
#July 27, 2023

###############################

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

ListOfRatDEResults<-list(DEResults_GSE14720, DEResults_GSE153043)

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
#'data.frame':	35078 obs. of  7 variables:


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

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (even the two chronic stress conditions)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))



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
#   0     1     2     3     4 
#5323  6832  4565 17833   525 

#Let's try running a meta-analysis using genes that were found in at least 2 sets of differential expression results
#Since there are 3 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
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
NumberOfComparisons=4
CutOffForNAs=2
#I have 4 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4 
# 5323  6832  4565 17833   525 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	12155 obs. of  7 variables:
#   $ Rat_EntrezGene.ID              : chr  "498097" "114521" "360576" "24628" ...
# $ Mouse_EntrezGene.ID            : chr  "68980" "21665" "237858" "18591" ...
# $ MouseVsRat_EntrezGene.ID       : chr  "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# $ GSE89692EarlyLifeStressVS_Ctrl : num  -0.01747 0.0231 0.09934 0.01031 0.00505 ...
# $ GSE116416SeparationInsecurity  : num  0.00262 0.05223 NA -0.02839 -0.01151 ...
# $ GSE14720SeparationVS_Ctrl      : num  NA -0.1799 0.0552 0.1398 -0.1989 ...
# $ GSE153043EarlyLifeStressVS_Ctrl: num  0.13014 0.05731 0.39797 -0.02019 -0.00453 ...
# NULL

str(metaOutput)

# num [1:12155, 1:6] 0.00138 0.03614 0.2114 -0.00689 -0.0041 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:12155] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)

# Log2FC_estimate         SE      pval       CI_lb      CI_ub Number_Of_Comparisons
# 68980_498097      0.001381562 0.03726999 0.9704300 -0.07166629 0.07442941                     3
# 21665_114521      0.036142383 0.02724126 0.1845909 -0.01724951 0.08953427                     4
# 237858_360576     0.211399764 0.13358674 0.1135376 -0.05042543 0.47322496                     3
# 18591_24628      -0.006887192 0.04196391 0.8696353 -0.08913495 0.07536057                     4
# 116870_64520     -0.004096038 0.02413680 0.8652453 -0.05140330 0.04321122                     4
# 57869_312258      0.025035450 0.02129893 0.2398218 -0.01670969 0.06678059                     4
tail(metaOutput)

# Log2FC_estimate         SE       pval        CI_lb      CI_ub Number_Of_Comparisons
# 240066_362845         0.03521227 0.03936707 0.37107585 -0.041945770 0.11237030                     3
# 240041_100362641      0.06027072 0.03001554 0.04464503  0.001441337 0.11910010                     3
# 71640_100233177       0.06349881 0.03164877 0.04481791  0.001468360 0.12552926                     3
# 232853_308320         0.07478279 0.03059410 0.01451120  0.014819463 0.13474612                     3
# 101197_312301        -0.03811848 0.02848127 0.18077587 -0.093940755 0.01770379                     3
# 243308_304336         0.02972951 0.05680613 0.60073040 -0.081608455 0.14106748                     3
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
# num [1:12155, 1:7] 0.00138 0.03614 0.2114 -0.00689 -0.0041 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:12155] "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 8
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate         SE         pval      CI_lb
# 18417_84588               84588               18417      -0.3117494 0.06273241 6.712186e-07 -0.4347027
# 22784_366568             366568               22784       0.1800466 0.03641759 7.655717e-07  0.1086694
# 21802_24827               24827               21802      -0.3031851 0.07158143 2.280158e-05 -0.4434821
# 66727_293711             293711               66727      -0.6231745 0.15176764 4.023796e-05 -0.9206336
# 100340_362619            362619              100340       0.3233799 0.07890078 4.157208e-05  0.1687372
# 17136_29409               29409               17136      -0.2469745 0.06082236 4.894952e-05 -0.3661841
# CI_ub Number_Of_Comparisons         FDR MouseVsRat_EntrezGene.ID Mouse_Symbol
# 18417_84588   -0.1887962                     4 0.004652762              18417_84588       Cldn11
# 22784_366568   0.2514237                     4 0.004652762             22784_366568      Slc30a3
# 21802_24827   -0.1628881                     3 0.092384401              21802_24827         Tgfa
# 66727_293711  -0.3257154                     3 0.093319597             66727_293711       Plaat5
# 100340_362619  0.4780226                     3 0.093319597            100340_362619      Smpdl3b
# 17136_29409   -0.1277648                     4 0.093319597              17136_29409          Mag
# Mouse_Genetic.Location Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 18417_84588                 Chr3  cM                                 Chr3:31204069-31218473(+)
# 22784_366568                Chr5  cM                                 Chr5:31243450-31265581(-)
# 21802_24827                 Chr6  cM                                 Chr6:86172205-86252701(+)
# 66727_293711               Chr19  cM                                  Chr19:7589906-7617007(+)
# 100340_362619               Chr4  cM                               Chr4:132460277-132484558(-)
# 17136_29409                 Chr7  cM                                 Chr7:30598601-30614298(-)
# Mouse_Name Rat_Symbol Rat_Genetic.Location
# 18417_84588                                              claudin 11     Cldn11             Chr2 q24
# 22784_366568  solute carrier family 30 (zinc transporter), member 3    Slc30a3             Chr6 q14
# 21802_24827                        transforming growth factor alpha       Tgfa             Chr4 q34
# 66727_293711                  phospholipase A and acyltransferase 5     Plaat5             Chr1 q43
# 100340_362619         sphingomyelin phosphodiesterase, acid-like 3B    Smpdl3b             Chr5 q36
# 17136_29409                          myelin-associated glycoprotein        Mag             Chr1 q21
# Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 18417_84588                                                        NA
# 22784_366568                                                       NA
# 21802_24827                                                        NA
# 66727_293711                                                       NA
# 100340_362619                                                      NA
# 17136_29409                                                        NA
# Rat_Name
# 18417_84588                                      claudin 11
# 22784_366568              solute carrier family 30 member 3
# 21802_24827                transforming growth factor alpha
# 66727_293711          phospholipase A and acyltransferase 5
# 100340_362619 sphingomyelin phosphodiesterase, acid-like 3B
# 17136_29409                  myelin-associated glycoprotein
############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:17)]
# [1] "Cldn11"   "Slc30a3"  "Tgfa"     "Plaat5"   "Smpdl3b"  "Mag"      "Iqgap3"   "Klk6"     "Dpp4"    
# [10] "Thop1"    "Gnaz"     "Henmt1"   "Gpr37"    "Ttc23"    "Ppp1r14a" "Slc39a4"  "Pex7"    
metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:17)]
# [1] "18417"  "22784"  "21802"  "66727"  "100340" "17136"  "404710" "19144"  "13482"  "50492"  "14687" 
# [12] "66715"  "14763"  "67009"  "68458"  "72027"  "18634" 

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -1 to 1, but there are a few with Log2FC as big as -3-3

MakeForestPlots<-function(EntrezIDAsCharacter){
  
  pdf(paste("ForestPlot_", EntrezIDAsCharacter, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
  
  mtext(paste(EntrezIDAsCharacter), line=-1.5, cex=2)
  dev.off()
}


#Example Usage:
MakeForestPlots("17153") 
