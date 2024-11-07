#Example code: Exploring a Dataset (with more detail)
#Toni Duan & Megan Hagenauer
#Sept 27 2024 - we accidentally excluded this dataset during the original analysis (misread abstract)

##################################

#Example Dataset:
#GSE124387 

#############

library(gemma.R)
library(tidyverse)

#############

#Reading in a summarized experiment object for the filtered data:
SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE124387", type = 'se', consolidate="average")
#filter=TRUE no longer works
#Error in rbindlist(l, use.names, fill, idcol) : 
# Item 2 has 22 columns, inconsistent with item 1 which has 13 columns. To fill missing columns use fill=TRUE.

#############

#I guess maybe we could just use the DE results provided by Gemma? There isn't much that we would change from their pipeline anyway.
