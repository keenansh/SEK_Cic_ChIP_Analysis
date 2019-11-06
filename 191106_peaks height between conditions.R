#Determine how height of the peaks change across conditions 
#How many peaks have a fold change in height greater than 2? 

library(GenomicRanges)
library(rtracklayer)

#Load in consensus peaks list from MACS2 
load('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/consensus_peaks_cictoOreR_q0.01.gr')

#Load in basecount files 
dest = '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Basecounts/'
file.names <- list.files(dest) %>% .[grep('.basecount',.)] %>% gsub('.basecount','',.)

basecount.ranges <- list()
for (i in 1:length(file.names)){
  basecount.ranges[[i]] <- get(load(paste0(dest,file.names)[i]))
}
names(basecount.ranges) <- file.names %>% gsub('.basecount','',.)

#Calculate the maximum height of each peak in each condition, based on read count
#load in max.read.matrix function 
source('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Rscripts/max.reads.matrix.R')

r=lapply(basecount.ranges,max.reads.matrix,peaks=consensus_peaks)
fin.reads=matrix(unlist(r), ncol = length(file.names), byrow = FALSE)
colnames(fin.reads)<- file.names

#Find which ones have a height less than half compared to the WT CicsfGFP signal
fol.cha=matrix(data=NA,nrow=length(consensus_peaks),ncol=ncol(fin.reads)-1)
for (i in 1:length(consensus_peaks)){
  fol.cha[i,]=fin.reads[i,2:ncol(fin.reads)]/fin.reads[i,1]
}

percent.diff=matrix(data=NA,nrow=1,ncol=ncol(fin.reads)-1)
for (i in 1:8){
  w=which(fol.cha[,i]<=0.5)
  s=sum(as.logical(w))
  percent.diff[1,i]=s/103*100
}