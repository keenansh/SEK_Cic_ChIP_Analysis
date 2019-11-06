#Calculating z-scores for heatmaps for my merged Bam datasets
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

#load in our consensus peaks, noise peaks, and genome broken into 10 bp windows

load('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/consensus_peaks_cictoOreR_q0.01.gr')
load('/Volumes/SEK_Drive/ChIP-Seq/190129_ChIP_Opto_rep3/peaks_lists/Dmel.dm6.10bp.windows') 
load('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/noise_peaks.gr') 

#We also need to load in the basecounts (remember already in CPM)

dest = '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Basecounts/'
file.names = list.files(dest) %>% .[grep('.basecount',.)]

basecount.ranges <- list()
for (i in 1:length(file.names)){
  basecount.ranges[[i]] <- get(load(paste0(dest,file.names[i])))
}
names(basecount.ranges)<- file.names

#Create a heat map matrix for each condition to show z-scores around a given Cic target

Heatmap.matrix <- list()

for (i in 1:length(file.names)){
  
  #Determine the mean and stdev of reads for each condition within the noise regions (determine background signal)
  b <- basecount.ranges[[1]]
  noise.index=subjectHits(findOverlaps(noise.peaks,b))
  mu_CPM=mean(b$score[noise.index])
  sd_CPM=sd(b$score[noise.index])
  
  #generate a matrix to record basecounts 1kb upstream and downstream of all cic binding sites
  mapmat=matrix(nrow=201,ncol=length(consensus_peaks))
  
  #Find which window of the peak has the max score, so we can center it around the max 
  f=findOverlaps(consensus_peaks,b)
  for (j in 1:length(consensus_peaks)){
    w=which(queryHits(f)==j)
    s=subjectHits(f[w])
    m=which.max(b[s]$score)
    mapmat[201,j]=s[m]
  }
    
  #fill in the rest of the matrix by calculating z score with the noise peak mean and stdev
  for (j in 1:length(consensus_peaks)) {
    k=mapmat[201,j]
    kp=seq((k-100),(k+100))
    mapmat[,j]=(mcols(b[kp])$score-mu_CPM)/sd_CPM
    }
  Heatmap.matrix[[i]]<- mapmat
}
names(Heatmap.matrix) <- file.names

#Add in the gene annotations and save as an excel file, which can be plotted in GraphPad Prism 
for (i in 1:length(file.names)){
  colnames(Heatmap.matrix[[i]])=consensus_peaks$SYMBOL
  write.csv(Heatmap.matrix[[i]],paste0('~/Desktop/ChIP_opto_paper/heatmap_excel_files/',file.names,'.csv'),row.names = FALSE)
}