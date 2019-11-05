# get 'basecounts' from ChIP-seq under different ERK activating conditions,
#convert to bigwig, save everything.

library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)

#Load in previously established 10 bp windows and function to calculate basecounts within these windows
load('/Volumes/SEK_Drive/ChIP-Seq/190129_ChIP_Opto_rep3/peaks_lists/Dmel.dm6.10bp.windows') #dm6.10bp- meaning 10 bp window
source('/Volumes/SEK_Drive/ChIP-Seq/Useful_R_Functions/PWM_subfunctions.R', chdir = TRUE)

#first load the bam files into R and pull out the GRanges for each of the reads)
dest <- '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Merged_Bam_Files_BWAMEM/'
file.names = list.files(dest) %>% gsub('.bam','',.)

#Load in and parse your genomic ranges for each condition
#Using the basecount function, determines how many reads overlap with each 10 bp section of the genome.
#Normalize the score per million base pairs

Basecount.ranges <- list()

for (i in 1:length(file.names)){
  Cic.ranges <- paste0(dest,file.names[i]) %>% readGAlignments() %>% GRanges()
  Basecount.ranges[[i]] <- basecount(dm6.10bp,Cic.ranges)
  score(Basecount.ranges[[i]]) = score(Basecount.ranges[[i]])/length(Cic.ranges)*1e+6
}
names(Basecount.ranges) <- file.names

#Save files as bigWigs to visualize in the genome browser

dest1<- '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Basecounts/'

for (i in 1:length(file.names)){
  s <- Basecount.ranges[[i]]
  save(s, file = paste0(dest1, file.names[i],'.basecount'))
  export(s, con = paste0(dest1,file.names[i],'.bigWig'),format='bw')
}
