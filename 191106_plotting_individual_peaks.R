#In this script, I will create a file to plot individual peaks near a selected gene
#Resulting data frame can be imported into Graphpad Prism and plotted. 

library(dplyr)
library(GenomicRanges)

dest = '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Basecounts/'
file.names <- list.files(dest) %>% .[grep('.basecount',.)]

#load in basecounts
basecount.ranges <- list()
for (i in 1:length(file.names)){
  basecount.ranges[[i]] <- get(load(paste0(dest,file.names)[i]))
}
names(basecount.ranges) <- file.names %>% gsub('.basecount','',.)
 
#Specify the chromosome and position where the gene of interest is located
#chr3R:30847931-30857074 - tailless
p=GRanges('chr3R',IRanges(30847931,30857074))

#find overlaps and count total reads within a given window
peak.height <- list()
for (i in 1:length(file.names)){
  b <- basecount.ranges[[i]]
  f <- findOverlaps(p,b)
  peak.height[[i]] <- b$score[subjectHits(f)]
}

#Creqte into a data frame that can easily be plotted
d <- as.data.frame(peak.height)
colnames(d)=names(basecount.ranges)
rownames(d) = start(ranges(b[subjectHits(f)]))

write.csv(d,file='/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/Individual_peaks/tll.csv')

