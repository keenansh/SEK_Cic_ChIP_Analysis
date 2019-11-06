#Generating noise peaks to later calculate z-score

library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)

#Load in consensus peaks from MACS2, 10 bp windows of genome, and helpful functions

load('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/consensus_peaks_cictoOreR_q0.01.gr')
load('/Volumes/SEK_Drive/ChIP-Seq/Useful peaks lists/Dmel.dm6.10bp.windows')
source('/Volumes/SEK_Drive/ChIP-Seq/Useful_R_Functions/PWM_subfunctions.R', chdir = TRUE)
seqlevels(dm6.10bp,pruning.mode='coarse')<-c('chr2L','chr2R','chr3L','chr3R','chrX','chr4')

#resize all peaks to be 20 kb. The trim function gets rid of regions that fall off the chromosome 
noise.peaks=trim(resize(consensus_peaks,fix='center',width=20000))

#Use the reduce function. This function looks for overlaps in peaks and then merges them to get a simplified set 
noise.peaks=reduce(noise.peaks)

#Now we want to use the 'gaps' function, which will find all the gaps between our ranges
noise.peaks=gaps(noise.peaks)

#only keep strand * 
noise.peaks = noise.peaks[strand(noise.peaks) %in% '*']

#The findOverlaps function finds the indices of 10 bp windows overlapping the set of noise ranges
#this gives us our noise regions and we want to sample N of the noise regions 
#N will be 25,000 in this case (number of background regions we want)
n=subjectHits(findOverlaps(noise.peaks,dm6.10bp))
set.seed(8)
s=sample(seq(1,length(n),1),25000)
s=sort(n[s])
#s is a sample of 25,000 of the noise regions 

# this models a set of noise peaks whose centers are the regions we sampled
# above, and whose widths are randomly distributed over the mean and sd of the
# log2 widths of the true peaks list. The log2 widths of the peaks approximates
# a normal distribution.

wi=rnorm(25000,mean=mean(log2(width(consensus_peaks))),sd=sd(log2(width(consensus_peaks))))
noise.peaks=trim(resize(dm6.10bp[s], width=2^wi,fix='center'))

# sometimes the peaks we've randomly picked overlap with one another.
# Reduce them to single intervals rather than count them twice.

noise.peaks=reduce(noise.peaks[!as.logical(countOverlaps(noise.peaks,consensus_peaks))])

#Save this as the official noise regions
save(noise.peaks, file = '/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/noise_peaks.gr')

