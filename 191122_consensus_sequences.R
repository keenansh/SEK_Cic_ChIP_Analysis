# Getting sequencs of cic peaks and determining how many Cic binding sites are within each peak. 

library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(stringr)

#load in your cic peaks (consensus_peaks)
#Lets also expands the peaks a bit because some of them may be cut off

load('/Volumes/SEK_Drive/ChIP-Seq/CicsfGFP_Replicates_Analysis/consensus_peaks_cictoOreR_q0.01.gr')
start(consensus_peaks)=start(consensus_peaks)-250
end(consensus_peaks)=end(consensus_peaks)+250

cic_seqs = getSeq(Dmelanogaster, consensus_peaks)
cic_seqs = as.character(cic_seqs)

#We can specify our own consensus sequences from the literature.  
#But the hannomers function can also tell us the top sites in out data. 

source('/Volumes/SEK_Drive/ChIP-Seq/Useful_R_Functions/hannonmers.dm6.R', chdir = TRUE)

cic.8mers = hannonmers(8, consensus_peaks,1)
cic.8mer.pval = cic.8mers[[1]]
cic.8mer.pval = cic.8mer.pval[order(cic.8mer.pval, decreasing = FALSE)]

#Define what are Cic consensus sequences: 
#For now we are looking at the strongest optimal consensus sequence 

opt.CBS = c('TGAATGAA','TTCATTCA','TCAATGAA','TTCATTGA')

#scan each binding site and calculate total strong binding sites 

opt.counts = sapply(opt.CBS, str_count, string=cic_seqs)
opt.counts = cbind(opt.counts, 'Total'=rowSums(opt.counts))

hist(opt.counts[,'Total'])

#we can also scan to check for weak binding sites in addition to strong binding sites 

all.CBS = c('TGAATGAA','TTCATTCA','TCAATGAA','TTCATTGA', 'TGAACGAA','TTCGTTCA','TTCTTTGA',
            'TCAAAGAA','TTTATTGA','TCAATAAA','TGAATGCA','TGCATTCA','TTAATGAA','TTCATTAA',
            'TCCATGAA','TTCATGAA','TGAATGGA','TCAATTCA',
            'TTCATTCT','AGAATGAA','TGAATGAG','CTCATTCA','TGAGTGAA','TTCACTCA','TCAATTAA','TTAATTGA',
            'TCAATCAA','TTGATTGA','TCAATGCA','TGCATTGA',
            'TCAATGAT','ATCATTGA','CCAATGAA','TTCATTGG','TCAAGGAA','TTCCTTGA',
            'TAAATGAA','TTCATTTA')

all.counts = sapply(all.CBS, str_count, string=cic_seqs)
all.counts = cbind(all.counts, 'Total'=rowSums(all.counts))

hist(all.counts[,'Total'])
length(which(all.counts[,'Total']==0))

