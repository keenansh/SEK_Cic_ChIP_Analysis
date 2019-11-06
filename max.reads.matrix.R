max.reads.matrix<- function(condition,peaks)
{
  
  fin.reads=c(1:length(consensus_peaks))
  
  for (i in 1:length(peaks)){
    f=findOverlaps(peaks[i],condition)
    s=subjectHits(f)
    s2=(condition[s]$score)
    w=which.max(s2)
    fin.reads[i]=s2[w]
  }
  
  return(fin.reads)
  
}

