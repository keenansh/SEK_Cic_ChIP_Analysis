# these are subfunctions that are useful for the analysis of PWMs

get.all.PWM.instances = function(PWM, pct.match, genome, chrlist)
{
	# this function takes as input the PWM of interest, the percent match limit (in the format "80%", including the quotes), the BSgenome object of your genome of interest, and a list of chromosomes you want to query (or 'all'). It outputs a genomic ranges object that contains the location of each instance of the PWM.
	motif = PWM
	motif.gr = GRangesList()
	cat('extracting motif locations from the genome of interest. This triggers warnings, they can be ignored.\n')
	chrlist = which(seqnames(genome)%in%chrlist)
	
	for (i in 1:length(chrlist))
	{
		watson = matchPWM(motif,genome[[chrlist[i]]],min.score = pct.match)
		gr = GRanges(seqnames = rep(seqnames(genome)[chrlist[i]],length(watson)),IRanges(start(watson),end(watson)),strand = rep('+',length(watson)))
		crick = matchPWM(reverseComplement(motif),genome[[chrlist[i]]],min.score = pct.match)
		gr = append(gr,GRanges(seqnames = rep(seqnames(genome)[chrlist[i]],length(crick)),IRanges(start(crick),end(crick)),strand = rep('-',length(crick))))
		motif.gr[[i]] = gr
	}
	motif.gr = unlist(motif.gr)
	seqlevels(motif.gr) = seqlevels(genome)
	seqinfo(motif.gr) = seqinfo(genome)
	seqlevels(motif.gr) = seqlevelsInUse(motif.gr)
	genome(motif.gr) = unique(genome(genome))
	motif.gr
}

basecount = function(windows,data,chrlist = c('chr2L','chr2R','chr3L','chr3R','chr4','chrX'))
{
	# this is an instance of the general basecount function that takes as input a set of bins (in GRanges format), a GRanges object containing count data, and a list of chromosomes of interest. It returns a GRanges object identical to the input bins, with a score column indicating the average number of instances where a feature in the count data was present within the bin. It does not normalize to anything.
	require(GenomicRanges)
	require(GenomicAlignments)
	binnedAverage = function(bins,numvar)
	{
		stopifnot(is(bins,"GRanges"))
		stopifnot(is(numvar,"RleList"))
		stopifnot(identical(seqlevels(bins),names(numvar)))
		bins_per_chrom = split(ranges(bins),seqnames(bins))
		means_list = lapply(names(numvar),
			function(seqname){
				views = Views(numvar[[seqname]],
				bins_per_chrom[[seqname]])
			viewMeans(views)	
			})
		new_mcol = unsplit(means_list,as.factor(seqnames(bins)))
		mcols(bins)[['basecount']] = new_mcol
		bins
	}
		
	good.uns = chrlist
	windows = windows[seqnames(windows)%in%good.uns]
	seqlevels(windows) = seqlevelsInUse(windows)
	cov = coverage(data)
	cov = cov[match(seqlevels(windows),names(cov))]
	
	bcount = binnedAverage(windows,cov)
	bcount = mcols(bcount)[,1]
	score(windows) = bcount
	return(windows)
	
}

count.PWM.in.ROI = function(ROI.gr,motif.basecount,flank,exclude.empty)
{
	# this function takes as input a GRanges of peaks (ROI.gr), the output of the basecount function (motif.basecount), the size in basepairs to flank the region of interest, and a logical that specifies whether the exclude regions of interest that are empty or do not have any overlap with a motif. It returns a matrix of counts where rows are basepairs and columns are peaks.
	
	require(GenomicRanges)
	siglist = ROI.gr
	motif.bc = motif.basecount
	if(missing(exclude.empty)){exclude.empty = FALSE}
	if(!all(width(siglist)==1))
	{
		cat(paste0('\nRegions of interest will be resized to the ',as.numeric(flank),'bp flanking the center of the region.\n'))
		siglist = resize(siglist,width = 1, fix = 'center')
	}
	
	sigindex = subjectHits(findOverlaps(siglist,motif.bc))
	sigmat = sapply(c(1:length(sigindex)),function(x) seq(sigindex[x]-flank,sigindex[x]+flank,1))
	countmat = sapply(c(1:ncol(sigmat)),function(x) score(motif.bc[sigmat[,x]]))
	countmat[,which(strand(siglist)=='-')] = rev(countmat[,which(strand(siglist)=='-')])
	
	if(exclude.empty){countmat = countmat[,!colSums(countmat)<6];cat('\nexcluding peaks with fewer than 6bp overlap with a motif.\n')}
	countmat
}

count.ATAC.feature.around.motif = function(ROI.gr,ATAC.basecount,flank,motif.GRanges,exclude.chr,orient)
{
	siglist = ROI.gr
	motif.gr = motif.GRanges
	nucs = ATAC.basecount
	nucflank = round(flank/width(nucs[1]))
	if(missing(orient)){orient = FALSE}
	
	if(all(width(siglist)==1)){
		#if the input Regions of Interest were supplied as single base pair ranges, this resizes them to represent the broader regions of interest.
		broad.siglist = siglist
		start(broad.siglist) = start(broad.siglist)-flank
		end(broad.siglist) = end(broad.siglist)+flank
	}else{
		#if the input Regions of Interest were supplied as not single base pair ranges, they are resized to a uniform size based on the flanking basepair input variable, centered around their original centers.
		broad.siglist = resize(siglist,width = 1, fix = 'center')
		start(broad.siglist) = start(broad.siglist)-flank
		end(broad.siglist) = end(broad.siglist)+flank
	}
	motif.center = resize(motif.gr,width = 1, fix = 'center')
	motif.center = motif.center[!seqnames(motif.center)%in%exclude.chr]
	motif.center = motif.center[unique(queryHits(findOverlaps(motif.center,broad.siglist)))]
	
	nuc.scores = score(nucs)
	nucindex = subjectHits(findOverlaps(motif.center,nucs))
	nucmat = sapply(c(1:length(nucindex)),function(x) seq(nucindex[x]-nucflank,nucindex[x]+nucflank,1))
	nuccount = sapply(c(1:ncol(nucmat)),function(x) nuc.scores[nucmat[,x]])
	if(orient){
	nuccount[,which(strand(motif.center) == '-')] = rev(nuccount[,which(strand(motif.center) == '-')])
	}
	nuccount
}

count.ATAC.feature.within.peaks = function(ROI.gr,ATAC.basecount,flank)
{
	siglist = ROI.gr
	nucs = ATAC.basecount
	nucflank = round(flank/width(nucs[1]))
	
	siglist = resize(siglist,width = 1, fix = "center")
	
	nuc.scores = score(nucs)
	nucindex = subjectHits(findOverlaps(siglist,nucs))
	nucmat = sapply(c(1:length(nucindex)),function(x) seq(nucindex[x]-nucflank,nucindex[x]+nucflank,1))
	nuccount = sapply(c(1:ncol(nucmat)),function(x) nuc.scores[nucmat[,x]])
	nuccount
}

get.edgeR.norm.factors = function(list.of.GRanges,ROI)
{
	require(edgeR)
	dge = NULL
	for(i in 1 : length(list.of.GRanges))
	{
		x = countOverlaps(ROI,list.of.GRanges[[i]])
		dge = cbind(dge,x)
	}
	colnames(dge) = seq(1,ncol(dge),1)
	cds = DGEList(dge)
	return(cds$samples$lib.size)
}
