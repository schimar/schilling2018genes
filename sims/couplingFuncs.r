# functions to calculate recombination and coupling statistics 

library(rhdf5)
#library(scales)
#library(ape)
#library(spatstat)




calcRecomb <- function(nLoci, map, nChrom= 4,...){
	# calculate recombination rate from numberLoci, map, and number of chromosomes
	if(nLoci <= nChrom){
		recomb <- 0.5 
		return(recomb)}
	else {
		nNeighbors <- nLoci - 1
		nNeighborsSame <- nLoci - nChrom
        nNeighborsDiff <- nNeighbors - nNeighborsSame
        sameChromAvgDist <- map/nLoci
        diffChromDist <- 0.5
        recomb <- ((nNeighborsSame * sameChromAvgDist) + (nNeighborsDiff * diffChromDist)) / nNeighbors
		recomb[which(recomb < 0 )] <- NaN
		recomb[which(recomb > 0.5)] <- 0.5
		return(recomb)
	}
}

singleLocusEq <- function(s, m, singleS= T){
	# function to calculate the equilibrium frequencies for (vectors of) single loci  (with singleS = T for vector of s values, and F when using one single s value (e.g. sMax))
	peq <- (m * (0.5 + 0.75* s) - 0.125* s - (0.125* (sqrt(s^2 - (4* m* s^2) + (4* m^2) * ((2 + s)^2))))) / (-1* (0.25 + m)*s )
	clineWidth <- (peq - (1 - peq))
	if (singleS == T) {
		return(unlist(list(peq, abs(clineWidth))))
	}
	else {
		out <- list()
		out[[1]] <- peq
		out[[2]] <- abs(clineWidth)
		return(out)
	}
}

calcLe <- function(m, pBar, s) {
	# calculate effective number of loci (L_e)
	sStar <- (m * (pBar - 0.5)) / ((0.25 - 0.25* pBar)* pBar + m * (0.5 + pBar * (pBar - 1.5)))
	Le <- sStar / s
	outp <- list(sStar, Le, s, pBar)
	names(outp) <- c('sStar', 'Le', 's', 'pBar')
	return(outp)
}



calcMaxEffMig <- function(s, m, L) {
	# calculate maximum effective migration for a given s vector, mutation rate and number of (selected) loci in one generation
	# NOTE:  only for single population
	
	# calculate equilibrium frequencies with given s   (here: sBar or meanS) 
	peq <- singleLocusEq(s, m, singleS= F)[[1]]
	q <- 1 - peq

	# calc fitnesses and maxEffMig
	randomResFit <- (1 + (peq * s))^L
	randomImmFit <- (1 + (q * s))^L
	patchAvgFit <- (randomImmFit * m) + ((1 - m) * randomResFit)
	#
	maxEffMig <- (randomImmFit * m) / patchAvgFit
	return(list(peq, maxEffMig))
}


calcGWCtime <- function(em, maxEffMig, endAllopatry){
# calculate the time of genome-wide congealingi, i.e. (last time step where effective migration was > maxEffMig ) + 1
	nGen <- unique(em[,1]) 
	emVec <- em[,2]
	t <- length(nGen)
	if (emVec[t] >= 0.5 * maxEffMig[t]) {
		# the effective migration rate can still be considered close to random, thus gwcTime not reached
		gwcTime <- NA
		gwcReached <- 0
	}
	else {
		gwcIndex <- tail(which(emVec >= maxEffMig), n= 1) + 1 
		if (length(gwcIndex) < 1) {
			gwcTime <- endAllopatry
		}
		else {
			gwcTime <- nGen[gwcIndex]
		}
		gwcReached <- 1
	}
	out <- list(gwcTime, gwcReached)
	names(out) <- c('gwcTime', 'gwcReached')
	return(out)
}


geomean = function(x, na.rm=TRUE){
	# function to calculate the geometric mean of vector x
	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



calcMeanS <- function(sSel, gen, mutDist, sConst) { 
	# calculate mean s for each generation in a given run
	if (mutDist == 3) {      # constant s
		meanS <- rep(sConst, gen)
		return(meanS)
	}
	else {
		meanS <- lapply(sSel, mean)		# geomean
		return(unlist(meanS))
	}
}



readCCobj <- function(run, path = '.', ...) {
	# function to read a single data set (as input for ccStats)
	path5 <- paste('/runs/', run, sep= '')
	fstTmp <- as.data.frame(h5read(paste(path, '/fst_AFt.h5', sep= ''), name= path5)[[1]])
	H5close()
	aftsTmp <- as.data.frame(h5read(paste(path, '/afts_AFt.h5', sep= ''), name= path5)[[1]])
	colnames(fstTmp) <- c("nGen", "locusID", "Fst", "allele_frequencies", "S_MAX1", "S_MAX0", "chromosomeMembership", "MAP", "locType")
	colnames(aftsTmp) <- c("nGen", "locusID", "AFpatch0", "AFpatch1", "is_reversed_locus", "locType", "AF", "AFdiff")
	H5close()
	#
	LDselTmp <- h5read(paste(path, '/LDselAvg_AFt.h5', sep= ''), name= path5)[[1]]
	H5close()
	LDneuTmp <- h5read(paste(path, '/LDneutAvg_AFt.h5', sep= ''), name= path5)[[1]]
	H5close()
	effMig <- h5read(paste(path, '/effMig_AFt.h5', sep= ''), name= path5)[[1]]
	colnames(effMig) <- c("nGen", "eme0", "eme1", "nVariableLoci", "nRes", "nImm", paste('V', seq(7, 26,1)))
	H5close()
	#dXY <- h5read(paste(path, '/dXY_AFt.h5', sep= ''), name= path5)[[1]]
	#colnames(dXY) <- c('nGen', 'dXY', 'deme0', 'deme1')
	H5close()
	#
	out <- list(fstTmp, aftsTmp, LDselTmp, LDneuTmp, effMig)	#, dXY)
	names(out) <- c('fst', 'afts', 'LDsel', 'LDneut', 'effMig')	#, 'dXY')
	return(out)
}



ccStats.3 <- function(run, df, ccObj, maf= 25e-4, nChrom= 4) {    #fst, afts, LDsel, LDneut, effMig, run, maf= 25e-4, nChrom= 4) {
	# function to calculate coupling/congealing stats for a given run 
	fst <- ccObj$fst
	afts <- ccObj$afts
	LDsel <- ccObj$LDsel #emRed <- em[em$nGen %in% unique(fst$nGen),]
	#LDsel <- as.data.frame(ccObj$LDsel)[ccObj$LDsel[,1] %in% names(ccObj$phiObs),]
	LDneut <- ccObj$LDneut
	effMig <- ccObj$effMig[,1:2]
	params <- df[which(df$run == run),]
	m = params$sd_move
	fstSpl <- split(fst, fst$nGen)
	aftsSpl <- split(afts, afts$nGen)
	nLoci <- unlist(lapply(lapply(fstSpl, dim), '[', 1))
	gen <- length(aftsSpl)
	fstSplS <- lapply(fstSpl, subset, locType == 1)
	aftsSplS <- lapply(aftsSpl, subset, locType == 1)
	nLociS <- as.numeric(unlist(lapply(fstSplS, nrow)))
	p <- lapply(fstSplS, '[[', 4)
	sSel <- lapply(fstSplS, '[[', 6)
	# 
	sConst <- params$deme0_constant_s
	mutDist <- params$mutation_distribution
	#
	meanS <- calcMeanS(sSel, mutDist, gen, sConst)

	# calc Kruuk's phi and phiObs 
	PHIs <- calcPHIs.3(aftsSpl, fstSpl, maf= maf, mapL= params$total_map_length, nChrom= params$nchromosomes, meanS= meanS)     # NOTE: we use a MAF threshold here! 
	#
	nLoci <- as.data.frame(cbind(nLoci, nLociS, PHIs$nLoci, PHIs$nLociS))
	colnames(nLoci) <- c('nLoci', 'nLociS', 'nLocimaf', 'nLociSmaf')
	recomb <- PHIs$recomb
	#  afDiffs
	afDiff_s <- lapply(aftsSplS, '[[', 8)
	afDiff_n <- lapply(lapply(aftsSpl, subset, locType == 0), '[[', 8)
	avgAFdiff <- unlist(lapply(lapply(lapply(aftsSpl, '[[', 8), abs), mean))    # mean(abs(allAFdiffsPerGen))
	#
	# calc single locus exp. under full coupling (sMax) 
	sMaxlist <- lapply(PHIs$sMax, singleLocusEq, m= m, singleS= T)
	pHatsMax <- unlist(lapply(sMaxlist, '[[', 1))
	cWsMax <- unlist(lapply(sMaxlist, '[[', 2))
	# calc single locus exp. using actual selection coefficients (sSel == S_MAX0)
	sBarlist <- lapply(lapply(PHIs$sBar, singleLocusEq, m= m, singleS= F), unlist)
	pHatsBar <- unlist(lapply(sBarlist, '[[', 1))
	cWsBar <- unlist(lapply(sBarlist, '[[', 2))
	
	# equilibrium freq & cline width for all sSel
	clineWidthAllS <- lapply(sSel, singleLocusEq, m= m, singleS= F)
	pBarAllS <- lapply(lapply(clineWidthAllS, '[[', 1), unlist)
	pBar <- unlist(lapply(pBarAllS, mean))
	clineWallS <- lapply(lapply(clineWidthAllS, '[[', 2), unlist)

	# get Fst values (s, n & total) per generation
	fstSel <- unlist(lapply(lapply(lapply(fstSpl, function(x)x[x$locType == 1,]), '[[', 3), mean))
	fstNeut <- unlist(lapply(lapply(lapply(fstSpl, function(x)x[x$locType == 0,]), '[[', 3), mean))
	fstAll <- unlist(lapply(lapply(fstSpl, '[[', 3), mean))
	FSTs <- as.data.frame(cbind(fstSel, fstNeut, fstAll))
	names(FSTs) <- c('FSTsel', 'FSTneut', 'FSTtot')
	# get effective s and Le   
	sSLe <- calcLe(m, pBar, meanS)
	sStarLeS <- as.data.frame(do.call(cbind, sSLe))
	names(sStarLeS) <- c('sStar', 'Le', 's', 'pBar')
		# 
	LDsell <- LDsel[, c(1,4)]
	LDneutl <- LDneut[, c(1,4)]
	#FSTout <- as.data.frame(matrix(unlist(FSTs),ncol=3, byrow=TRUE, dimnames= list(NULL, c('FSTneut', 'FSTsel', 'FSTtot'))))
	# maxEffMig with meanS (mean of sMax[i])
	maxEffMigMeanS <- calcMaxEffMig(meanS, m, unlist(lapply(fstSpl, length)))[[2]]
	gwcTimeMeanS <- calcGWCtime(effMig, maxEffMigMeanS, params$end_period_allopatry)
##### output
	out <- list(FSTs, LDsell, LDneutl, afDiff_s, afDiff_n, avgAFdiff, pHatsBar, cWsBar, meanS, sStarLeS, m, effMig, unlist(maxEffMigMeanS), gwcTimeMeanS, clineWallS, pBarAllS, nLoci, maf, recomb, PHIs$kphiMeanS, PHIs$thetaMeanS)
	names(out) <- c('FSTs', 'LDsel', 'LDneut', 'afDiffS', 'afDiffN', 'avgAFdiffs', 'pHatsBar', 'cWsBar', 'meanS', 'sStarLeS', 'sd_move', 'effMig', 'maxEffMigMeanS', 'gwcTimeMeanS', 'cWallS', 'pBarAllS', 'nLoci', 'maf', 'recomb', 'kphiMeanS', 'thetaMeanS')
	return(out)
}



plotStatic.2 <- function(ccObj, run, data, i= '', ...) {
	# function to create the per-run plots with LD, Fst, PHIs, Le and me
	close.screen(all.screens= T)
	par(oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#par(bg = "white") # erase.screen() will appear not to work if the background color is transparent 
	split.screen(c(2,4))  
	#mtext(run, outer = TRUE )
	mainList <- c("neutral sites", "selected sites")
	LDlist <- list(ccObj$LDneut, ccObj$LDsel)
	param <- data[which(data$run == run),]
	m <- param$sd_move
	s <- param$mean_s
	mutdist <- param$mutation_distribution  #[which(data$run == run)]
	tsFreq <- param$ts_sampling_frequency
	xLab <- paste('sampling times (tsFreq =', tsFreq, ')')
	#
	# LD
	screen(1)
	plot(1:length(ccObj$phiObs), LDlist[[1]][,2], xlim= c(0, length(ccObj$phiObs)+2), ylim= c(0, 1), main= 'LD', ylab= 'Avg LD', xlab= xLab, type= 'l', col= 'deepskyblue1')
	points(1:length(ccObj$phiObs), LDlist[[2]][,2], col= 'firebrick1', type= 'l')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('blue', 'red'))
	#
	# Fst 
	screen(2)
	cols <- c('deepskyblue1', 'firebrick1', 'black') #t(rainbow(3))
	plot(ccObj$FSTs$FSTneut, type= 'l', col= cols[1], ylim= c(0,1), xlim= c(0, length(ccObj$phiObs)+2), ylab= expression('avg F'[st]), main= expression('avg F'[st]), xlab= xLab)
	points(ccObj$FSTs$FSTsel, col= cols[2], type= 'l')
	points(ccObj$FSTs$FSTtot, col= cols[3], type= 'l', lty= 1)
	legend('bottomright', legend= c('neutral', 'selected', 'total'), fill= cols, cex= 0.8)
	#
	# PHIs per generation
	screen(3) #, new= F)
	plot(log10(ccObj$phiObs), ylab= expression(paste('log'[10]~ phi)), xlab= xLab, type= 'l', ylim= c(-2, max(log10(ccObj$kphisMax))+ 1), col= 'grey70')  # xlim=  c(0, length(ccObj$phiObs)+100)
	points(log10(ccObj$kphisMax), type= 'l', col= 'black')
	legend('topleft', legend= c(expression(paste(phi[Kruuk])), expression(paste(phi, ' ', bar(s)))), fill= c('black', 'grey70'), cex= 0.8)
	text(x= 0.8*length(ccObj$meanS), y= 0.6*max(log10(ccObj$kphisMax)), paste('s = ', s))
	text(x= 0.8*length(ccObj$meanS), y= 0.5*max(log10(ccObj$kphisMax)), paste('m = ', m))

	#
	# cline widths and PHIs
	screen(4)
	cWallS <- lapply(ccObj$cWallS, unlist)
	phiOncW <- mapply(rep, ccObj$phiObs, times= unlist(lapply(cWallS, length)))
	#yLim <- max(unlist(ccObj$clineWidthSmax))
	plot(log10(ccObj$kphisMax), ccObj$pHatsMax, type= 'l', ylim= c(-0.45,  1), xlim= c(-2, max(log10(ccObj$kphisMax))), xlab= expression(paste('log'[10], ' ', phi)), ylab= expression(paste('p'[i0]~'- p'[i1])), col= 'black')
	#abline(h= 0, lty= 3)
	#}
	points(log10(unlist(phiOncW)), unlist(cWallS), pch= '.', col= 'grey70') #type= 'l')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', cex= 1.4, col= 'firebrick1')
	points(log10(unlist(ccObj$phiObs)), unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', cex= 1.4, col= 'deepskyblue1')
	#points(log10(unlist(ccObj$phiObs)), unlist(ccObj$avgAFDiff), type= 'l', col= 'black')

	legend('bottomright', legend= c(expression(paste(phi[Kruuk], ' ~ ', 'peq '[sMax])), expression(paste(phi, ' ~ ', bar(p), ' all s')), expression(paste(phi, ' ~ avg p S')), expression(paste(phi, ' ~ avg p N'))), fill= c('black', 'grey70', 'firebrick1', 'deepskyblue1'), cex= 0.8)
	#
	# plot Le and me
	screen(5)
	#plot(1, type="n", xlab= 'gen / 1e3', ylab="", ylim=c(-0.1, 0.1), xlim=c(0, (ccObj$sStarLeS)))
	emig <- ccObj$effMig[,2]
	emig[which(emig == 0)] <- 1e-10
	plot(log10(abs(ccObj$sStarLeS$Le)), type= 'l', xlab= xLab, ylab= expression(paste('log'[10], '.')) , ylim= c(-10, 7))#, ylim= c(min(log10(emig)), max(log10(abs(ccObj$sStarLeS$Le)))))      # ylab was expression(paste('log'[10]~'L'[e]))
	abline(v= which.min(log10(abs(ccObj$sStarLeS$Le))), lty= 2, col= 'black')
	points(log10(emig), col= 'deepskyblue1', type= 'l')
	points(log10(ccObj$sBar), type= 'l', col= 'orange')
	points(log10(ccObj$sMax), type= 'l', col= 'firebrick1')
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, lty= 2, col= 'green')
	legend('bottomright', legend= c(expression('L'[e]), expression('m'[e]), 'gwcTime meanS', expression(bar(s)), 'sMax'), lty= c(1, 1, 2), col= c('black', 'deepskyblue1', 'green', 'orange', 'firebrick1'), cex= 0.8)
	text(x= 0.6*length(ccObj$meanS), y= 0.6*max(log10(emig)), paste("gwc time = ", ccObj$gwcTimeMeanS$gwcTime))
	#
	# plot effMig  
	screen(6)
	plot(ccObj$maxEffMigSbar, col= 'orange', type= 'l', ylim= c(0, 0.25), ylab= 'migration rate', xlab= xLab)
	points(ccObj$effMig[,2], col= 'black', type= 'l')
	points(ccObj$maxEffMigMeanS, col= 'green', type= 'l')
	points(ccObj$recomb, type= 'l', col= 'cyan')

	abline(h= m, col= 'firebrick1', lty= 3)
	abline(h= s, col= 'deepskyblue1', lty= 3)
	abline(v= ccObj$gwcTimeSbar$gwcTime/tsFreq, col= 'orange', lty= 2)
	abline(v= ccObj$gwcTimeMeanS$gwcTime/tsFreq, col= 'green', lty= 2)
	endAllo <- data$end_period_allopatry[which(data$run == run)]
	if (ccObj$gwcTimeMeanS$gwcTime == endAllo) {
		text(x= 0.3* length(ccObj$meanS), y= 1.1*max(ccObj$effMig[,2]), expression(paste('gwcTime == endPallo')))
	}
	text(x= 0.8* length(ccObj$meanS), y= 0.4*max(ccObj$effMig[,2]), paste('mutdist = ', mutdist)) 
	text(x= 0.3* length(ccObj$meanS), y= 1.3*max(ccObj$effMig[,2]), paste('endPallo = ', param$end_period_allopatry))
	cols <- c('orange', 'green', 'black', 'cyan', 'firebrick1', 'deepskyblue1', 'orange', 'green')
	legend('topright', legend= c(expression(paste('maxEffMig ', bar(s))), 'maxEffMig meanS', expression(paste('m'[e0])), 'recomb.', 'm', 's', expression(paste('gwcTime ', bar(s))), 'gwcTime meanS'), col= cols, lty= c(1,1,1,3,3,2,2), pch= c(NA,NA,NA,NA,NA,NA,NA), pt.cex= 2, seg.len= 0.9, cex= 0.7)
	#
	screen(7)
	maflgnd <- paste('MAF (', ccObj$maf, ')', sep= '')
	plot(ccObj$nLoci$nLoci, type= 'l', ylab= 'numbers of loci', ylim= c(0, 1.1*max(ccObj$nLoci$nLoci)))
	points(ccObj$nLoci$nLociS, type= 'l', col= 'red')
	points(ccObj$nLoci$nLocimaf, type= 'l', col= 'blue')
	points(ccObj$nLoci$nLociSmaf, type= 'l', col= 'yellow')
	legend('topleft', legend= c('total', paste('total', maflgnd, sep= ' '), 'selected', paste('selected ',maflgnd, sep= ' ')), fill= c('black', 'blue', 'red', 'yellow'), cex= 0.8)

	#
	screen(8)
	plot(unlist(lapply(lapply(ccObj$afDiffS, abs), mean)), pch= '.', col= 'firebrick1', ylab= 'allele freq diffs', cex= 1.4)
	points(unlist(lapply(lapply(ccObj$afDiffN, abs), mean)), pch= '.', col= 'deepskyblue1')
	legend('bottomright', legend= c('neutral', 'selected'), fill= c('deepskyblue1', 'firebrick1'))
	#
	close.screen(all= T)
	title(paste(run, ' (',i, ')', sep= ''), outer= T)
}


wrapH5static <- function(data, setname, path= '/media/schimar/dapperdata/bu2s/h5/', maf= 25e-4, sleep= 0,...) {
	# function to read individual runs (from vector of runs), calculate CC and plotStatic
	#
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		ccTmp <- ccStats.3(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		#plotStatic.1(ccTmp, run, data)
		plotStatic.2(ccTmp, run, data, i= i)
		#
		Sys.sleep(sleep)
		H5close()
	}
}


calcPHIs.3 <- function(aftsSpl, fstSpl, maf= 25e-4, mapL= 100, nChrom= 4, meanS= NULL) {
	#afts <- subset(afts, fst$allele_frequencies >= maf)
	#fst <- subset(fst, fst$allele_frequencies >= maf)
	aftsSplmaf <- lapply(aftsSpl, function(x)x[x$AF >= 0.1,])
	fstSplmaf <- lapply(fstSpl, function(x)x[x$allele_frequencies >= 0.1,])
	nLoci <- unlist(lapply(lapply(fstSplmaf, dim), '[', 1))
	# 
	fstSplS <- lapply(fstSplmaf, subset, locType == 1)
	aftsSplS <- lapply(aftsSplmaf, subset, locType == 1)
	nLociS <- as.numeric(unlist(lapply(fstSplS, nrow)))
	p <- lapply(fstSplS, '[[', 4)
	sSel <- lapply(fstSplS, '[[', 6)
	# calc recombination 
	recomb <- unlist(lapply(nLociS, calcRecomb, map= 0.02*mapL, nChrom= nChrom))
	# calc sMax
	sMax <- unlist(lapply(sSel, function(x) prod(1 + x)-1))
	# calc Kruuk's phi 
	kphiSmax <- ((unlist(nLociS) - 1) * (sMax / recomb))
	#
	kphiMeanS <- ((unlist(nLociS) -1) * (meanS / recomb))

	# Barton's theta w/ meanS
	thetaMeanS <- meanS / recomb
	# Barton's phi w/ sMax
	phiBarsMax <- sMax / recomb

	#
	out <- list(sMax, kphiMeanS, kphiSmax, nLociS, nLoci, recomb, thetaMeanS, phiBarsMax)
	names(out) <- c('sMax', 'kphiMeans', 'kphiSmax', 'nLociS', 'nLoci', 'recomb', 'thetaMeanS', 'phiBarsMax')
	return(out)
}


xtractPhis <- function(data, folder= '.', path= '.', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains phiObs and kphismax
	#
	runs <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, path)
		ccTmp <- ccStats.3(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		avgAFdiffS <- unlist(lapply(lapply(ccTmp$afDiffS, abs), mean))
		avgAFdiffN <- unlist(lapply(lapply(ccTmp$afDiffN, abs), mean))
		cWallS <- lapply(ccTmp$cWallS, unlist)
		runs[[i]] <- list(avgAFdiffS, avgAFdiffN, cWallS, ccTmp$pBarAllS, ccTmp$sStarLeS$Le, ccTmp$kphiMeanS, ccTmp$thetaMeanS) 
		names(runs)[i] <- run
		names(runs[[i]]) <- c('afDiffS', 'afDiffN', 'cWallS', 'pBarAllS', 'Le', 'kphiMeanS', 'thetaMeanS')
		#phiObs[[i]] <- ccTmp$phiObs
		#names(phiObs)[i] <- run
		#kphisMax[[i]] <- ccTmp$kphisMax
		#names(kphisMax)[i] <- run
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		


xtractLD <- function(data, setname, folder, path= '/media/schimar/FLAXMAN/h5/', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains effMig, Le and gwcTime  
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, setname, folder, path)
		#ccTmp <- ccStats.2(data, ccObjTmp$fst, ccObjTmp$afts, ccObjTmp$LDsel, ccObjTmp$LDneut, ccObjTmp$effMig, run, maf= maf)
		#ccTmp <- ccStats.2(run= run, df= df, ccObj= ccObjTmp, maf= maf)
		#
		runs[[i]] <- list(ccObjTmp$LDneut, ccObjTmp$LDsel)	#, ccObjTmp$dXY)   

		names(runs)[i] <- run
		names(runs[[i]]) <- c('LDneut', 'LDsel')	#, 'dXY')
	}
	#out <- list(phiObs, kphisMax)
	#names(out) <- c('phiObs', 'kphisMax')
	return(runs)
}		

xtractFst <- function(data, folder= '.', path= '.', maf= 25e-4, ...) {
	# function to read individual runs (from vector of runs), calculate CC and create new list (of length(data)) that contains phiObs and kphismax
	#
	runs <- list()
	#kphisMax <- list()
	for (i in 1:dim(data)[1]){
		run <- data$run[i]
		path5 <- paste('/runs/', run, sep= '')
		#
		ccObjTmp <- readCCobj(run, path)
		ccTmp <- ccStats.3(run= run, df= df, ccObj= ccObjTmp, maf= maf)

		runs[[i]] <- ccTmp		 
		names(runs)[i] <- run

	}
	return(runs)
}		





