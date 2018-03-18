
source("funcs.r")

# read the genotypes for respective chromosomes
chr2 <- read.table('hmel2.5.30chr2est.txt', header= T)
chr7 <- read.table('hmel2.5.30chr7est.txt', header= T)
chr10 <- read.table('hmel2.5.30chr10est.txt', header= T)
chr18 <- read.table('hmel2.5.30chr18est.txt', header= T)
chr21 <- read.table('hmel2.5.30chr21est.txt', header= T)

chrl <- list(chr2, chr7, chr10, chr18, chr21)

chromnum <- c(2, 7, 10, 18, 21)
chromnames <- paste('chrom', chromnum, sep= '')
names(chrl) <- paste0('chr', chromnum)


# read the individual info for all 30 included inds

ids <- read.csv('kronforst2013ids.txt', sep= '\t', header= T)
ids <- ids[c(1:10, 13:32),]	# throw out h665 and i02_210
ids[] <- lapply(ids, function(x) if(is.factor(x)) factor(x) else x) # drop old factors 


# get allele frequency differences from genotype estimates

pls <- pML(chrl, ids)

pDifls <- list()
for(i in 1:length(pls)){
	pDifls[[i]] <- getAFdiff(pls[[i]])
}

# get the 99th quantile 

pQls <- list()
for(i in 1:length(pDifls)){
	pQls[[i]] <- getPquant(pDifls[[i]], 0.99)
}


# run correlation analyses for each chromosome and species pair (i.e. 15 runs)

cgmr2 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 2, pops= c('cgal', 'mros'))
cgmr7 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 7, pops= c('cgal', 'mros'))
cgmr10 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 10, pops= c('cgal', 'mros'))
cgmr18 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 18, pops= c('cgal', 'mros'))
cgmr21 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 21, pops= c('cgal', 'mros'))

cgmrls <- list(cgmr2, cgmr7, cgmr10, cgmr18, cgmr21)
names(cgmrls) <- chromnames

pamr2 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 2, pops= c('pach', 'mros'))
pamr7 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 7, pops= c('pach', 'mros'))
pamr10 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 10, pops= c('pach', 'mros'))
pamr18 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 18, pops= c('pach', 'mros'))
pamr21 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 21, pops= c('pach', 'mros'))

pamrls <- list(pamr2, pamr7, pamr10, pamr18, pamr21) 
names(pamrls) <- chromnames

cgpa2 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 2, pops= c('cgal', 'pach'))
cgpa7 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 7, pops= c('cgal', 'pach'))
cgpa10 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 10, pops= c('cgal', 'pach'))
cgpa18 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 18, pops= c('cgal', 'pach'))
cgpa21 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 21, pops= c('cgal', 'pach'))

cgpals <- list(cgpa2, cgpa7, cgpa10, cgpa18, cgpa21)
names(cgpals) <- chromnames

Qldls <- list(cgmrls, pamrls, cgpals)
names(Qldls) <- c('cgmr', 'pamr', 'cgpa')



# get stats for distributions of correlation coefficients

library(moments)

datMu <- getQrStats(Qldls[[3]], mean)
write.table(datMu, 'pgStats/afDiffs/cgpaLogMean.txt')
rm(datMu)

datMed <- getQrStats(Qldls[[1]], median)
write.table(datMed, 'pgStats/afDiffs/cgmrLogMed.txt')
rm(datMed)


datSd <- getQrStats(Qldls[[3]], sd)
write.table(datSd, 'pgStats/afDiffs/cgpaLogSd.txt')
rm(datSd)

#datSkew <- getQrStats(Qldls[[3]], skewness)
#write.table(datSkew, 'pgStats/afDiffs/cgpaLogSkew.txt')
#rm(datSkew)

#datKurt <- getQrStats(Qldls[[3]], kurtosis)
#write.table(datKurt, 'pgStats/afDiffs/cgpaLogKurt.txt')
#rm(datKurt)



# plot density curves for correlation coefficients  (figures 6 & 7, and resp. figures in the supplementary materials)

# within chromosomes

plotDensCurve(cgmrls[[1]], 1, fname= 'figs/dc_cgmr2_cg.pdf')
plotDensCurve(cgmrls[[1]], 2, fname= 'figs/dc_cgmr2_mr.pdf')
plotDensCurve(pamrls[[1]], 1, fname= 'figs/dc_pamr2_pa.pdf')
plotDensCurve(pamrls[[1]], 2, fname= 'figs/dc_pamr2_mr.pdf')
plotDensCurve(cgpals[[1]], 1, fname= 'figs/dc_cgpa2_cg.pdf')
plotDensCurve(cgpals[[1]], 2, fname= 'figs/dc_cgpa2_pa.pdf')

plotDensCurve(cgmrls[[2]], 1, fname= 'figs/dc_cgmr7_cg.pdf')
plotDensCurve(cgmrls[[2]], 2, fname= 'figs/dc_cgmr7_mr.pdf')
plotDensCurve(pamrls[[2]], 1, fname= 'figs/dc_pamr7_pa.pdf')
plotDensCurve(pamrls[[2]], 2, fname= 'figs/dc_pamr7_mr.pdf')
plotDensCurve(cgpals[[2]], 1, fname= 'figs/dc_cgpa7_cg.pdf')
plotDensCurve(cgpals[[2]], 2, fname= 'figs/dc_cgpa7_pa.pdf')

plotDensCurve(cgmrls[[3]], 1, fname= 'figs/dc_cgmr10_cg.pdf')
plotDensCurve(cgmrls[[3]], 2, fname= 'figs/dc_cgmr10_mr.pdf')
plotDensCurve(pamrls[[3]], 1, fname= 'figs/dc_pamr10_pa.pdf')
plotDensCurve(pamrls[[3]], 2, fname= 'figs/dc_pamr10_mr.pdf')
plotDensCurve(cgpals[[3]], 1, fname= 'figs/dc_cgpa10_cg.pdf')
plotDensCurve(cgpals[[3]], 2, fname= 'figs/dc_cgpa10_pa.pdf')

plotDensCurve(cgmrls[[4]], 1, fname= 'figs/dc_cgmr18_cg.pdf')
plotDensCurve(cgmrls[[4]], 2, fname= 'figs/dc_cgmr18_mr.pdf')
plotDensCurve(pamrls[[4]], 1, fname= 'figs/dc_pamr18_pa.pdf')
plotDensCurve(pamrls[[4]], 2, fname= 'figs/dc_pamr18_mr.pdf')
plotDensCurve(cgpals[[4]], 1, fname= 'figs/dc_cgpa18_cg.pdf')
plotDensCurve(cgpals[[4]], 2, fname= 'figs/dc_cgpa18_pa.pdf')

plotDensCurve(cgmrls[[5]], 1, fname= 'figs/dc_cgmr21_cg.pdf')
plotDensCurve(cgmrls[[5]], 2, fname= 'figs/dc_cgmr21_mr.pdf')
plotDensCurve(pamrls[[5]], 1, fname= 'figs/dc_pamr21_pa.pdf')
plotDensCurve(pamrls[[5]], 2, fname= 'figs/dc_pamr21_mr.pdf')
plotDensCurve(cgpals[[5]], 1, fname= 'figs/dc_cgpa21_cg.pdf')
plotDensCurve(cgpals[[5]], 2, fname= 'figs/dc_cgpa21_pa.pdf')



# between chromosomes

plotDensCurveOtherChr(cgmrls[[1]], 1, fname= 'figs/dcOth_cgmr2_cg.pdf')
plotDensCurveOtherChr(cgmrls[[1]], 2, fname= 'figs/dcOth_cgmr2_mr.pdf')
plotDensCurveOtherChr(pamrls[[1]], 1, fname= 'figs/dcOth_pamr2_pa.pdf')
plotDensCurveOtherChr(pamrls[[1]], 2, fname= 'figs/dcOth_pamr2_mr.pdf')
plotDensCurveOtherChr(cgpals[[1]], 1, fname= 'figs/dcOth_cgpa2_cg.pdf')
plotDensCurveOtherChr(cgpals[[1]], 2, fname= 'figs/dcOth_cgpa2_pa.pdf')

plotDensCurveOtherChr(cgmrls[[2]], 1, fname= 'figs/dcOth_cgmr7_cg.pdf')
plotDensCurveOtherChr(cgmrls[[2]], 2, fname= 'figs/dcOth_cgmr7_mr.pdf')
plotDensCurveOtherChr(pamrls[[2]], 1, fname= 'figs/dcOth_pamr7_pa.pdf')
plotDensCurveOtherChr(pamrls[[2]], 2, fname= 'figs/dcOth_pamr7_mr.pdf')
plotDensCurveOtherChr(cgpals[[2]], 1, fname= 'figs/dcOth_cgpa7_cg.pdf')
plotDensCurveOtherChr(cgpals[[2]], 2, fname= 'figs/dcOth_cgpa7_pa.pdf')

plotDensCurveOtherChr(cgmrls[[3]], 1, fname= 'figs/dcOth_cgmr10_cg.pdf')
plotDensCurveOtherChr(cgmrls[[3]], 2, fname= 'figs/dcOth_cgmr10_mr.pdf')
plotDensCurveOtherChr(pamrls[[3]], 1, fname= 'figs/dcOth_pamr10_pa.pdf')
plotDensCurveOtherChr(pamrls[[3]], 2, fname= 'figs/dcOth_pamr10_mr.pdf')
plotDensCurveOtherChr(cgpals[[3]], 1, fname= 'figs/dcOth_cgpa10_cg.pdf')
plotDensCurveOtherChr(cgpals[[3]], 2, fname= 'figs/dcOth_cgpa10_pa.pdf')

plotDensCurveOtherChr(cgmrls[[5]], 1, fname= 'figs/dcOth_cgmr21_cg.pdf')
plotDensCurveOtherChr(cgmrls[[5]], 2, fname= 'figs/dcOth_cgmr21_mr.pdf')
plotDensCurveOtherChr(pamrls[[5]], 1, fname= 'figs/dcOth_pamr21_pa.pdf')
plotDensCurveOtherChr(pamrls[[5]], 2, fname= 'figs/dcOth_pamr21_mr.pdf')
plotDensCurveOtherChr(cgpals[[5]], 1, fname= 'figs/dcOth_cgpa21_cg.pdf')
plotDensCurveOtherChr(cgpals[[5]], 2, fname= 'figs/dcOth_cgpa21_pa.pdf')




# Figure 8 (1x2 figure with s and n (different sets which could represent different time points)

selSets <- list()
selSets[[1]] <- list(cgpals[[1]]$cgal$LDmat, cgpals[[1]]$cgal$LDmatN)
selSets[[2]] <- list(pamrls[[1]]$mros$LDmat, pamrls[[1]]$mros$LDmatN)
selSets[[3]] <- list(pamrls[[5]]$mros$LDmat, pamrls[[5]]$mros$LDmatN)
selSets[[4]] <- list(cgmrls[[5]]$mros$LDmat, cgmrls[[5]]$mros$LDmatN)

sVec <- rbind(unlist(lapply(selSets, '[[', 1)))
lsets <- unlist(lapply(lapply(selSets, '[[', 1), length))
sites <- c(rep(1, lsets[1]), rep(2, lsets[2]), rep(3, lsets[3]), rep(4, lsets[4]))
sDF <- data.frame(sites= sites, logR2= sVec)

myplot <- ggplot(sDF, aes(logR2, colour=sites)) + coord_cartesian(xlim= c(-8, 0), ylim= c(0, 0.4)) + geom_density(alpha= 0.2, show.legend= F, adjust= 1, size= 1.1) + theme_minimal() #+ scale_color_manual(values= c(x= cols[1], y= cols[2], z= cols[3]))
ggsave(myplot, filename= 'dc_sSets.pdf')

##
nVec  <- rbind(unlist(lapply(selSets, '[[', 2)))
nDF <- data.frame(sites= sites, logR2= sVec)

myplot <- ggplot(nDF, aes(logR2, colour=sites)) + coord_cartesian(xlim= c(-8, 0), ylim= c(0, 0.4)) + geom_density(alpha= 0.2, show.legend= F, adjust= 1, size= 1.1) + theme_minimal() #+ scale_color_manual(values= c(x= cols[1], y= cols[2], z= cols[3]))
ggsave(myplot, filename= 'dc_nSets.pdf', width= 7, height= 7)



#####################################

# Fst correlation coefficients (supplement)




