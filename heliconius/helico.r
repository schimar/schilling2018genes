
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

pamr2 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 2, pops= c('pach', 'mros'))
pamr7 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 7, pops= c('pach', 'mros'))
pamr10 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 10, pops= c('pach', 'mros'))
pamr18 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 18, pops= c('pach', 'mros'))
pamr21 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 21, pops= c('pach', 'mros'))

cgpa2 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 2, pops= c('cgal', 'pach'))
cgpa7 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 7, pops= c('cgal', 'pach'))
cgpa10 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 10, pops= c('cgal', 'pach'))
cgpa18 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 18, pops= c('cgal', 'pach'))
cgpa21 <- getQr_allChrSub(chrls= chrl, pQls= pQls, tChr= 21, pops= c('cgal', 'pach'))



