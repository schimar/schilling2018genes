


####  
source("couplingFuncs.r")

df <- read.table("~/flaxmans/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
names(df) <- tolower(names(df))
#

cols <- c("#c7485a", "#237fa9", "#4a921e")

################################
# sets for manuscript (ms) (in order in figure - ABC, DEF in par(mfrow= c(2,3)))

sA <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.01 & df$mean_s == 0.005),]
sB <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.01 & df$mean_s == 0.01),]
sC <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.01 & df$mean_s == 0.02),]	# tsFreq = 172
#
sD <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.1 & df$mean_s == 0.005),]
sE <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.1 & df$mean_s == 0.01),]
sF <- df[which(df$mutation_distribution == 0 & df$sd_move == 0.1 & df$mean_s == 0.02),]		# tsFreq = 500
#


sAsym <- sA[which(sA$end_period_allopatry == -1),]
sBsym <- sB[which(sB$end_period_allopatry == -1),]
sCsym <- sC[which(sC$end_period_allopatry == -1),]

sDsym <- sD[which(sD$end_period_allopatry == -1),]
sEsym <- sE[which(sE$end_period_allopatry == -1),]
sFsym <- sF[which(sF$end_period_allopatry == -1),]

#symls <- list(sAsym, sBsym, sCsym, sDsym, sEsym, sFsym)


############
# get Fst for set A (figure 1 B)


sAfst <- xtractFst(sAsym[1,])
fobj <- sAfst[[1]]

sels <- list()
neuts <- list()

for(i in 1:length(fobj)){
	sels[[i]] <- fobj[[1]]$FSTs$FSTsel
	neuts[[i]] <- fobj[[1]]$FSTs$FSTneut
}

source('as.data.frame.list.R')
DFs <- as.data.frame(sels)
DFs <- t(DFs)

DFn <- as.data.frame(neuts)
DFn <- t(DFn)

sMed <- apply(DFs, 1, median, na.rm= T)
nMed <- apply(DFn, 1, median, na.rm= T)

#
## figure 1 B

par(mfrow= c(1,2))
xLab <- expression(paste("generations (x 10" ^"3",")"))

plot(1:length(sMed), type= 'n', cex.lab= 1.8, xlab= xLab, ylab= expression(F[ST]), xlim= c(0, 620), ylim= c(0,1), axes= F)
#
axis(1, at= c(0, 100, 200, 300, 400, 500, 600), cex.axis= 1.4)
axis(2, at= c(0, 0.5, 1.0), cex.axis= 1.4)


lines(sMed, lty= 1, cex= 1.2, col= cols[1], lwd= 2.0)
lines(nMed, lty= 1, cex= 1.2, col= cols[2], lwd= 2.0)

box()

#
## figure 1 C

plot(log10(sAfst[[1]]$thetaMeanS), type= 'l', ylab= expression(paste('log'[10], ' metric')), cex.lab= 1.8, cex.axis= 1.4, lwd= 2, xlab= xLab, ylim= c(-2, 5))
lines(log10(sAfst[[1]]$sStarLeS$Le), lwd= 2, lty= 2)

legend('topleft', legend= c(expression(paste(phi)), expression('L'[e])), lty= c(1,2,1,3), cex= 2)

#lines(log10(sAfst[[1]]$effMig[,2]), col= cols[3], lwd= 2)
#lines(log10(sAfst[[1]]$maxEffMigSbar), col= cols[3], lwd= 2.2, lty= 3)
#, expression('m'[e]) 

############
# get thetas and phis for each respective set (A - F)

psA <- xtractPhis(sAsym[1,], maf= 0.025)
psB <- xtractPhis(sBsym, maf= 0.025)
psC <- xtractPhis(sCsym, maf= 0.025)
##
psD <- xtractPhis(sDsym, maf= 0.025)
psE <- xtractPhis(sEsym, maf= 0.025)
psF <- xtractPhis(sFsym, maf= 0.025)

psls <- list(psA, psB, psC, psD, psE, psF)


	
##	nonlinear least squares fit of avg.deltaP~theta (Barton's cc)
# can either be run with a loop across all six sets or individually (which might be better for starting values for each set)

fitlsTheta <- list()

# simply run this for all psls[[1:6]] and change the list integer respectively (here and below where it's written to fitlsTheta)
	ps <- psls[[1]]
	
	theta <- unlist(lapply(ps, '[[', 2))
	deltaPn <- unlist(lapply(ps, '[[', 4))
	deltaPs <- unlist(lapply(ps, '[[', 3))
	
	x <- as.numeric(log10(theta))
	ys <- as.numeric(deltaPs)
	yn <- as.numeric(deltaPn)
	
	fitS <- nls(ys~1/(1 + exp(-a * (x-b))), start= list(a= 3, b= -0.5), na.action= na.omit)
	
	#fitN <- nls(yn~1/(1 + exp(-a * (x-b))), start= list(a= 1, b= 2), na.action= na.omit)
	fitN <- nls(yn~1/(1 + exp(-a * (x-b))), start= list(a= 2, b= 0.5), na.action= na.omit)

	fitlsTheta[[1]] <- list(fitS, fitN)
##}


	parS <- coef(fitS)
	parN <- coef(fitN)
	
	#ys2 <- sigmoid(parS, x)
	#yn2 <- sigmoid(parN, x)

sTimeHighestSlope <- c(-0.6311687, -0.6861583, -0.6405026, 1.5250583, -0.2157905, -0.2408711)
nTimeHighestSlope <- c(0.01410054, 0.33841506, 0.85201709, 86.25805769, 0.36454080, 0.74927278)
#

# figure 4
pdf('thetaActual.pdf', width= 16, height= 9)
par(mfrow= c(2,3))

for(i in 1:length(psls)){
	plot(1:10, type= 'n', xlim= c(-2, 0.75), ylim= c(0, 1), cex.lab= 1.8, xlab= expression(paste('log'[10], ' ', theta)), ylab= expression(paste(bar(p)[2], '- ', bar(p)[1])), cex.axis= 1.4)
	thetas <- psls[[i]]
	for(j in length(phis)){
	theta <- unlist(lapply(phis, '[[', 2))
	deltaPs <- as.numeric(unlist(lapply(phis, '[[', 3)))
	deltaPn <- as.numeric(unlist(lapply(phis, '[[', 4)))


		points(log10(phis[[j]]$theta_sBar), phis[[j]]$afDiffS, col= cols[1], pch= 20, cex= 1.3)
		points(log10(phis[[j]]$theta_sBar), phis[[j]]$afDiffN, col= cols[2], pch= 20)
		abline(v= sTimeHighestSlope[i], col= 'grey60', lwd= 1.3)
		abline(v= nTimeHighestSlope[i], col= 'grey60', lwd= 1.3)
	}
}
dev.off()

############
# get LD for each respective set (A - F)

ldsA <- xtractLD(sAsym, maf= 0.025)
ldsB <- xtractLD(sBsym, maf= 0.025)
ldsC <- xtractLD(sCsym, maf= 0.025)

ldsD <- xtractLD(sDsym, maf= 0.025)
ldsE <- xtractLD(sEsym, maf= 0.025)
ldsF <- xtractLD(sFsym, maf= 0.025)

ldls <- list(ldsA, ldsB, ldsC, ldsD, ldsE, ldsF)

#	nonlinear least squares fit of LD

lds <- list()
ldn <- list()
for(i in 1:length(ldls)){
	ldObj <- ldls[[i]]
	ldsels <- list()
	ldneuts <- list()
	for(j in 1:length(ldObj)){
		ldsels[[j]] <- ldObj[[j]]$LDsel[,4]
		ldneuts[[j]] <- ldObj[[j]]$LDneut[,4]
	}
	
	ldDFs <- as.data.frame(ldsels)
	ldDFs <- t(ldDFs)
	ldDFn <- as.data.frame(ldneuts)
	ldDFn <- t(ldDFn)
	
	sMed <- apply(ldDFs, 1, median, na.rm= T)
	nMed <- apply(ldDFn, 1, median, na.rm= T)
	lds[[i]] <- as.numeric(sMed) 
	ldn[[i]] <- as.numeric(nMed)
}

# can either be run with a loop across all six sets or individually (which might be better for starting values for each set)
fitld <- list()
#for(i in 1:length(lds)){

# simply run this for all ldls[[1:6]]
	ys <- lds[[1]]
	yn <- ldn[[1]]
	x <- 1:length(ys)
	
	fitS <- nls(ys~z/(1 + exp(-a * (x-b))), start= list(z= 1, a= 2, b= 10), na.action= na.omit)
	fitN <- nls(yn~z/(1 + exp(-a * (x-b))), start= list(z= 1, a= 0.5, b= 30), na.action= na.omit)
	fitld[[1]] <- list(fitS, fitN)
#}


ldsHighSlope <- c(81.12, 16.927, 28.774, NA, 161.3, 37.05)
ldnHighSlope <- c(231.6, 80.5, 166, NA, 234.6, 90.6)





pdf('figs/LDblack_ablne.pdf', width= 16, height= 9)

par(mfrow= c(2,3))	#, mar= c(5,4,4,8) + 0.1)

for(i in 1:length(ldls)){
	ldObj <- ldls[[i]]
	ldsels <- list()
	ldneuts <- list()
	for(j in 1:length(ldObj)){
		ldsels[[j]] <- ldObj[[j]]$LDsel[,4]
		ldneuts[[j]] <- ldObj[[j]]$LDneut[,4]
	}
	
	ldDFs <- as.data.frame(ldsels)
	ldDFs <- t(ldDFs)
	ldDFn <- as.data.frame(ldneuts)
	ldDFn <- t(ldDFn)
	
	sMed <- apply(ldDFs, 1, median, na.rm= T)
	nMed <- apply(ldDFn, 1, median, na.rm= T)
	qUps <- apply(ldDFs, 1, quantile, prob= 0.95, type= 8, na.rm= T)
	qUpn <- apply(ldDFn, 1, quantile, prob= 0.95, type= 8, na.rm= T)
	qLos <- apply(ldDFs, 1, quantile, prob= 0.05, type= 8, na.rm= T)
	qLon <- apply(ldDFn, 1, quantile, prob= 0.05, type= 8, na.rm= T)
	xaxs <- list(c(0, 200, 400, 600), c(0, 50, 100, 150), c(0, 10, 20, 30, 40, 50), c(0, 500, 1000, 1500), c(0, 100, 200, 300, 400), c(0, 20, 40, 60, 80))

	xax  <- xaxs[[i]]*1000/tsFreq[i]	
	xLab <- expression(paste("generations (x 10" ^"3)"))		#, tsFreq[i], ')')
	#else{

	#
	plot(x= 1:length(sMed), type= 'n', cex.axis= 1.4, cex.lab= 1.8, xlab= xLab, ylab= 'LD', xlim= c(0, length(sMed)), ylim= c(0,1), axes= F)
	abline(v= ldsHighSlope[i], col= 'grey70', lwd= 1.3)
	abline(v= ldnHighSlope[i], col= 'grey70', lwd= 1.3)

	axis(1, at= xax, labels= xaxs[[i]], cex.axis= 1.4)
	axis(2, at= c(0, 0.5, 1.0), cex.axis= 1.4)

	# nGen[i]/tsFreq[i]
	lines(qUps, lty= 1, cex= 1.2, col= 'grey50')
	lines(qLos, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])
	
	lines(qUpn, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])
	lines(qLon, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])
	lines(sMed, lty= 1, cex= 1.2)	#, col= cols[1])
	lines(nMed, lty= 1, cex= 1.2)	#, col= cols[2])
	box()
	}

dev.off()



