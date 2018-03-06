


####  
source("couplingFuncs.r")

df <- read.table("~/flaxmans/bu2s/runs/paramsALL.txt", header= T, sep= '\t')
names(df) <- tolower(names(df))
#



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

symls <- list(sAsym, sBsym, sCsym, sDsym, sEsym, sFsym)





############
# get LD for each respective set (A - F)

ldsA <- xtractLD(sAsym, setname= 'sM', folder= 'sM', path= paths)
ldsB <- xtractLD(sBsym, setname= 'sm', folder= 'sm', path= paths)
ldsC <- xtractLD(sCsym, setname= 'Sm', folder= 'Sm', path= paths)

ldsD <- xtractLD(sDsym, setname= 'sM', folder= 'sM', path= paths)
ldsE <- xtractLD(sEsym, setname= 'sM', folder= 'sM', path= paths)
ldsF <- xtractLD(sFsym, setname= 'sM', folder= 'sM', path= paths)

ldls <- list(ldsA, ldsB, ldsC, ldsD, ldsE, ldsF)

#LDls <- list()

#for(i in 1:length(ldls)){
#	for(j in 1:length(ldls[[i]])){


############
# get Fst for A - F




paths <- "/media/schimar/dapperdata/bu2s/h5/"


sEfst <- xtractFst(sEsym, setname= 'sM', folder= 'sM', path= paths)

sAfst <- xtractFst(sAsym, setname= 'sM', folder= 'sM', path= paths)


cols <- c("#c7485a", "#237fa9", "#4a921e")

# for each of the 6 sets (A-F), run the 
fobj <- sAfst

sels <- list()
neuts <- list()

for(i in 1:length(fobj)){
	sels[[i]] <- fobj[[1]]$FSTs$FSTsel
	neuts[[i]] <- fobj[[1]]$FSTs$FSTneut
}


#source('../bu2s_utils/as.data.frame.list.R')
DFs <- as.data.frame(sels)
DFs <- t(DFs)

DFn <- as.data.frame(neuts)
DFn <- t(DFn)

sMed <- apply(DFs, 1, median, na.rm= T)
nMed <- apply(DFn, 1, median, na.rm= T)
qUps <- apply(DFs, 1, quantile, prob= 0.95, type= 8, na.rm= T)
qUpn <- apply(DFn, 1, quantile, prob= 0.95, type= 8, na.rm= T)
qLos <- apply(DFs, 1, quantile, prob= 0.05, type= 8, na.rm= T)
qLon <- apply(DFn, 1, quantile, prob= 0.05, type= 8, na.rm= T)



#
xLab <- expression(paste("generations (x 10" ^"3",")"))


par(mfrow= c(1,2))

plot(1:length(sMed), type= 'n', cex.lab= 1.8, xlab= xLab, ylab= expression(F[ST]), xlim= c(0, 620), ylim= c(0,1), axes= F)
#
#axis(1, at= xax, labels= xaxs[[i]], cex.axis= 1.4)
axis(1, at= c(0, 100, 200, 300, 400, 500, 600), cex.axis= 1.4)
axis(2, at= c(0, 0.5, 1.0), cex.axis= 1.4)


lines(sMed, lty= 1, cex= 1.2, col= cols[1], lwd= 2.0)
lines(nMed, lty= 1, cex= 1.2, col= cols[2], lwd= 2.0)
#lines(qUps, lty= 1, cex= 1.2, col= 'grey50')
#lines(qLos, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])

#lines(qUpn, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])
#lines(qLon, lty= 1, cex= 1.2, col= 'grey50')	#cols[1])

box()
###


plot(log10(sAfst[[1]]$phiObs), type= 'l', ylab= expression(paste('log'[10], ' metric')), cex.lab= 1.8, cex.axis= 1.4, lwd= 2, xlab= xLab, ylim= c(-2, 5))
lines(log10(sAfst[[1]]$sStarLeS$Le), lwd= 2, lty= 2)

legend('topleft', legend= c(expression(paste(phi)), expression('L'[e])), lty= c(1,2,1,3), cex= 2)

#lines(log10(sAfst[[1]]$effMig[,2]), col= cols[3], lwd= 2)
#lines(log10(sAfst[[1]]$maxEffMigSbar), col= cols[3], lwd= 2.2, lty= 3)
#, expression('m'[e]) 


############
# get thetas and phis for each respective set (A - F)

psA <- xtractPhis(sAsym, setname= 'sM', folder= 'sM', maf= 0.025, path= paths)
psB <- xtractPhis(sBsym, setname= 'sm', folder= 'sm', maf= 0.025, path= paths)
psC <- xtractPhis(sCsym, setname= 'Sm', folder= 'Sm', maf= 0.025, path= paths)
##
psD <- xtractPhis(sDsym, setname= 'sM', folder= 'sM', maf= 0.025, path= paths)
psE <- xtractPhis(sEsym, setname= 'sM', folder= 'sM', maf= 0.025, path= paths)
psF <- xtractPhis(sFsym, setname= 'sM', folder= 'sM', maf= 0.025, path= paths)



