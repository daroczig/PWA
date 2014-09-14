# Files read-in -----------------------------------------------------
library(tcltk)

# setwd("~/Documents/Munka/K3/2013 - IKER2/PWV")
 setwd("C:/K3/Projektek/2013 - IKER2/PWV")

# Cuff-measurements source file
brachialis <- read.delim("iker2_blood_pressure.txt", dec=",")
brachialis <- brachialis[-c(4,5)]
colnames(brachialis)<-c("id","sbp","dbp")

# Open recordings for one patient
filename <- tk_choose.files()
filename <- strsplit(basename(filename)," - ")
fname <- substring(filename,1, nchar(filename)-4)
id <- unique(substr(fname, 1,6))

# file sorrend: _b,_c,_f,_r

sbpb <- brachialis$sbp[brachialis$id==id]
dbpb <- brachialis$dbp[brachialis$id==id]
ppb <- sbpb-dbpb
mbpb <- dbpb+(ppb)/3




# Open the brachial recording
adat <- read.table(as.character(filename[1]), header=F, sep = "\t", dec=".")
c <- proc(adat, ekg=T, 0.8, 1000, dt)

ptt.b <- foot(adat$V2, c, d=T)
if(out(ptt.b,1.6)>0){ ptt.b <- ptt.b[-c(out(ptt.b,1.6))]}
ptt.b <- mean(ptt.b)

# correction for wrong cycles
#c <- c[-c(1,4)]
#c <- c[-length(c)]

pa <- avg(adat$V2, c, rm=1.2)
dev.off()
pac <- cal(pa, sbpb,0,dbpb, central=F)
pressures.b <- list("ID"=id, "SBPb"=round(sbpb, digits = 2), "DBPb"=round(dbpb, digits = 2),
"PPb"=round(ppb, digits = 2), "MBPb"=round(mbpb, digits = 2), "MBPbI"=round(mean(pac), digits = 2))
plot(pac, type="l", bty="l", las=1, xlab="time[ms]", ylab="pressure[mmHg]", lwd=2, main=fname[1])
brach.s <- pac


# Open the carotid recording
adat <- read.table(as.character(filename[2]), header=F, sep = "\t", dec=".")
c <- proc(adat, ekg=T, 0.9, 1000, dt)
# correction for wrong cycles
# c <- c[-c(1,11)]
ptt.c <- foot(adat$V2, c, d=T)
if(out(ptt.c,1.6)>0){ptt.c <- ptt.c[-c(out(ptt.c,1.6))]}
ptt.c <- mean(ptt.c)

pa <- avg(adat$V2, c, rm=1)
pac <- cal(pa, sbpb,mean(pac),dbpb, central=T)
dev.off()

pressures <- c(pressures.b, pwa(pac,abs=F)) # semi-automatic P1 detection
pressures <- c(pressures.b, pwa(pac,abs=T)) # manual P1 detection

wfp.l <- wfp(pac, pressures$incisura, pressures$end)
pressures<-c(pressures,"tau"=round(wfp.l[3],digits = 2))
car.s <- pac


# Open the femoral recording
adat <- read.table(as.character(filename[3]), header=F, sep = "\t", dec=".")
c <- proc(adat, ekg=T, 0.9, 1000, dt)
ptt.f <- foot(adat$V2, c, d=T)
if(out(ptt.f,1.6)>0){ptt.f <- ptt.f[-c(out(ptt.f,1.6))]}
ptt.f <- mean(ptt.f)

pa <- avg(adat$V2, c, rm=1.2)
pac <- cal(pa, sbpb,pressures$MBPbI,dbpb, central=T)
dev.off()
plot(pac, type="l", bty="l", las=1, xlab="time[ms]", ylab="pressure[mmHg]", lwd=2, main=fname[3])
fem.s <- pac


# Open the radial recording
adat <- read.table(as.character(filename[4]), header=F, sep = "\t", dec=".")
c <- proc(adat, ekg=T, 0.9, 1000, dt)
# correction for wrong cycles
#c <- c[-length(c)]
ptt.r <- foot(adat$V2, c, d=T)
if(out(ptt.r,1.6)>0){ptt.r <- ptt.r[-c(out(ptt.r,1.6))]}
ptt.r <- mean(ptt.r)

pa <- avg(adat$V2, c, rm=1.2)
pac <- cal(pa, sbpb,pressures$MBPbI,dbpb, central=T)
dev.off()
plot(pac, type="l", bty="l", las=1, xlab="time[ms]", ylab="pressure[mmHg]", lwd=2, main=fname[4])
rad.s <- pac

ptt <- list("Brachial TT"=round(ptt.b, digits = 2), "Carotid TT"=round(ptt.c, digits = 2), "Femoral TT"=round(ptt.f, digits = 2), "Radial TT"=round(ptt.r, digits = 2))
cat("\014") # Clear screen
print(paste("Results for subject", id))
pressures<-c(pressures,ptt)
t(pressures)


#Save summary data
longest <- max(length(brach.s), length(car.s), length(fem.s), length(rad.s))
summary <- data.frame(matrix(ncol = 4, nrow = longest))
summary[1] <- c(brach.s, rep(length.out=(longest-length(brach.s)), NA))
summary[2] <- c(car.s, rep(length.out=(longest-length(car.s)), NA))
summary[3] <- c(fem.s, rep(length.out=(longest-length(fem.s)), NA))
summary[4] <- c(rad.s, rep(length.out=(longest-length(rad.s)), NA))
colnames(summary)=c("brachial","carotid","femoral","radial")

write.table(summary, file=paste(id,"_waves.csv"), sep = ";", dec=",", col.names = T, qmethod = "double", na="NA", row.names=FALSE)
write.table(pressures, file=paste(id,"_SUM.csv"), sep = ";", dec=",", col.names = T, qmethod = "double", na="-999", row.names=FALSE)
