# ZERO CROSSING FUNCTION --------------------------------------------
zc <- function(p_d, inc, start){

  if (isTRUE(inc)){
    s<-which.min(p_d[start:(start+300)])+start
    r<-which.max(p_d[s:(s+100)])+s
    print("Incisura detected")
  }else{
    s <- start-100
    r<-which.max(p_d[s:(s+99)])+s
  }

  rm<-which.min(p_d[r:(r+50)])+r
  z<-which.min(abs(p_d[r:rm]))+r

  #   abline(v=r, col="red")
  #    abline(v=z, col="blue")
  #    abline(v=rm, col="green")
  print(z)
  return (z)
  # end of function --
}

# OUTLIER DETECTION FUNCTION ----------------------------------------
out <- function (x,z){
  mx <- scale(x, scale=T, center=T)
  i <- which(mx>=z | mx<=-z)
  return(i)
# end of function --
}

# RMSE FUNCTION -----------------------------------------------------
rmse <- function(act, pred){
  act<-na.omit(act)
  pred<-na.omit(pred)

  longest <- max(length(act), length(pred))
  shortest <- min(length(act), length(pred))
  null <- rep(length.out=(longest-shortest), 0)

  if (length(act)>length(pred)){
    pred1<-c(pred,null)
    error<-act-pred1
    rmse <- sqrt(mean(error^2))
  }

  if (length(pred)>length(act)){
    act1<-c(act,null)
    error<-act1-pred
    rmse <- sqrt(mean(error^2))
  }

  if (length(pred)==length(act)){
    error<-act-pred
    rmse <- sqrt(mean(error^2))
  }
  return(rmse)
# end of function --
}

# PEAK DETECTOR FUNCTION --------------------------------------------
peak <- function(p,ecg,sr,t,dt){
  # p - pressure signal
  # ecg - ecg signal
  # sr - sampling rate (default = 1000 Hz)
  # t - threshold (default = 0.7)
  # dt - detrend (y/n)

  if(exists("t")){
    t <- t
  }else{
    t <- 0.7
  }

  if(exists("sr")){
    sr <- sr
  }else{
    sr <- 1000
  }

  if (isTRUE(dt)){
    p1 <- resid(lm(diff(p,150) ~ seq(along = diff(p,150)))) # detrending
  }else{
    p1 <- p
  }
  p2 <- ma(p1,80) # smoothing

  if (length(ecg)>100){
    # Detection based on ECG
    print("Peak detection based on ECG")
    ecg1 <- resid(lm(ecg ~ seq(along = ecg))) # ECG detrending X
    ecg2 <- diff(ecg1, 10) # ECG differencing (-10 data points) Y

    inc <- sr/10
    window <- length(ecg2)
    i <- 1
    n <- 0
    ciklus <- 0
    R <- 0

    ecg2[ecg1<=max(ecg2)*t] <- 0
    while (i<=window-1){
      if(ecg2[i]>0){
        n <- n+1
       # R <- ecg1[i:(i+inc)]
       # r.max <- max(R)
        ciklus[n] <- as.vector(which.max(ecg1[i:(i+inc)]))
        ciklus[n] <- ciklus[n]+i
        if(n>1){
          if (ciklus[n]-ciklus[n-1]<600){# minimum cycle length: 600 ms
            n<-n-1
          }
        }
        i <- i+inc
      }
      i <- i+1
    }
  }else{
    # Detection based on pressure alone
    print("Peak detection based on pressure")
    k <- max(p2, na.rm=T)*t # threshold

    p2[which(p2<k)] <- 0
    p2[which(p2>k)] <- k
    p3 <- diff(p2)
    p3[which(p3!=k)] <- 0
    p3[which(p3==k)] <- 1
    ciklus <- which(p3==1)
  }
  return (ciklus)
# end of function --
}

# AVERAGING FUNCTION ------------------------------------------------
avg <- function(p, c, rm=0){
  # p - raw pressure signal
  # c - cycle detector ouput (borders)
  # rm - outlier removal threshold (default = 0)

  if (rm!=0){
    rm <- rm
  }

  # cycle number
  cn <- length(c)-1

  # dataframe initialization
  h <- max(diff(c)) # maximal cycle length
  press <- data.frame(matrix(ncol = h+1, nrow = cn))

  # dataframe fill-up
  for (i in 1:cn){
    press[i,] <- p[c[i]:(c[i]+h)]
  }

  # average aortic waveform
  press_mean <- apply(press,2,mean,na.rm=T)

  # RMSE from the mean curve
  press.error <- 0
  for (i in 1:cn){
    press.error[i] <- rmse(as.numeric(press[i,]), press_mean) # root mean square error
  }

  # Handling with outlier cycles
  if (rm!=0){
    outliers <- out (press.error, rm)
    if(length(outliers)>0){
      press.c<-press[-outliers,]
      press_mean <- apply(press.c,2,mean,na.rm=T)
      cat("\014") # Clear screen
      print(paste(length(outliers),"outliers removed, averaged from", cn-length(outliers),"cycles"))
    }
  }else{
    cat("\014") # Clear screen
    print("Keep outliers")
  }
  return(press_mean)
# end of function --
}

# BP CALIBRATION FUNCTION -------------------------------------------
cal <- function (inp, max, mean, min, central){
  # inp - data file
  # central - T: central; F: peripheral
  mini <- min(inp[1:500])

  if(central){
    print("Central waveform scaled")
    pr.ce <- min+((inp-mini)*((mean-min)/(mean(inp)-mini)))
  }else{
    print("Peripheral waveform scaled")
    pr.ce <- min+((inp-mini)*((max-min)/(max(inp)-mini)))
  }
  return(pr.ce)
# end of function --
}

# PWA FUNCTION ------------------------------------------------------
pwa <- function(p, abs){
  if (is.na(abs)){abs <- F}

  # 4. derivative
  p.dif <- diff(p, 20, 4)

  ps <- max(p, na.rm = T) # systolic pressure
  ps.i <- min(which(p==ps)) # systole time index
  pd <- min(p[1:ps.i], na.rm = T) # diastolic pressure
  pd.i <- min(which(p[1:ps.i]==pd)) # diastole time index
  inc.i <- as.numeric(zc(p.dif,inc=T,250)) # incisura time index based on zero crossing function
  pd.i2<-max(which(p==min(p[inc.i:length(p)], na.rm = T))) # end of cycle (diastole)
  pp <- ps-pd # pulse pressure

  plot(p, type="l", bty="l", las=1, xlab="time[ms]", ylab="pressure[mmHg]", lwd=2)
  text(700,100,"Click around P1 and hit ESC")
  p1.i<-identify(p, plot=F) # identify and ESC
  if (isTRUE(abs)){
  p1.i2 <- p1.i
  print("User-defined inflexion")
  }else{
    to<-p1.i+100-1
    from<-p1.i-100

  if(to>ps.i & p1.i<ps.i){
    to <- ps.i
  }

  if(from<ps.i & p1.i>ps.i){
    from <- ps.i
  }

  if(from<80){
    from <- 80
  }


  if(p1.i-100<ps.i & p1.i>ps.i){p1.i <- ps.i+100}

  p1.i2 <- zc(p.dif,inc=F,p1.i) # p1 detection based on zero crossing
  #p1.i2 <- zc((p.dif[(from-80):(to-80)]),0,0) # p1 detection based on zero crossing
  p1.i2 <- as.numeric(substr(names(p1.i2),2,4))
  if(p1.i>ps.i){p1.i2 <- p1.i2-80}
  } # else for absolute detection end

  p1<-as.numeric(p[p1.i2]-pd) # p1 pressure measured from diastole

  plot(p, type="l", bty="l", las=1, xlab="time[ms]", ylab="pressure[mmHg]", lwd=2)

  points(ps.i, ps, col="red", pch=19)
  points(pd.i, pd, col="blue", pch=19)
  points(inc.i, p[inc.i], col="darkgreen", pch=19)
  points(p1.i2, p[p1.i2], col="magenta", pch=19)
  points(pd.i2, p[pd.i2], col="darkblue", pch=19)

  # Augmentation index
  if(p1<=pp & p1.i2<ps.i){
    ap <- pp-p1
    if(ap<=12) print("B type")
    if(ap>12) print("A type")
  }else{
    ap <- p1-pp
    print("C type")
  }
  aix <- 100*ap/pp

  # 75/min HR normalized Aix
  hr <- 60/((pd.i2-pd.i)/1000)
  aix.75 <- aix-(0.39*(75-hr))

  # SEVR determination
  tti<-sum(p[pd.i:inc.i]) # systolic AUC
  dti<-sum(p[inc.i:pd.i2]) # distolic AUC
  sevr<-dti/tti # aka Buckberg-index

  # Form factor
  ff <- mean(p, na.rm = T)/pp
  pressures <- list("SBPc"=round(ps,digits=2), "DBPc"=round(pd,digits=2), "PPc"=round(pp,digits=2), "P1"=round(p1,digits=2), "tP1"=round(p1.i2, digits = 2), "AP"=round(ap,digits=2), "AiX"=round(aix,digits=2), "AiX@75"=round(aix.75,digits=2), "Form"=round(ff,digits=2), "T"=round((pd.i2-pd.i)/1000,digits=2), "HR"=round(hr,digits=2), "incisura"=round(inc.i,digits=2), "start"=round(pd.i, digits = 2) , "end"=round(pd.i2, digits = 2),"TTI"=round(tti, digits = 2), "DTI"=round(dti, digits = 2),"SEVR"=round(sevr, digits = 2))

# end of function --
}

# WAVE-FREE PERIOD DETECTOR FUNCTION --------------------------------
wfp <- function(p, inc, end){
  # p - averaged, scaled pressure wave
  # inc - location of incisura (end-systole)
  # end - end of the pressure wave

  # wave-free period beginning-end
  wfp.i<-c(inc+(0.25*(end-inc)), end-5)

  # plot
  xx <- c(wfp.i[1],wfp.i[1]:wfp.i[2],wfp.i[2])
  yy <- c(pac[wfp.i[2]],pac[wfp.i[1]:wfp.i[2]],pac[wfp.i[2]])
  polygon(xx, yy, col="pink", lty="blank")

  da <- sum(pac[wfp.i[1]:wfp.i[2]])
  tau <- da/(pac[wfp.i[1]]-pac[wfp.i[2]])/2000

  # wave-free period
  wfp<-p[wfp.i[1]:wfp.i[2]]

  return(as.numeric(c(wfp.i, tau)))
# end of function --
}

# INTERSECTING TANGENTS FUNCTION ------------------------------------
foot <- function (bp, c, d){
  # bp - raw blood pressure recording
  # c - ECG-based cycle delimiter
  # d - draw (y/n)
  p.s <- sgolayfilt(bp, n=21) # smoothing filter
  p.s.d <- diff(p.s, 20) # differentiation by 20 points
  l <- length(c)
  m <- 0
  d.i <- 0
  dp.max <- 0
  for (i in 1:l){
    m <- which.max(p.s[c[i]:(c[i]+300)])+c[i]
    d.i <- c(d.i, which.min(p.s[c[i]:m])+c[i]) # most negative values after ECG and before systole
    dp.max <- c(dp.max, which.max(p.s.d[c[i]:(c[i]+300)])+c[i]) # maximum of the derivative
  }
  # correction of first zeros
  d.i <- d.i[-1]
  dp.max <- dp.max[-1]

  dp.max <- dp.max+20 # correction for the lag due to differentiation
  d.min <- p.s[d.i]

  koef <- data.frame(matrix(ncol = 2, nrow = l)) # dataframe initialization
  colnames(koef) <- c("a","m")
  isect <- 0

  for(i in 1:l){
    k <- 1
    for (ii in 1:15){
      # pick the best fitting regressional line
      if(as.numeric(summary(lm(p.s[(dp.max[i]-ii):(dp.max[i]+ii)]~seq(2*ii+1)))$r.squared)>k){
        k <- ii
      }
      # regression coefficients for each cardiac cycle
      koef[i,]<- as.numeric(coef(lm(p.s[(dp.max[i]-k):(dp.max[i]+k)]~seq(2*k+1))))
      isect[i] <- (d.min[i]-koef[i,1])%/%koef[i,2] # integer definition of the instersection
    }
  }
  it <- dp.max+isect # location of the intersection relative to the maximum first derivative
  ptt <- it-c

  if(isTRUE(d)){
    # Drawings if requested
    par(mfrow=c(3,1), mar=c(2,2,2,2))
    plot(adat$V1, type="l", col="darkgreen", bty="l")
    points(c,adat$V1[c(c)], pch=16, col="red")

    plot(bp, type="l", bty="l")
    points(it,bp[c(it)], pch=4, col="blue")
    points(c,bp[c(c)], pch=16, col="red")

    plot(ptt, pch=15, col="blue", bty="l", xaxt="n", ylim=c(0.98*min(ptt), 1.02*max(ptt)))
    abline(h=mean(ptt)+(1.96*sd(ptt)), lty=2, col="red")
    abline(h=mean(ptt)-(1.96*sd(ptt)), lty=2, col="red")
    abline(h=mean(ptt), lty=2, col="darkblue")
  }
return(ptt)
}

# PROCESS FUNCTION --------------------------------------------------
proc <- function(adat, ekg, thres, freq, dt){
  # adat - raw, unprocessed data file
  # central - peripheral or central pressure processing
  # ecg - pressure/ECG-based cycle recognition
  # thres - peak detector threshold
  # freq - smapling rate
  # dt - detrending
  p<-adat$V2
  ecg<-adat$V1

  if (is.na(thres)){
    thres <- 0.9
  }

  if (is.na(freq)){
    freq <- 1000
  }

  if (isTRUE(dt)){
    dt <- T
  }

  if (isTRUE(ekg)){
    # Call with ECG, 1000 Hz, treshold=0.7, detrending before analysis
    c <- peak(p, ecg, freq, thres, dt)
  }else{
    # Call without ECG, 1000 Hz, treshold=0.6, detrending before analysis
    c <- peak(p, ecg="F", freq, thres, dt)
  }

  # Some plots
  par(mfrow=c(2,1), mar=c(4,4,2,2))
  plot(ecg, type="l", col="darkgreen", bty="l")
  abline(v=c, col="red")

  plot(p, type="l", bty="l")
  abline(v=c, col="blue")
  return (c)

  # end function --
}

# FFR FUNCTION --------------------------------------------------
proc <- ffr(pp, pd, draw){
  # raw proximal pressure waveform
  # raw distal pressure waveform
  # draw the results (y/n)

  # Savitzky-Golay smoothing
  pp<-sgolayfilt(na.omit(pp),p=3,n=21)
  pd<-sgolayfilt(na.omit(pd),p=3,n=21)

  auto <- ccf(pp, pd, plot=F, lag.max = 100) # cross correlation
  l <- auto$lag[which(auto$acf==max(auto$acf))] # optimal lag
  print (paste("Lag =", l))

  pp1 <- tail(pp,length(pd)-l)
  pd1 <- head(pd,length(pd)-l)

  cikl <- peak(pd1, 0, 1000, 0.7, dt)
  # Ciklusok szama
  cn <- length(cikl)
  h <- max(diff(cikl)) # a leghosszabb ciklus
  dia<-0
  for (i in 1:cn){
    if (cikl[i]-200<0) {st<-0}else{st<-cikl[i]-200}
    d<-which.min(pd1[st:cikl[i]])
    dia<-c(dia, st+d)
  }
  dia<-dia[-1]

  # dataframe initialization
  pres_m <- data.frame(matrix(ncol = 4, nrow = cn))
  colnames(pres_m) <- c("prox_mean","dist_mean", "FFR", "iFR")

  for (i in 1:cn){
    # incisura detektalas
    inc.i.p <- zc(diff(pp1[dia[i]:(dia[i]+h)], 20, 4),250,0)
    inc.i.d <- zc(diff(pd1[dia[i]:(dia[i]+h)], 20, 4),250,0)

    # ciklus-vege diastole
    pp.d.i<-which(pp1==min(pp1[inc.i.p:h], na.rm = T))
    pd.d.i<-which(pd1==min(pd1[inc.i.d:h], na.rm = T))

    pres_m[i,1] <- mean(pp1[dia[i]:(dia[i]+h)])
    pres_m[i,2] <- mean(pd1[dia[i]:(dia[i]+h)])
    if (i==1){
      pres_m[i,3] <- NA
    }else{
      pres_m[i,3] <- pres_m[i,2]/pres_m[i,1]
    }
    pres_m[i,4] <- wfp(pp1[dia[i]:(dia[i]+h)], inc.i.p, pp.d.i)
  }


# Files read-in -----------------------------------------------------


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
