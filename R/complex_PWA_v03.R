#' Zero crossing function
#' @param p_d p_d
#' @param inc inc
#' @param start start
#' @return vector
#' @export
zc <- function(p_d, inc, start) {

  if (isTRUE(inc)) {
      s <- which.min(p_d[start:(start + 300)]) + start
      r <- which.max(p_d[s:(s + 100)]) + s
      message("Incisura detected")
  } else {
      s <- start - 100
      r <- which.max(p_d[s:(s + 99)]) + s
  }

  rm <- which.min(p_d[r:(r + 50)]) + r
  z  <- which.min(abs(p_d[r:rm])) + r

  z
}


#' Root-mean-square error function
#' @param act act
#' @param pred pred
#' @return rmse
#' @export
rmse <- function(act, pred) {

    act  <- na.omit(act)
    pred <- na.omit(pred)

    longest  <- max(length(act), length(pred))
    shortest <- min(length(act), length(pred))
    null     <- rep(length.out = longest - shortest, 0)

    if (length(act) > length(pred)) {
        pred1 <- c(pred,null)
        error <- act-pred1
        rmse  <- sqrt(mean(error^2))
    }

    if (length(pred) > length(act)) {
        act1  <- c(act, null)
        error <- act1 - pred
        rmse  <- sqrt(mean(error^2))
    }

    if (length(pred) == length(act)){
        error <- act - pred
        rmse  <- sqrt(mean(error^2))
    }

    rmse
}


#' Peak detector functoin
#' @param p pressure signal
#' @param ecg ecg signal
#' @param sr sampling rate (Hz)
#' @param t threshold
#' @param dt detrend
#' @return value
#' @export
peak <- function(p,ecg,sr = 1000, t = 0.7, dt = TRUE){

    if (isTRUE(dt)) {
        p1 <- resid(lm(diff(p, 150) ~ seq(along = diff(p, 150)))) # detrending
    } else {
        p1 <- p
    }
    p2 <- ma(p1, 80) # smoothing

    if (length(ecg)>100){

        ## Detection based on ECG
        print("Peak detection based on ECG")
        ecg1 <- resid(lm(ecg ~ seq(along = ecg))) # ECG detrending X
        ecg2 <- diff(ecg1, 10) # ECG differencing (-10 data points) Y

        inc <- sr / 10
        window <- length(ecg2)
        i <- 1
        n <- 0
        ciklus <- 0
        R <- 0

        ecg2[ecg1 <= max(ecg2) * t] <- 0
        while (i <= window - 1){
            if (ecg2[i] > 0) {
                n <- n + 1
                ## R <- ecg1[i:(i+inc)]
                ## r.max <- max(R)
                ciklus[n] <- as.vector(which.max(ecg1[i:(i + inc)]))
                ciklus[n] <- ciklus[n] + i
                if (n > 1) {
                    if (ciklus[n] - ciklus[n - 1] < 600){# minimum cycle length: 600 ms
                        n <- n - 1
                    }
                }
                i <- i + inc
            }
            i <- i + 1
        }
    } else {

        ## Detection based on pressure alone
        message("Peak detection based on pressure")
        k <- max(p2, na.rm = TRUE) * t # threshold

        p2[which(p2 < k)] <- 0
        p2[which(p2 > k)] <- k
        p3 <- diff(p2)
        p3[which(p3 != k)] <- 0
        p3[which(p3 == k)] <- 1
        ciklus <- which(p3 == 1)

    }

    ciklus

}


#' Averaging functoin
#' @param p raw pressure signal
#' @param c cycle detector ouput (borders)
#' @param rm outlier removal threshold
#' @return value
#' @export
avg <- function(p, c, rm = 0){

    ## cycle number
    cn <- length(c) - 1

    ## dataframe initialization
    h <- max(diff(c)) # maximal cycle length
    press <- data.frame(matrix(ncol = h + 1, nrow = cn))

    ## dataframe fill-up
    for (i in 1:cn){
        press[i, ] <- p[c[i]:(c[i] + h)]
    }

    ## average aortic waveform
    press_mean <- apply(press, 2, mean, na.rm = TRUE)

    ## RMSE from the mean curve
    press.error <- 0
    for (i in 1:cn){
        press.error[i] <- rmse(as.numeric(press[i, ]), press_mean) # root mean square error
    }

  ## Handling with outlier cycles
    if (rm != 0){
        outliers <- out(press.error, rm)
        if (length(outliers) > 0) {
            press.c <- press[-outliers, ]
            press_mean <- apply(press.c, 2, mean, na.rm = TRUE)
            ## cat("\014") # Clear screen
            print(paste(length(outliers),"outliers removed, averaged from", cn-length(outliers),"cycles"))
        }
    } else {
        ## cat("\014") # Clear screen
        print("Keep outliers")
    }
    return(press_mean)
}


#' BP calibration function
#' @param inp data file
#' @param max max
#' @param mean mean
#' @param min min
#' @param central central or peripheral
#' @return value
#' @export
cal <- function (inp, max, mean, min, central = TRUE){

    mini <- min(inp[1:500])

    if (central) {
        print("Central waveform scaled")
        pr.ce <- min + ((inp - mini) * ((mean - min) / (mean(inp) - mini)))
    } else {
        print("Peripheral waveform scaled")
        pr.ce <- min + ((inp - mini) * ((max - min) / (max(inp) - mini)))
    }
  pr.ce
}


#' Pulse Wave Analysis
#' @param p p
#' @param abs abs
#' @return value
#' @export
pwa <- function(p, abs) {

    if (is.na(abs))
        abs <- F

    ## 4. derivative
    p.dif <- diff(p, 20, 4)

    ps    <- max(p, na.rm = TRUE) # systolic pressure
    ps.i  <- min(which(p == ps)) # systole time index
    pd    <- min(p[1:ps.i], na.rm = TRUE) # diastolic pressure
    pd.i  <- min(which(p[1:ps.i] == pd)) # diastole time index
    inc.i <- as.numeric(zc(p.dif, inc = TRUE, 250)) # incisura time index based on zero crossing function
    pd.i2 <- max(which(p == min(p[inc.i:length(p)], na.rm = TRUE))) # end of cycle (diastole)
    pp    <- ps - pd # pulse pressure

    plot(p, type = "l", bty = "l", las = 1, xlab = "time[ms]", ylab = "pressure[mmHg]", lwd = 2)
    text(700, 100, "Click around P1 and hit ESC")
    p1.i <- identify(p, plot=FALSE) # identify and ESC

    if (isTRUE(abs)) {

        p1.i2 <- p1.i
        print("User-defined inflexion")

    }else{

        to <- p1.i + 100 - 1
        from <- p1.i - 100

        if (to > ps.i & p1.i < ps.i)
            to <- ps.i

        if (from < ps.i & p1.i > ps.i)
            from <- ps.i

        if (from<80)
            from <- 80

        if (p1.i - 100 < ps.i & p1.i > ps.i)
            p1.i <- ps.i + 100

        p1.i2 <- zc(p.dif, inc = FALSE, p1.i) # p1 detection based on zero crossing
        ##p1.i2 <- zc((p.dif[(from-80):(to-80)]),0,0) # p1 detection based on zero crossing
        p1.i2 <- as.numeric(substr(names(p1.i2), 2, 4))
        if (p1.i > ps.i)
            p1.i2 <- p1.i2 - 80
    } # else for absolute detection end

    p1 <- as.numeric(p[p1.i2] - pd) # p1 pressure measured from diastole

    plot(p, type = "l", bty = "l", las = 1, xlab = "time[ms]", ylab = "pressure[mmHg]", lwd = 2)

    points(ps.i, ps, col = "red", pch = 19)
    points(pd.i, pd, col = "blue", pch = 19)
    points(inc.i, p[inc.i], col = "darkgreen", pch = 19)
    points(p1.i2, p[p1.i2], col = "magenta", pch = 19)
    points(pd.i2, p[pd.i2], col = "darkblue", pch = 19)

    ## Augmentation index
    if (p1 <= pp & p1.i2 < ps.i) {

        ap <- pp - p1
        if (ap <= 12)
            print("B type")
        if (ap > 12)
            print("A type")

    } else {
        ap <- p1 - pp
        print("C type")
    }

    aix <- 100 * ap / pp

    ## 75/min HR normalized Aix
    hr <- 60 / ((pd.i2 - pd.i) / 1000)
    aix.75 <- aix - (0.39 * (75 - hr))

    ## SEVR determination
    tti  <- sum(p[pd.i:inc.i]) # systolic AUC
    dti  <- sum(p[inc.i:pd.i2]) # distolic AUC
    sevr <- dti / tti # aka Buckberg-index

    ## Form factor
    ff <- mean(p, na.rm = TRUE) / pp
    pressures <- list(
        "SBPc"     = round(ps, digits = 2),
        "DBPc"     = round(pd, digits = 2),
        "PPc"      = round(pp, digits = 2),
        "P1"       = round(p1, digits = 2),
        "tP1"      = round(p1.i2, digits = 2),
        "AP"       = round(ap, digits = 2),
        "AiX"      = round(aix, digits = 2),
        "AiX@75"   = round(aix.75, digits = 2),
        "Form"     = round(ff, digits = 2),
        "T"        = round((pd.i2 - pd.i) / 1000, digits = 2),
        "HR"       = round(hr, digits = 2),
        "incisura" = round(inc.i,digits = 2),
        "start"    = round(pd.i, digits = 2),
        "end"      = round(pd.i2, digits = 2),
        "TTI"      = round(tti, digits = 2),
        "DTI"      = round(dti, digits = 2),
        "SEVR"     = round(sevr, digits = 2)
        )

    pressures

}


#' Wave-free period detector function
#' @param p averaged, scaled pressure wave
#' @param inc location of incisura (end-systole)
#' @param end end of the pressure wave
#' @return value
#' @export
wfp <- function(p, inc, end){

    ## wave-free period beginning-end
    wfp.i <- c(inc + (0.25 * (end - inc)), end - 5)

    ## plot
    xx <- c(wfp.i[1], wfp.i[1]:wfp.i[2], wfp.i[2])
    yy <- c(pac[wfp.i[2]], pac[wfp.i[1]:wfp.i[2]], pac[wfp.i[2]])
    polygon(xx, yy, col = "pink", lty = "blank")

    da <- sum(pac[wfp.i[1]:wfp.i[2]])
    tau <- da / (pac[wfp.i[1]] - pac[wfp.i[2]]) / 2000

    ## wave-free period
    wfp <- p[wfp.i[1]:wfp.i[2]]

    as.numeric(c(wfp.i, tau))

}


#' Intersecting tangents function
#' @param bp raw blood pressure recording
#' @param c ECG-based cycle delimiter
#' @param d draw
#' @return value
#' @export
foot <- function (bp, c, d){

    p.s    <- sgolayfilt(bp, n = 21) # smoothing filter
    p.s.d  <- diff(p.s, 20) # differentiation by 20 points
    l      <- length(c)
    m      <- 0
    d.i    <- 0
    dp.max <- 0

    for (i in 1:l) {
        m <- which.max(p.s[c[i]:(c[i] + 300)]) + c[i]
        d.i <- c(d.i, which.min(p.s[c[i]:m]) + c[i]) # most negative values after ECG and before systole
        dp.max <- c(dp.max, which.max(p.s.d[c[i]:(c[i] + 300)]) + c[i]) # maximum of the derivative
  }

    ## correction of first zeros
    d.i <- d.i[-1]
    dp.max <- dp.max[-1]

    dp.max <- dp.max + 20 # correction for the lag due to differentiation
    d.min <- p.s[d.i]

    koef <- data.frame(matrix(ncol = 2, nrow = l)) # dataframe initialization
    colnames(koef) <- c("a","m")
    isect <- 0

    for(i in 1:l) {
        k <- 1
        for (ii in 1:15){

            ## pick the best fitting regressional line
            if (as.numeric(summary(lm(p.s[(dp.max[i] - ii):(dp.max[i] + ii)] ~ seq(2 * ii + 1)))$r.squared) > k) {
                k <- ii
            }

            ## regression coefficients for each cardiac cycle
            koef[i,]<- as.numeric(coef(lm(p.s[(dp.max[i] - k):(dp.max[i] + k)] ~ seq(2 * k + 1))))
            isect[i] <- (d.min[i] - koef[i, 1]) %/% koef[i, 2] # integer definition of the instersection
        }
    }

    it <- dp.max + isect # location of the intersection relative to the maximum first derivative
    ptt <- it-c

    if (isTRUE(d)) {

        ## Drawings if requested
        par(mfrow = c(3, 1), mar = c(2, 2, 2, 2))
        plot(adat$V1, type = "l", col = "darkgreen", bty = "l")
        points(c, adat$V1[c(c)], pch = 16, col = "red")

        plot(bp, type = "l", bty = "l")
        points(it, bp[c(it)], pch = 4, col = "blue")
        points(c, bp[c(c)], pch = 16, col = "red")

        plot(ptt, pch = 15, col = "blue", bty = "l", xaxt = "n", ylim = c(0.98 * min(ptt), 1.02 * max(ptt)))
        abline(h = mean(ptt) + (1.96 * sd(ptt)), lty = 2, col = "red")
        abline(h = mean(ptt) - (1.96 * sd(ptt)), lty = 2, col = "red")
        abline(h = mean(ptt), lty = 2, col = "darkblue")
    }

    ptt

}


#' Process function
#' @param adat raw, unprocessed data file
#' @param ekg pressure/ECG-based cycle recognition
#' @param thres peak detector threshold
#' @param freq smapling rate
#' @param dt detrending
#' @return value
#' @export
proc <- function(adat, ekg = TRUE, thres = 0.9, freq = 1000, dt = TRUE){

    p   <- adat$V2
    ecg <- adat$V1

    if (isTRUE(ekg)) {
        ## Call with ECG, 1000 Hz, treshold=0.7, detrending before analysis
        c <- peak(p, ecg, freq, thres, dt)
    } else {
        ## Call without ECG, 1000 Hz, treshold=0.6, detrending before analysis
        c <- peak(p, ecg = "F", freq, thres, dt)
    }

    ## Some plots
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
    plot(ecg, type = "l", col = "darkgreen", bty = "l")
    abline(v = c, col = "red")

    plot(p, type = "l", bty = "l")
    abline(v = c, col = "blue")

    return (c)

}


#' FFR function
#' @param pp raw proximal pressure waveform
#' @param pd raw distal pressure waveform
#' @param draw plot
#' @return value
#' @export
ffr <- function(pp, pd, draw) {

  ## Savitzky-Golay smoothing
    pp <- sgolayfilt(na.omit(pp), p = 3, n = 21)
    pd <- sgolayfilt(na.omit(pd), p = 3, n = 21)

    auto <- ccf(pp, pd, plot = FALSE, lag.max = 100) # cross correlation
    l <- auto$lag[which(auto$acf == max(auto$acf))]  # optimal lag
    print(paste("Lag =", l))

    pp1 <- tail(pp, length(pd) - l)
    pd1 <- head(pd, length(pd) - l)

    cikl <- peak(pd1, 0, 1000, 0.7, dt)

    ## Ciklusok szama
    cn <- length(cikl)

    h <- max(diff(cikl)) # a leghosszabb ciklus
    dia<-0
    for (i in 1:cn) {
        if (cikl[i] - 200 < 0) {
            st <- 0
        } else {
            st < -cikl[i] - 200
        }
        d   <- which.min(pd1[st:cikl[i]])
        dia <- c(dia, st + d)
  }

  dia <- dia[-1]

    ## dataframe initialization
    pres_m <- data.frame(matrix(ncol = 4, nrow = cn))
    colnames(pres_m) <- c("prox_mean","dist_mean", "FFR", "iFR")

    for (i in 1:cn) {

        ## incisura detektalas
        inc.i.p <- zc(diff(pp1[dia[i]:(dia[i] + h)], 20, 4),250,0)
        inc.i.d <- zc(diff(pd1[dia[i]:(dia[i] + h)], 20, 4),250,0)

        ## ciklus-vege diastole
        pp.d.i <- which(pp1 == min(pp1[inc.i.p:h], na.rm = TRUE))
        pd.d.i <- which(pd1 == min(pd1[inc.i.d:h], na.rm = TRUE))

        pres_m[i, 1] <- mean(pp1[dia[i]:(dia[i] + h)])
        pres_m[i, 2] <- mean(pd1[dia[i]:(dia[i] + h)])
        if (i == 1) {
            pres_m[i, 3] <- NA
        } else {
            pres_m[i, 3] <- pres_m[i, 2] / pres_m[i, 1]
        }

        pres_m[i,4] <- wfp(pp1[dia[i]:(dia[i] + h)], inc.i.p, pp.d.i)

    }

    pres_m

}
