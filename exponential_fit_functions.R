## This code was built by Choongwon Jeong

## A function to do a single-pulse exponential fit
single_pulse = function(dmat, mindist) {
  dist.vec = as.vector(dmat$dist.cM)
  ld.vec = as.vector(dmat$weighted.LD)
  fil.vec = (dist.vec >= mindist) & (dist.vec < Inf)
  affine.term = ld.vec[length(ld.vec)]

  dist.vec.2 = dist.vec[fil.vec]; ld.vec.2 = ld.vec[fil.vec]
  model.dat = nls(ld.vec.2 ~ exp(a -0.01 * b * dist.vec.2) + affine.term, start = list(a=-10, b=10), control=list(maxiter=5000))
  return(model.dat)
}


## Single-pulse fitting wrap-up script for jackknifing
single_pulse_wrap = function(dlist, mindist) {
  tms = list(); tm.amp = c(); tm.decay = c()
  n.contig = length(dlist)
  for (i in 1:n.contig) {
    tms[[i]] = single_pulse(dlist[[i]], mindist)
    tm.amp = c(tm.amp, exp(summary(tms[[i]])$coefficient[1,1]))
    tm.decay = c(tm.decay, summary(tms[[i]])$coefficient[2,1])
  }
  tm.amp.sd = sqrt(sum((tm.amp[2:n.contig] - tm.amp[1])**2) * (n.contig-2)/(n.contig-1))
  tm.decay.sd = sqrt(sum((tm.decay[2:n.contig] - tm.decay[1])**2) * (n.contig-2)/(n.contig-1))
  tout = list()
  tout[[1]] = tms; tout[[2]] = tm.amp; tout[[3]] = tm.decay
  tout[[4]] = c(tm.amp[1], tm.amp.sd, tm.decay[1], tm.decay.sd)
  names(tout) = c("model", "amp", "decay", "est")
  return(tout)
}

## A function to do a double-pulse exponential fit
two_pulse = function(dmat, mindist) {
  dist.vec = as.vector(dmat$dist.cM)
  ld.vec = as.vector(dmat$weighted.LD)
  fil.vec = (dist.vec >= mindist) & (dist.vec < Inf)
  affine.term = ld.vec[length(ld.vec)]

  dist.vec.2 = dist.vec[fil.vec]; ld.vec.2 = ld.vec[fil.vec]
  model.dat = nls(ld.vec.2 ~ exp(a -0.01 * b * dist.vec.2) + exp(c -0.01 * d * dist.vec.2) + affine.term, start = list(a=-10, b=10, c=-10, d=100), control=list(maxiter=5000))
  return(model.dat)
}

## Single-pulse fitting wrap-up script for jackknifing
two_pulse_wrap = function(dlist, mindist) {
  tms = list(); tm.amp1 = c(); tm.decay1 = c(); tm.amp2 = c(); tm.decay2 = c()
  n.contig = length(dlist)
  for (i in 1:n.contig) {
    tms[[i]] = two_pulse(dlist[[i]], mindist)
    tm.amp1 = c(tm.amp1, exp(summary(tms[[i]])$coefficient[1,1]))
    tm.decay1 = c(tm.decay1, summary(tms[[i]])$coefficient[2,1])
    tm.amp2 = c(tm.amp2, exp(summary(tms[[i]])$coefficient[3,1]))
    tm.decay2 = c(tm.decay2, summary(tms[[i]])$coefficient[4,1])
  }
  tm.amp1.sd = sqrt(sum((tm.amp1[2:n.contig] - tm.amp1[1])**2) * (n.contig-2)/(n.contig-1))
  tm.decay1.sd = sqrt(sum((tm.decay1[2:n.contig] - tm.decay1[1])**2) * (n.contig-2)/(n.contig-1))
  tm.amp2.sd = sqrt(sum((tm.amp2[2:n.contig] - tm.amp2[1])**2) * (n.contig-2)/(n.contig-1))
  tm.decay2.sd = sqrt(sum((tm.decay2[2:n.contig] - tm.decay2[1])**2) * (n.contig-2)/(n.contig-1))
  tout = list()
  tout[[1]] = tms; tout[[2]] = tm.amp1; tout[[3]] = tm.decay1; tout[[4]] = tm.amp2; tout[[5]] = tm.decay2
  tout[[6]] = c(tm.amp1[1], tm.amp1.sd, tm.decay1[1], tm.decay1.sd)
  tout[[6]] = c(tout[[6]], tm.amp2[1], tm.amp2.sd, tm.decay2[1], tm.decay2.sd)
  names(tout) = c("model", "amp1", "decay1", "amp2", "decay2", "est")
  return(tout)
}

## A function for single-pulse plotting
plot_single_pulse = function(dlist, min.dist, max.dist) {
  dmat = dlist[[1]]
  dist.vec = as.vector(dmat$dist.cM)
  filvec = (dist.vec >= min.dist & dist.vec <= max.dist)

  y1 = min(dmat$weighted.LD[filvec]); y2 = max(dmat$weighted.LD[filvec])
  y.interval = 0.5*10**(-1*ceiling(log10(1/y2)))
  ymin = min(c(0,floor(y1/y.interval) * y.interval))
  ymax = ceiling(y2/y.interval) * y.interval
  yticks = seq(ymin, ymax, y.interval)

  par(mar=c(5.1, 5.1, 2.1, 2.1)); par(cex.axis=1.3, cex.lab=1.3)
  xlabvec = "Distance (cM)"
  ylabvec = expression(paste("Weighted LD (x 10"^-4, ")"))
  plot(dmat$dist.cM[filvec], dmat$weighted.LD[filvec], xlim=c(0,max.dist), ylim=c(ymin,ymax), xlab=xlabvec, ylab="", yaxt="n", pch=4, col="skyblue")
  axis(2, at=yticks, label=yticks*10**4, las=2, cex.axis=1.1)
  mtext(2, text=ylabvec, line=3.5, cex=1.3)

  dv = as.vector(dmat$dist.cM)
  dv2 = dv[dv >= min.dist & dv < Inf]
  m1 = single_pulse_wrap(dlist, mindist=min.dist)
  lines(dv2, predict(m1$model[[1]], list(x=dv2)), col="red", lwd=3)
  l1 = round(m1$est[3], digits=2); l1.sd = round(m1$est[4], digits=2)
  lg1 = paste("1 pulse (", l1, " \u00b1 ", l1.sd, "); dist=", min.dist, " cM", sep="")
  legend("topright", legend=lg1, lty=1, lwd=2, col="red", cex=1.3, y.intersp=1.3)
}

## A function for two-pulse plotting
plot_two_pulse = function(dlist, min.dist, max.dist) {
  dmat = dlist[[1]]
  dist.vec = as.vector(dmat$dist.cM)
  filvec = (dist.vec >= min.dist & dist.vec <= max.dist)

  y1 = min(dmat$weighted.LD[filvec]); y2 = max(dmat$weighted.LD[filvec])
  y.interval = 0.5*10**(-1*ceiling(log10(1/y2)))
  ymin = min(c(0,floor(y1/y.interval) * y.interval))
  ymax = ceiling(y2/y.interval) * y.interval
  yticks = seq(ymin, ymax, y.interval)

  par(mar=c(5.1, 5.1, 2.1, 2.1)); par(cex.axis=1.3, cex.lab=1.3)
  xlabvec = "Distance (cM)"
  ylabvec = expression(paste("Weighted LD (x 10"^-4, ")"))
  plot(dmat$dist.cM[filvec], dmat$weighted.LD[filvec], xlim=c(0,max.dist), ylim=c(ymin,ymax), xlab=xlabvec, ylab="", yaxt="n", pch=4, col="skyblue")
  axis(2, at=yticks, label=yticks*10**4, las=2, cex.axis=1.1)
  mtext(2, text=ylabvec, line=3.5, cex=1.3)

  dv = as.vector(dmat$dist.cM)
  dv2 = dv[dv >= min.dist & dv < Inf]
  m1 = single_pulse_wrap(dlist, mindist=min.dist)
  m2 = two_pulse_wrap(dlist, mindist=min.dist)
  lines(dv2, predict(m1$model[[1]], list(x=dv2)), col="red", lwd=3)
  lines(dv2, predict(m2$model[[1]], list(x=dv2)), col="orange", lwd=3)
  l1 = round(m1$est[3], digits=2); l1.sd = round(m1$est[4], digits=2)
  l2 = round(m2$est[3], digits=2); l2.sd = round(m2$est[4], digits=2)
  l3 = round(m2$est[7], digits=2); l3.sd = round(m2$est[8], digits=2)
  lg1 = paste("1 pulse (", l1, " \u00b1 ", l1.sd, ")", sep="")
  lg2 = paste("2 pulses (", l2, " \u00b1 ", l2.sd, "; ", l3, " \u00b1 ", l3.sd, ")", sep="")
  legend("topright", legend=c(lg1, lg2), lty=1, lwd=2, col=c("red","orange"), cex=1.3, y.intersp=1.2)
}
