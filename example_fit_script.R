## This code was built by Choongwon Jeong

# rm(list=ls())
setwd("~/Software/github/TCLamnidis/DATES_multi_wave/")

source("exponential_fit_functions.R")

contigs = c("all", paste("", 1:22, sep=""))

mindist.val = 1.0  ## Set the minimal distance (cM) to be used

fn = "test/LOM070"
fn1 = paste(fn, ".LD.txt.gz", sep="")  ## LD table file
td1 = read.table(gzfile(fn1), header=T)  ## Import LD table file

## Split LD table file into each run (main run + 22 jackknife runs)
d1 = list()
for (i in 1:length(contigs)) {
  contig.vec = as.vector(td1$contig)
  for (j in 1:length(contigs)) d1[[j]] = td1[contig.vec == contigs[j],]
}


pdf(paste(fn, ".pdf", sep=""), height=8, width=8, colormodel="cmyk")
par(mar=c(5.1, 4.6, 2.1, 2.1))
plot_single_pulse(d1, min.dist=mindist.val, max.dist=20)
plot_single_pulse(d1, min.dist=0.5, max.dist=20)
plot_single_pulse(d1, min.dist=1.0, max.dist=20)
plot_two_pulse(d1, min.dist=0.5, max.dist=20)
plot_two_pulse(d1, min.dist=1.0, max.dist=20)
dev.off()




