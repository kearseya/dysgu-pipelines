# Uses the copyNumber package to segment read depth data. The segments are saved
# along with a heatmap, genome plots and individual chromosome plots.

library(copynumber)

args = commandArgs(trailingOnly=TRUE) # 1st working dir, 2nd reftype, 3rd plots (yes/no)

if (length(args) < 2){
  stop("requires: working directory, assmebly type, dataset name")
}

print(args)
setwd(dir = args[1])
assembly.type <- args[2]
dataset.tag <- args[3]

cyto.band.path <- sprintf("./%s_cytoBand.txt", assembly.type)
print(cyto.band.path)
cyto.band <- read.csv(cyto.band.path, sep="\t", header=F)
cyto.band <- cyto.band[cyto.band[,1] %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                            "chr20", "chr21", "chr22", "chrX", "chrY"),]

copynumbers.all.path <- sprintf("./cnpinter/%s/copyNumber_io/copynumbers.all.csv", dataset.tag)
print(copynumbers.all.path)
data <- read.csv(copynumbers.all.path, header=T, check.names=F)
#print (head(data))
cn <- colnames(data)
samples <- (cn[3:length(cn)])


if (length(args) >= 4) {
  plot_from_library <- args[4]
} else {
  plot_from_library <- "yes"
}

if (length(args) == 6) {
  gamma_value <- args[5]
  kmin_value <- args[6]
} else {
  gamma_value <- 1000
  kmin_value <- 5
}

### winsorize
custom.winsorize <- function (data, pos.unit = "bp", arms = NULL, method = "mad", 
          tau = 2.5, k = 25, gamma = 40, iter = 1, assembly = "hg38", 
          digits = 4, return.outliers = FALSE, save.res = FALSE, file.names = NULL, 
          verbose = TRUE) 
{
  stopifnot(pos.unit %in% c("bp", "kbp", "mbp"))
  stopifnot(method %in% c("mad", "pcf"))
  # if (!assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16", "mm7", 
  #                      "mm8", "mm9")) {
  #   stop("assembly must be one of hg38, hg19, hg18, hg17 or hg16", 
  #        call. = FALSE)
  # }
  stopifnot(class(data) %in% c("matrix", "data.frame", "character"))
  isfile <- class(data) == "character"
  if (!isfile) {
    stopifnot(ncol(data) >= 3)
    chrom <- data[, 1]
    pos <- data[, 2]
    nSample <- ncol(data) - 2
    sample.names <- colnames(data)[-c(1:2)]
  }
  else {
    f <- file(data, "r")
    head <- scan(f, nlines = 1, what = "character", quiet = TRUE, 
                 sep = "\t")
    if (length(head) < 3) {
      stop("Data in file must have at least 3 columns", 
           call. = FALSE)
    }
    sample.names <- head[-c(1:2)]
    nSample <- length(sample.names)
    chrom.pos <- read.table(file = data, sep = "\t", header = TRUE, 
                            colClasses = c(rep(NA, 2), rep("NULL", nSample)), 
                            as.is = TRUE)
    chrom <- chrom.pos[, 1]
    pos <- chrom.pos[, 2]
  }
  if (is.factor(chrom)) {
    chrom <- as.character(chrom)
  }
  num.chrom <- numericChrom(chrom)
  nProbe <- length(num.chrom)
  if (!is.numeric(pos)) {
    stop("input in data column 2 (posistions) must be numeric", 
         call. = FALSE)
  }
  if (is.null(arms)) {
    arms <- getArms(num.chrom, pos, pos.unit, cyto.band)
  }
  else {
    if (length(arms) != nProbe) {
      stop("'arms' must be the same length as number of rows in data", 
           call. = FALSE)
    }
  }
  num.arms <- numericArms(num.chrom, arms)
  arm.list <- unique(num.arms)
  nArm <- length(arm.list)
  if (!save.res) {
    wins.data <- matrix(nrow = 0, ncol = nSample)
    if (return.outliers) {
      wins.outliers <- matrix(nrow = 0, ncol = nSample)
    }
  }
  else {
    if (is.null(file.names)) {
      dir.res <- "Wins_res"
      if (!dir.res %in% dir()) {
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res, "/", "wins.data.txt", 
                            sep = ""), paste(dir.res, "/", "wins.outliers.txt", 
                                             sep = ""))
    }
    else {
      if (length(file.names) < 2) {
        stop("'file.names' must be of length 2", call. = FALSE)
      }
    }
  }
  for (c in 1:nArm) {
    probe.c <- which(num.arms == arm.list[c])
    wins.data.c <- matrix(nrow = length(probe.c), ncol = 0)
    if (return.outliers || save.res) {
      wins.outliers.c <- matrix(nrow = length(probe.c), 
                                ncol = 0)
    }
    if (!isfile) {
      arm.data <- data[probe.c, -c(1:2), drop = FALSE]
    }
    else {
      arm.data <- read.table(f, nrows = length(probe.c), 
                             sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", 
                                                                            nSample)))
    }
    if (any(!sapply(arm.data, is.numeric))) {
      stop("input in data columns 3 and onwards (copy numbers) must be numeric", 
           call. = FALSE)
    }
    for (i in 1:nSample) {
      y <- arm.data[, i]
      na <- is.na(y)
      use.y <- y[!na]
      ywins <- rep(NA, length(y))
      outliers <- rep(NA, length(y))
      wins <- switch(method, mad = madWins(use.y, tau = tau, 
                                           k = k, digits = digits), pcf = pcfWins(use.y, 
                                                                                  tau = tau, k = k, gamma = gamma, iter = iter, 
                                                                                  digits = digits))
      ywins[!na] <- wins$ywin
      outliers[!na] <- wins$outliers
      ywins <- round(ywins, digits = digits)
      wins.data.c <- cbind(wins.data.c, ywins)
      if (return.outliers || save.res) {
        wins.outliers.c <- cbind(wins.outliers.c, outliers)
      }
    }
    if (!save.res) {
      wins.data <- rbind(wins.data, wins.data.c)
      if (return.outliers) {
        wins.outliers <- rbind(wins.outliers, wins.outliers.c)
      }
    }
    else {
      if (c == 1) {
        wd <- file(file.names[1], "w")
        wo <- file(file.names[2], "w")
      }
      write.table(data.frame(chrom[probe.c], pos[probe.c], 
                             wins.data.c, stringsAsFactors = FALSE), file = wd, 
                  col.names = if (c == 1) 
                    c("chrom", "pos", sample.names)
                  else FALSE, row.names = FALSE, quote = FALSE, 
                  sep = "\t")
      write.table(data.frame(chrom[probe.c], pos[probe.c], 
                             wins.outliers.c, stringsAsFactors = FALSE), file = wo, 
                  col.names = if (c == 1) 
                    c("chrom", "pos", sample.names)
                  else FALSE, row.names = FALSE, quote = FALSE, 
                  sep = "\t")
    }
    if (verbose) {
      chr <- unique(chrom[probe.c])
      a <- unique(arms[probe.c])
      cat(paste("winsorize finished for chromosome arm ", 
                chr, a, sep = ""), "\n")
    }
  }
  if (isfile) {
    close(f)
  }
  if (!save.res) {
    wins.data <- data.frame(chrom, pos, wins.data, stringsAsFactors = FALSE)
    colnames(wins.data) <- c("chrom", "pos", sample.names)
    if (return.outliers) {
      wins.outliers <- data.frame(chrom, pos, wins.outliers, 
                                  stringsAsFactors = FALSE)
      colnames(wins.outliers) <- c("chrom", "pos", sample.names)
      return(list(wins.data = wins.data, wins.outliers = wins.outliers))
    }
    else {
      return(wins.data)
    }
  }
  else {
    close(wd)
    close(wo)
    cat(paste("winsorized data were saved in file", file.names[1]), 
        sep = "\n")
    cat(paste("outliers were saved in file", file.names[2]), 
        sep = "\n")
    return(invisible(NULL))
  }
}

### pcf
custom.pcf <- function (data, pos.unit = "bp", arms = NULL, Y = NULL, kmin = 5, 
                 gamma = 40, normalize = TRUE, fast = TRUE, assembly = "hg38", 
                 digits = 4, return.est = FALSE, save.res = FALSE, file.names = NULL, 
                 verbose = TRUE) 
{
  if (!pos.unit %in% c("bp", "kbp", "mbp")) {
    stop("pos.unit must be one of bp, kbp and mbp", call. = FALSE)
  }
  # if (!assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16", "mm7", 
  #                      "mm8", "mm9")) {
  #   stop("assembly must be one of hg38, hg19, hg18, hg17 or hg16", 
  #        call. = FALSE)
  # }
  isfile.data <- class(data) == "character"
  if (!isfile.data) {
    data <- pullOutContent(data, what = "wins.data")
    stopifnot(ncol(data) >= 3)
    chrom <- data[, 1]
    position <- data[, 2]
    nSample <- ncol(data) - 2
    sampleid <- colnames(data)[-c(1:2)]
  }
  else {
    f <- file(data, "r")
    head <- scan(f, nlines = 1, what = "character", quiet = TRUE, 
                 sep = "\t")
    if (length(head) < 3) {
      stop("Data in file must have at least 3 columns", 
           call. = FALSE)
    }
    sampleid <- head[-c(1:2)]
    nSample <- length(sampleid)
    chrom.pos <- read.table(file = data, sep = "\t", header = TRUE, 
                            colClasses = c(rep(NA, 2), rep("NULL", nSample)), 
                            as.is = TRUE)
    chrom <- chrom.pos[, 1]
    position <- chrom.pos[, 2]
  }
  if (is.factor(chrom)) {
    chrom <- as.character(chrom)
  }
  num.chrom <- numericChrom(chrom)
  nProbe <- length(num.chrom)
  if (!is.numeric(position)) {
    stop("input in data column 2 (posistions) must be numeric", 
         call. = FALSE)
  }
  if (is.null(arms)) {
    arms <- getArms(num.chrom, position, pos.unit, cyto.band)
  }
  else {
    stopifnot(length(arms) == nProbe)
  }
  num.arms <- numericArms(num.chrom, arms)
  arm.list <- unique(num.arms)
  nArm <- length(arm.list)
  if (!is.null(Y)) {
    stopifnot(class(Y) %in% c("matrix", "data.frame", "character"))
    isfile.Y <- class(Y) == "character"
    if (!isfile.Y) {
      ncol.Y <- ncol(Y)
      nrow.Y <- nrow(Y)
    }
    else {
      f.y <- file(Y, "r")
      ncol.Y <- length(scan(f.y, nlines = 1, what = "character", 
                            quiet = TRUE, sep = "\t"))
      nrow.Y <- nrow(read.table(file = Y, sep = "\t", header = TRUE, 
                                colClasses = c(NA, rep("NULL", ncol.Y - 1)), 
                                as.is = TRUE))
    }
    if (nrow.Y != nProbe || ncol.Y != nSample + 2) {
      stop("Input Y does not represent the same number of probes and samples as found in input data", 
           call. = FALSE)
    }
  }
  gamma0 <- gamma
  sd <- rep(1, nSample)
  if (nProbe < 1e+05 && normalize) {
    for (j in 1:nSample) {
      if (!isfile.data) {
        sample.data <- data[, j + 2]
      }
      else {
        cc <- rep("NULL", nSample + 2)
        cc[j + 2] <- "numeric"
        sample.data <- read.table(file = data, sep = "\t", 
                                  header = TRUE, colClasses = cc)[, 1]
      }
      sd[j] <- getMad(sample.data[!is.na(sample.data)], 
                      k = 25)
    }
  }
  pcf.names <- c("chrom", "pos", sampleid)
  seg.names <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", 
                 "n.probes", "mean")
  segments <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(segments) <- seg.names
  if (return.est) {
    pcf.est <- matrix(nrow = 0, ncol = nSample)
  }
  if (save.res) {
    if (is.null(file.names)) {
      dir.res <- "pcf_results"
      if (!dir.res %in% dir()) {
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res, "/", "estimates.txt", 
                            sep = ""), paste(dir.res, "/", "segments.txt", 
                                             sep = ""))
    }
    else {
      if (length(file.names) < 2) {
        stop("'file.names' must be of length 2", call. = FALSE)
      }
    }
  }
  yest <- any(return.est, save.res)
  for (c in 1:nArm) {
    probe.c <- which(num.arms == arm.list[c])
    pos.c <- position[probe.c]
    this.arm <- unique(arms[probe.c])
    this.chrom <- unique(chrom[probe.c])
    segments.c <- data.frame(matrix(nrow = 0, ncol = 7))
    colnames(segments.c) <- seg.names
    if (yest) {
      pcf.est.c <- matrix(nrow = length(probe.c), ncol = 0)
    }
    if (!isfile.data) {
      arm.data <- data[probe.c, -c(1:2), drop = FALSE]
    }
    else {
      arm.data <- read.table(f, nrows = length(probe.c), 
                             sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", 
                                                                            nSample)))
    }
    if (any(!sapply(arm.data, is.numeric))) {
      stop("input in data columns 3 and onwards (copy numbers) must be numeric", 
           call. = FALSE)
    }
    if (!is.null(Y)) {
      if (!isfile.Y) {
        arm.Y <- Y[probe.c, -c(1:2), drop = FALSE]
      }
      else {
        arm.Y <- read.table(f.y, nrows = length(probe.c), 
                            sep = "\t", colClasses = c(rep("NULL", 2), 
                                                       rep("numeric", nSample)))
      }
      if (any(!sapply(arm.Y, is.numeric))) {
        stop("input in Y columns 3 and onwards (copy numbers) must be numeric", 
             call. = FALSE)
      }
    }
    for (i in 1:nSample) {
      sample.data <- arm.data[, i]
      obs <- !is.na(sample.data)
      obs.data <- sample.data[obs]
      if (yest) {
        yhat <- rep(NA, length(probe.c))
      }
      if (length(obs.data) > 0) {
        if (nProbe >= 1e+05 && normalize) {
          sd[i] <- getMad(obs.data, k = 25)
        }
        use.gamma <- gamma0
        if (normalize) {
          use.gamma <- gamma0 * (sd[i])^2
        }
        if (use.gamma == 0 || is.na(use.gamma)) {
          if (yest) {
            res <- list(Lengde = length(obs.data), sta = 1, 
                        mean = mean(obs.data), nIntervals = 1, 
                        yhat = rep(mean(obs.data)))
          }
          else {
            res <- list(Lengde = length(obs.data), sta = 1, 
                        mean = mean(obs.data), nIntervals = 1)
          }
        }
        else {
          if (!fast || length(obs.data) < 400) {
            res <- exactPcf(y = obs.data, kmin = kmin, 
                            gamma = use.gamma, yest = yest)
          }
          else {
            res <- selectFastPcf(x = obs.data, kmin = kmin, 
                                 gamma = use.gamma, yest = yest)
          }
        }
        seg.start <- res$sta
        seg.stop <- c(seg.start[-1] - 1, length(obs.data))
        seg.npos <- res$Lengde
        seg.mean <- res$mean
        nSeg <- res$nIntervals
        if (yest) {
          yhat[obs] <- res$yhat
        }
        pos.start <- pos.c[obs][seg.start]
        pos.stop <- pos.c[obs][seg.stop]
        if (any(!obs)) {
          nn <- findNN(pos = pos.c, obs = obs)
          new.res <- handleMissing(nn = nn, pos = pos.c, 
                                   obs = obs, pos.start = pos.start, pos.stop = pos.stop, 
                                   seg.npos = seg.npos)
          pos.start <- new.res$pos.start
          pos.stop <- new.res$pos.stop
          seg.npos <- new.res$seg.npos
          if (yest) {
            yhat[!obs] <- yhat[nn]
          }
        }
      }
      else {
        warning(paste("pcf is not run for sample ", i, 
                      " on chromosome arm ", this.chrom, this.arm, 
                      " because all observations are missing. NA is returned.", 
                      sep = ""), immediate. = TRUE, call. = FALSE)
        seg.start <- 1
        seg.stop <- length(pos.c)
        pos.start <- pos.c[seg.start]
        pos.stop <- pos.c[seg.stop]
        nSeg <- 1
        seg.mean <- NA
        seg.npos <- length(pos.c)
      }
      if (!is.null(Y)) {
        seg.mean <- rep(NA, nSeg)
        sample.y <- arm.Y[, i]
        for (s in 1:nSeg) {
          seg.mean[s] <- mean(sample.y[seg.start[s]:seg.stop[s]], 
                              na.rm = TRUE)
        }
      }
      seg.mean <- round(seg.mean, digits = digits)
      seg.arm <- rep(this.arm, nSeg)
      seg.chrom <- rep(this.chrom, nSeg)
      seg <- data.frame(rep(sampleid[i], nSeg), seg.chrom, 
                        seg.arm, pos.start, pos.stop, seg.npos, seg.mean, 
                        stringsAsFactors = FALSE)
      colnames(seg) <- seg.names
      segments.c <- rbind(segments.c, seg)
      if (yest) {
        yhat <- round(yhat, digits = digits)
        pcf.est.c <- cbind(pcf.est.c, yhat)
      }
    }
    if (save.res) {
      if (c == 1) {
        w1 <- file(file.names[1], "w")
        w2 <- file(file.names[2], "w")
      }
      write.table(data.frame(chrom[probe.c], pos.c, pcf.est.c, 
                             stringsAsFactors = FALSE), file = w1, col.names = if (c == 
                                                                                   1) 
                               pcf.names
                  else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(segments.c, file = w2, col.names = if (c == 
                                                         1) 
        seg.names
        else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
    segments <- rbind(segments, segments.c)
    if (return.est) {
      pcf.est <- rbind(pcf.est, pcf.est.c)
    }
    if (verbose) {
      cat(paste("pcf finished for chromosome arm ", this.chrom, 
                this.arm, sep = ""), "\n")
    }
  }
  if (isfile.data) {
    close(f)
  }
  if (!is.null(Y)) {
    if (isfile.Y) {
      close(f.y)
    }
  }
  if (save.res) {
    close(w1)
    close(w2)
    cat(paste("pcf-estimates were saved in file", file.names[1]), 
        sep = "\n")
    cat(paste("segments were saved in file", file.names[2]), 
        sep = "\n")
  }
  if (return.est) {
    pcf.est <- data.frame(chrom, position, pcf.est, stringsAsFactors = FALSE)
    colnames(pcf.est) <- pcf.names
    return(list(estimates = pcf.est, segments = segments))
  }
  else {
    return(segments)
  }
}

environment(custom.winsorize) <- asNamespace('copynumber')
# assignInNamespace("winsorize", custom.winsorize, ns = "copynumber")
environment(custom.pcf) <- asNamespace('copynumber')
# assignInNamespace("pcf", custom.pcf, ns = "copynumber")


# Get outliers
#wins.res <- winsorize(data=data,return.outliers=TRUE,verbose=FALSE)

# wnisorize
wins <- custom.winsorize(data=data,verbose=FALSE)
wins.path <- sprintf("./cnpinter/%s/copyNumber_io/winsorized.all.csv", dataset.tag)
write.csv(wins, file=wins.path)

# Segment
# Gamma for chromothripsis paper was 600, kmin was 10
# Also game gamma 4000, kmin 10
single.seg <- custom.pcf(data=wins, gamma=gamma_value, kmin=kmin_value, verbose=FALSE,
		  return.est=TRUE)

# Multipcf
#multi.seg <- multipcf(data=wins, gamma=500, 
#		      verbose=FALSE)

# Save segmented data
# print (head(single.seg))
# print(single.seg$segments)
segs.path <- sprintf("./cnpinter/%s/copyNumber_io/segmented.copynumber.csv", dataset.tag)
write.csv(single.seg$segments, file=segs.path)

if (plot_from_library == "yes") {
  graphs.path <- sprintf("./output/%s/R/graphs/", dataset.tag)
  
  # Plot a heatmap
  heatmap.path <- sprintf("./output/%s/R/heatmap.pdf", dataset.tag)
  pdf(heatmap.path)
  # upper.lim implies lower.lim is the same but negative
  plotHeatmap(segments=single.seg, upper.lim=1, colors=c("blue", "grey96", "red"),
  	    file.name=heatmap.path)
  dev.off()
  
  ### Save genome segmented plots
  for (i in 1:length(samples)){
  
  	samp <- samples[i]
  	
  	pdf(paste(graphs.path, samp, "/", samp, "_segs.pdf", sep=""), width=12, height=3)
  
  	plotGenome(data=data,segments=single.seg,sample=i,cex=3, h=0,
  		   ylim=c(-2, 3), ylab="Relative copy number", cex.lab=1.2,
  		   cex.lab=1.2,
  		   connect=FALSE)
  	dev.off()
  
  	# Save a plot for each chromosome also
  	for (j in 1:23){
  		samp <- samples[i]
  		pdf(paste(graphs.path, samp, "/", samp, "_chr", j, ".pdf", sep=""), width=12, height=3)
  
  		plotSample(data=data, segments=single.seg ,sample=i, chrom=c(j), cex=3, h=0,
  		   ylim=c(-2, 3), ylab="Relative copy number", cex.lab=1.2,
  		   cex.lab=1.2, plot.ideo=TRUE, assembly="hg19", connect=FALSE)
  		dev.off()
  	}
  }

  }