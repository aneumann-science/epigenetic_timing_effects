library(data.table)
library(readr)
library(dplyr)
library(lattice)
library(gridExtra)
library(karyoploteR)
source("qqunif.plot.R")
library(cowplot)

# Annotation data
annotation.data <- read_csv("datasets/humanmethylation450_ANNOTATION.csv")
annotation.data <- annotation.data[c("IlmnID","CHR","MAPINFO")]

### ADHD
load("results/meta_results_adhd.Rdata")
meta_results.data <- merge(meta_results.data, annotation.data, by.x = "cpg", by.y = "IlmnID")
adhd_birth.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_birth_cat","se_birth_cat","pval_birth_cat")]
adhd_school.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev")]
adhd_change.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat")]

names(adhd_birth.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(adhd_school.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(adhd_change.data) <- c("ID","CHROM","POS","BETA","SE","P")

### GPF
load("results/meta_results_gpf.Rdata")
meta_results.data <- merge(meta_results.data, annotation.data, by.x = "cpg", by.y = "IlmnID")
gpf_birth.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_birth_cat","se_birth_cat","pval_birth_cat")]
gpf_school.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev")]
gpf_change.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat")]

names(gpf_birth.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(gpf_school.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(gpf_change.data) <- c("ID","CHROM","POS","BETA","SE","P")

### Sleep
load("results/meta_results_sleep.Rdata")
meta_results.data <- merge(meta_results.data, annotation.data, by.x = "cpg", by.y = "IlmnID")
sleep_birth.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_birth_cat","se_birth_cat","pval_birth_cat")]
sleep_school.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev")]
sleep_change.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat")]

names(sleep_birth.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(sleep_school.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(sleep_change.data) <- c("ID","CHROM","POS","BETA","SE","P")

### BMI
load("results/meta_results_bmi.Rdata")
meta_results.data <- merge(meta_results.data, annotation.data, by.x = "cpg", by.y = "IlmnID")
bmi_birth.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_birth_cat","se_birth_cat","pval_birth_cat")]
bmi_school.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev")]
bmi_change.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat")]

names(bmi_birth.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(bmi_school.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(bmi_change.data) <- c("ID","CHROM","POS","BETA","SE","P")

### Asthma
load("results/meta_results_asthma.Rdata")
meta_results.data <- merge(meta_results.data, annotation.data, by.x = "cpg", by.y = "IlmnID")
asthma_birth.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_birth_cat","se_birth_cat","pval_birth_cat")]
asthma_school.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev")]
asthma_change.data <- meta_results.data[c("cpg","CHR","MAPINFO","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat")]

names(asthma_birth.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(asthma_school.data) <- c("ID","CHROM","POS","BETA","SE","P")
names(asthma_change.data) <- c("ID","CHROM","POS","BETA","SE","P")

### Figures
birth.list <- list(adhd_birth.data,gpf_birth.data,sleep_birth.data,bmi_birth.data,asthma_birth.data)
school.list <- list(adhd_school.data,gpf_school.data,sleep_school.data,bmi_school.data,asthma_school.data)
change.list <- list(adhd_change.data,gpf_change.data,sleep_change.data,bmi_change.data,asthma_change.data)

####### Lambda
lambda <- function(summary.data){
  z=qnorm(summary.data$P/2)
  lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
}
lapply(birth.list, lambda)
lapply(school.list, lambda)
lapply(change.list, lambda)

### Convert to GRanges and winsorize high p-values
convert <- function(summary.data) {
  summary.data <- summary.data[,c("ID","CHROM","POS","P")]
  # Winsorize
  summary.data$P[summary.data$P<1E-10] <- 1E-10
  summary.data$CHROM[summary.data$CHROM==23] <- "X"
  names(summary.data) <- c("SNP","CHR","BP","p")
  summary.gr <- toGRanges(summary.data)
  seqlevelsStyle(summary.gr) <- "UCSC"
  return(summary.gr)
}

birth_gr.list <- lapply(birth.list, convert)
school_gr.list <- lapply(school.list, convert)
change_gr.list <- lapply(change.list, convert)

####### QQ-plots
birth_p.list <- list("ADHD λ = 1.60"=birth_gr.list[[1]]$p, "GPF λ = 0.94"=birth_gr.list[[2]]$p, "Sleep λ = 0.91"=birth_gr.list[[3]]$p, "BMI λ = 0.88" = birth_gr.list[[4]]$p, "Asthma λ = 1.05" = birth_gr.list[[5]]$p)
school_p.list <- list("ADHD λ = 0.98"=school_gr.list[[1]]$p, "GPF λ = 0.95"=school_gr.list[[2]]$p, "Sleep λ = 0.93"=school_gr.list[[3]]$p, "BMI λ = 1.16" = school_gr.list[[4]]$p, "Asthma λ = 0.91" = school_gr.list[[5]]$p)
change_p.list <- list("ADHD λ = 0.90"=change_gr.list[[1]]$p, "GPF λ = 1.02"=change_gr.list[[2]]$p, "Sleep λ = 0.96"=change_gr.list[[3]]$p, "BMI λ = 1.13" = change_gr.list[[4]]$p, "Asthma λ = 1.10" = change_gr.list[[5]]$p)

birth <- qqunif.plot(birth_p.list, auto.key=list(corner=c(.95,.05)), main = "DNAm at Birth", par.settings = list(superpose.symbol=list(pch=20),fontsize=list(text=40,points=30)), xlim = 0:7, ylim = 0:10)
school <- qqunif.plot(school_p.list, auto.key=list(corner=c(.95,.05)), main = "DNAm at Childhood", par.settings = list(superpose.symbol=list(pch=20),fontsize=list(text=40,points=30)), xlim = 0:7, ylim = 0:10)
change <- qqunif.plot(change_p.list, auto.key=list(corner=c(.95,.05)), main = "Change from Birth to Childhood", par.settings = list(superpose.symbol=list(pch=20),fontsize=list(text=40,points=30)), xlim = 0:7, ylim = 0:10)

qq.list <- list(birth,school,change)

png(file="figures/qq_meta.png", width = 3000, height = 1500)
grid.arrange(grobs = qq.list, nrow = 1)
dev.off()

####### Manhattan
bonferroni_adhd <- 0.05/430327
bonferroni_gpf <- 0.05/372292
bonferroni_sleep <- 0.05/431159
bonferroni_bmi <- 0.05/435652
bonferroni_asthma <- 0.05/432728

png(file="figures/manhattan_birth.png", width = 4000, height = 4000, pointsize = 30)
kp_birth <- plotKaryotype(plot.type=4, chromosomes = c(paste0("chr",1:22)), cex = 0.9)
kpAddLabels(kp_birth, labels = "Asthma", srt=90, pos=3, r0=autotrack(1,5), cex = 2)
kpAddLabels(kp_birth, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(1,5))
kpAxis(kp_birth, ymin = 0, ymax=10, numticks = 11, r0=autotrack(1,5), side = 2)
kp_birth <- kpPlotManhattan(kp_birth, data=birth_gr.list[[5]], genomewideline = -log10(bonferroni_adhd), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(1,5))
kpAddLabels(kp_birth, labels = "BMI", srt=90, pos=3, r0=autotrack(2,5), cex = 2)
kpAddLabels(kp_birth, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(2,5))
kpAxis(kp_birth, ymin = 0, ymax=10, numticks = 11, r0=autotrack(2,5), side = 2)
kp_birth <- kpPlotManhattan(kp_birth, data=birth_gr.list[[4]], genomewideline = -log10(bonferroni_gpf), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(2,5))
kpAddLabels(kp_birth, labels = "Sleep", srt=90, pos=3, r0=autotrack(3,5), cex = 2)
kpAddLabels(kp_birth, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(3,5))
kpAxis(kp_birth, ymin = 0, ymax=10, numticks = 11, r0=autotrack(3,5), side = 2)
kp_birth <- kpPlotManhattan(kp_birth, data=birth_gr.list[[3]], genomewideline = -log10(bonferroni_sleep), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(3,5))
kpAddLabels(kp_birth, labels = "GPF", srt=90, pos=3, r0=autotrack(4,5), cex = 2)
kpAddLabels(kp_birth, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(4,5))
kpAxis(kp_birth, ymin = 0, ymax=10, numticks = 11, r0=autotrack(4,5), side = 2)
kp_birth <- kpPlotManhattan(kp_birth, data=birth_gr.list[[2]], genomewideline = -log10(bonferroni_bmi), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(4,5))
kpAddLabels(kp_birth, labels = "ADHD", srt=90, pos=3, r0=autotrack(5,5), cex = 2)
kpAddLabels(kp_birth, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(5,5))
kpAxis(kp_birth, ymin = 0, ymax=10, numticks = 11, r0=autotrack(5,5), side = 2)
kpAddMainTitle(kp_birth, main="DNAm-health associations at birth", cex = 4)
kp_birth <- kpPlotManhattan(kp_birth, data=birth_gr.list[[1]], genomewideline = -log10(bonferroni_asthma), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(5,5))
dev.off()


png(file="figures/manhattan_childhood.png", width = 4000, height = 4000, pointsize = 30)
kp_school <- plotKaryotype(plot.type=4, chromosomes = c(paste0("chr",1:22)), cex = 0.9)
kpAddLabels(kp_school, labels = "Asthma", srt=90, pos=3, r0=autotrack(1,5), cex = 2)
kpAddLabels(kp_school, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(1,5))
kpAxis(kp_school, ymin = 0, ymax=10, numticks = 11, r0=autotrack(1,5), side = 2)
kp_school <- kpPlotManhattan(kp_school, data=school_gr.list[[5]], genomewideline = -log10(bonferroni_adhd), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(1,5))
kpAddLabels(kp_school, labels = "BMI", srt=90, pos=3, r0=autotrack(2,5), cex = 2)
kpAddLabels(kp_school, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(2,5))
kpAxis(kp_school, ymin = 0, ymax=10, numticks = 11, r0=autotrack(2,5), side = 2)
kp_school <- kpPlotManhattan(kp_school, data=school_gr.list[[4]], genomewideline = -log10(bonferroni_gpf), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(2,5))
kpAddLabels(kp_school, labels = "Sleep", srt=90, pos=3, r0=autotrack(3,5), cex = 2)
kpAddLabels(kp_school, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(3,5))
kpAxis(kp_school, ymin = 0, ymax=10, numticks = 11, r0=autotrack(3,5), side = 2)
kp_school <- kpPlotManhattan(kp_school, data=school_gr.list[[3]], genomewideline = -log10(bonferroni_sleep), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(3,5))
kpAddLabels(kp_school, labels = "GPF", srt=90, pos=3, r0=autotrack(4,5), cex = 2)
kpAddLabels(kp_school, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(4,5))
kpAxis(kp_school, ymin = 0, ymax=10, numticks = 11, r0=autotrack(4,5), side = 2)
kp_school <- kpPlotManhattan(kp_school, data=school_gr.list[[2]], genomewideline = -log10(bonferroni_bmi), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(4,5))
kpAddLabels(kp_school, labels = "ADHD", srt=90, pos=3, r0=autotrack(5,5), cex = 2)
kpAddLabels(kp_school, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(5,5))
kpAxis(kp_school, ymin = 0, ymax=10, numticks = 11, r0=autotrack(5,5), side = 2)
kpAddMainTitle(kp_school, main="DNAm-health associations in childhood", cex = 4)
kp_school <- kpPlotManhattan(kp_school, data=school_gr.list[[1]], genomewideline = -log10(bonferroni_asthma), suggestive.col="orange", suggestive.lwd = 3,
                            genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(5,5))
dev.off()

png(file="figures/manhattan_change.png", width = 4000, height = 4000, pointsize = 30)
kp_change <- plotKaryotype(plot.type=4, chromosomes = c(paste0("chr",1:22)), cex = 0.9)
kpAddLabels(kp_change, labels = "Asthma", srt=90, pos=3, r0=autotrack(1,5), cex = 2)
kpAddLabels(kp_change, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(1,5))
kpAxis(kp_change, ymin = 0, ymax=10, numticks = 11, r0=autotrack(1,5), side = 2)
kp_change <- kpPlotManhattan(kp_change, data=change_gr.list[[5]], genomewideline = -log10(bonferroni_adhd), suggestive.col="orange", suggestive.lwd = 3,
                             genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(1,5))
kpAddLabels(kp_change, labels = "BMI", srt=90, pos=3, r0=autotrack(2,5), cex = 2)
kpAddLabels(kp_change, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(2,5))
kpAxis(kp_change, ymin = 0, ymax=10, numticks = 11, r0=autotrack(2,5), side = 2)
kp_change <- kpPlotManhattan(kp_change, data=change_gr.list[[4]], genomewideline = -log10(bonferroni_gpf), suggestive.col="orange", suggestive.lwd = 3,
                             genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(2,5))
kpAddLabels(kp_change, labels = "Sleep", srt=90, pos=3, r0=autotrack(3,5), cex = 2)
kpAddLabels(kp_change, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(3,5))
kpAxis(kp_change, ymin = 0, ymax=10, numticks = 11, r0=autotrack(3,5), side = 2)
kp_change <- kpPlotManhattan(kp_change, data=change_gr.list[[3]], genomewideline = -log10(bonferroni_sleep), suggestive.col="orange", suggestive.lwd = 3,
                             genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(3,5))
kpAddLabels(kp_change, labels = "GPF", srt=90, pos=3, r0=autotrack(4,5), cex = 2)
kpAddLabels(kp_change, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(4,5))
kpAxis(kp_change, ymin = 0, ymax=10, numticks = 11, r0=autotrack(4,5), side = 2)
kp_change <- kpPlotManhattan(kp_change, data=change_gr.list[[2]], genomewideline = -log10(bonferroni_bmi), suggestive.col="orange", suggestive.lwd = 3,
                             genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(4,5))
kpAddLabels(kp_change, labels = "ADHD", srt=90, pos=3, r0=autotrack(5,5), cex = 2)
kpAddLabels(kp_change, labels = "-log⏨(p)", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(5,5))
kpAxis(kp_change, ymin = 0, ymax=10, numticks = 11, r0=autotrack(5,5), side = 2)
kpAddMainTitle(kp_change, main="Change in DNAm-health associations from birth to childhood", cex = 4)
kp_change <- kpPlotManhattan(kp_change, data=school_gr.list[[1]], genomewideline = -log10(bonferroni_asthma), suggestive.col="orange", suggestive.lwd = 3,
                             genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(5,5))
dev.off()


system("convert +append figures/manhattan_birth.png figures/manhattan_childhood.png figures/manhattan.png")

# # Miami plot
# png(file="figures/manhattan.png", width = 4000, height = 6000, pointsize = 30)
# kp <- plotKaryotype(plot.type=4, chromosomes = c(paste0("chr",1:22),"chrX"))
# kpAddLabels(kp, labels = "ADHD", srt=90, pos=3, r0=autotrack(c(1,2),10))
# kpAddLabels(kp, labels = "-log⏨(p) birth", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(10,10,r=1,r1=0))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(10,10,r=1,r1=0), side = 2)
# kp <- kpPlotManhattan(kp, data=birth_gr.list[[1]], genomewideline = -log10(bonferroni_adhd), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(10,10,r=1,r1=0))
# kpAddLabels(kp, labels = "-log⏨(p) school age", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(2,10,r=0,r1=1))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(2,10,r=0,r1=1), side = 2)
# kp <- kpPlotManhattan(kp, data=school_gr.list[[1]], genomewideline = -log10(bonferroni_adhd), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#E1BE6A",chr2="#40B0A6",chr3="#E1BE6A",chr4="#40B0A6",chr5="#E1BE6A",chr6="#40B0A6",chr7="#E1BE6A",chr8="#40B0A6",chr9="#E1BE6A",chr10="#40B0A6",chr11="#E1BE6A",chr12="#40B0A6",chr13="#E1BE6A",chr14="#40B0A6",chr15="#E1BE6A",chr16="#40B0A6",chr17="#E1BE6A",chr18="#40B0A6",chr19="#E1BE6A",chr20="#40B0A6",chr21="#E1BE6A",chr22="#40B0A6",chrX="#E1BE6A"), r0=autotrack(2,10,r=0,r1=1))
# kpAddLabels(kp, labels = "GPF", srt=90, pos=3, r0=autotrack(c(3,4),10))
# kpAddLabels(kp, labels = "-log⏨(p) birth", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(8,10,r=1,r1=0))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(8,10,r=1,r1=0), side = 2)
# kp <- kpPlotManhattan(kp, data=birth_gr.list[[2]], genomewideline = -log10(bonferroni_gpf), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(8,10,r=1,r1=0))
# kpAddLabels(kp, labels = "-log⏨(p) school age", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(4,10,r=0,r1=1))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(4,10,r=0,r1=1), side = 2)
# kp <- kpPlotManhattan(kp, data=school_gr.list[[2]], genomewideline = -log10(bonferroni_gpf), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#E1BE6A",chr2="#40B0A6",chr3="#E1BE6A",chr4="#40B0A6",chr5="#E1BE6A",chr6="#40B0A6",chr7="#E1BE6A",chr8="#40B0A6",chr9="#E1BE6A",chr10="#40B0A6",chr11="#E1BE6A",chr12="#40B0A6",chr13="#E1BE6A",chr14="#40B0A6",chr15="#E1BE6A",chr16="#40B0A6",chr17="#E1BE6A",chr18="#40B0A6",chr19="#E1BE6A",chr20="#40B0A6",chr21="#E1BE6A",chr22="#40B0A6",chrX="#E1BE6A"), r0=autotrack(4,10,r=0,r1=1))
# kpAddLabels(kp, labels = "Sleep", srt=90, pos=3, r0=autotrack(c(5,6),10))
# kpAddLabels(kp, labels = "-log⏨(p) birth", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(6,10,r=1,r1=0))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(6,10,r=1,r1=0), side = 2)
# kp <- kpPlotManhattan(kp, data=birth_gr.list[[3]], genomewideline = -log10(bonferroni_sleep), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(6,10,r=1,r1=0))
# kpAddLabels(kp, labels = "-log⏨(p) school age", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(6,10,r=0,r1=1))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(6,10,r=0,r1=1), side = 2)
# kp <- kpPlotManhattan(kp, data=school_gr.list[[3]], genomewideline = -log10(bonferroni_sleep), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#E1BE6A",chr2="#40B0A6",chr3="#E1BE6A",chr4="#40B0A6",chr5="#E1BE6A",chr6="#40B0A6",chr7="#E1BE6A",chr8="#40B0A6",chr9="#E1BE6A",chr10="#40B0A6",chr11="#E1BE6A",chr12="#40B0A6",chr13="#E1BE6A",chr14="#40B0A6",chr15="#E1BE6A",chr16="#40B0A6",chr17="#E1BE6A",chr18="#40B0A6",chr19="#E1BE6A",chr20="#40B0A6",chr21="#E1BE6A",chr22="#40B0A6",chrX="#E1BE6A"), r0=autotrack(6,10,r=0,r1=1))
# kpAddLabels(kp, labels = "BMI", srt=90, pos=3, r0=autotrack(c(7,8),10))
# kpAddLabels(kp, labels = "-log⏨(p) birth", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(4,10,r=1,r1=0))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(4,10,r=1,r1=0), side = 2)
# kp <- kpPlotManhattan(kp, data=birth_gr.list[[4]], genomewideline = -log10(bonferroni_bmi), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(4,10,r=1,r1=0))
# kpAddLabels(kp, labels = "-log⏨(p) school age", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(8,10,r=0,r1=1))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(8,10,r=0,r1=1), side = 2)
# kp <- kpPlotManhattan(kp, data=school_gr.list[[4]], genomewideline = -log10(bonferroni_bmi), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#E1BE6A",chr2="#40B0A6",chr3="#E1BE6A",chr4="#40B0A6",chr5="#E1BE6A",chr6="#40B0A6",chr7="#E1BE6A",chr8="#40B0A6",chr9="#E1BE6A",chr10="#40B0A6",chr11="#E1BE6A",chr12="#40B0A6",chr13="#E1BE6A",chr14="#40B0A6",chr15="#E1BE6A",chr16="#40B0A6",chr17="#E1BE6A",chr18="#40B0A6",chr19="#E1BE6A",chr20="#40B0A6",chr21="#E1BE6A",chr22="#40B0A6",chrX="#E1BE6A"), r0=autotrack(8,10,r=0,r1=1))
# kpAddLabels(kp, labels = "Asthma", srt=90, pos=3, r0=autotrack(c(9,10),10))
# kpAddLabels(kp, labels = "-log⏨(p) birth", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(2,10,r=1,r1=0))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(2,10,r=1,r1=0), side = 2)
# kp <- kpPlotManhattan(kp, data=birth_gr.list[[5]], genomewideline = -log10(bonferroni_asthma), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#FFC20A",chr2="#0C7BDC",chr3="#FFC20A",chr4="#0C7BDC",chr5="#FFC20A",chr6="#0C7BDC",chr7="#FFC20A",chr8="#0C7BDC",chr9="#FFC20A",chr10="#0C7BDC",chr11="#FFC20A",chr12="#0C7BDC",chr13="#FFC20A",chr14="#0C7BDC",chr15="#FFC20A",chr16="#0C7BDC",chr17="#FFC20A",chr18="#0C7BDC",chr19="#FFC20A",chr20="#0C7BDC",chr21="#FFC20A",chr22="#0C7BDC",chrX="#FFC20A"), r0=autotrack(2,10,r=1,r1=0))
# kpAddLabels(kp, labels = "-log⏨(p) school age", srt=90, pos=3, side=2, label.margin=0.03, r0=autotrack(10,10,r=0,r1=1))
# kpAxis(kp, ymin = 0, ymax=10, numticks = 11, r0=autotrack(10,10,r=0,r1=1), side = 2)
# kp <- kpPlotManhattan(kp, data=school_gr.list[[5]], genomewideline = -log10(bonferroni_asthma), suggestive.col="orange", suggestive.lwd = 3,
#                       genomewide.col = "red", genomewide.lwd = 6, ymax = 10, points.cex = 0.8, points.col = c(chr1="#E1BE6A",chr2="#40B0A6",chr3="#E1BE6A",chr4="#40B0A6",chr5="#E1BE6A",chr6="#40B0A6",chr7="#E1BE6A",chr8="#40B0A6",chr9="#E1BE6A",chr10="#40B0A6",chr11="#E1BE6A",chr12="#40B0A6",chr13="#E1BE6A",chr14="#40B0A6",chr15="#E1BE6A",chr16="#40B0A6",chr17="#E1BE6A",chr18="#40B0A6",chr19="#E1BE6A",chr20="#40B0A6",chr21="#E1BE6A",chr22="#40B0A6",chrX="#E1BE6A"), r0=autotrack(10,10,r=0,r1=1))
# dev.off()

