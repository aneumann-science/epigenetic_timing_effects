library(psych)
library(corrplot)

### Load meta-analytic results
load("results/meta_results_adhd.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev
meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
# Select columns
meta_results_adhd.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig","z_birth_cat","z_school_cat_rev","z_birth_cat_abs","z_school_cat_rev_abs")]
meta_results.data <- NULL; gc()
names(meta_results_adhd.data) <- c("cpg","birth_adhd","school_adhd","birth_sig_adhd","school_sig_adhd","z_birth_cat_adhd","z_school_cat_rev_adhd","z_birth_cat_abs_adhd","z_school_cat_rev_abs_adhd")

load("results/meta_results_gpf.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev
meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
# Select columns
meta_results_gpf.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig","z_birth_cat","z_school_cat_rev","z_birth_cat_abs","z_school_cat_rev_abs")]
meta_results.data <- NULL; gc()
names(meta_results_gpf.data) <- c("cpg","birth_gpf","school_gpf","birth_sig_gpf","school_sig_gpf","z_birth_cat_gpf","z_school_cat_rev_gpf","z_birth_cat_abs_gpf","z_school_cat_rev_abs_gpf")

load("results/meta_results_sleep.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev
meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
# Select columns
meta_results_sleep.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig","z_birth_cat","z_school_cat_rev","z_birth_cat_abs","z_school_cat_rev_abs")]
meta_results.data <- NULL; gc()
names(meta_results_sleep.data) <- c("cpg","birth_sleep","school_sleep","birth_sig_sleep","school_sig_sleep","z_birth_cat_sleep","z_school_cat_rev_sleep","z_birth_cat_abs_sleep","z_school_cat_rev_abs_sleep")

load("results/meta_results_bmi.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev
meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
# Select columns
meta_results_bmi.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig","z_birth_cat","z_school_cat_rev","z_birth_cat_abs","z_school_cat_rev_abs")]
meta_results.data <- NULL; gc()
names(meta_results_bmi.data) <- c("cpg","birth_bmi","school_bmi","birth_sig_bmi","school_sig_bmi","z_birth_cat_bmi","z_school_cat_rev_bmi","z_birth_cat_abs_bmi","z_school_cat_rev_abs_bmi")

load("results/meta_results_asthma.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev
meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
# Select columns
meta_results_asthma.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig","z_birth_cat","z_school_cat_rev","z_birth_cat_abs","z_school_cat_rev_abs")]
meta_results.data <- NULL; gc()
names(meta_results_asthma.data) <- c("cpg","birth_asthma","school_asthma","birth_sig_asthma","school_sig_asthma","z_birth_cat_asthma","z_school_cat_rev_asthma","z_birth_cat_abs_asthma","z_school_cat_rev_abs_asthma")

meta_results_all.data <- Reduce(function(x, y) merge(x, y, by="cpg", all=T), list(meta_results_adhd.data,meta_results_gpf.data,meta_results_sleep.data,meta_results_bmi.data,meta_results_asthma.data))
meta_results_birth_sig.data <- meta_results_all.data[rowMeans(meta_results_all.data[c("birth_sig_adhd","birth_sig_gpf","birth_sig_sleep","birth_sig_bmi","birth_sig_asthma")], na.rm =T) > 0, ]
meta_results_school_sig.data <- meta_results_all.data[rowMeans(meta_results_all.data[c("school_sig_adhd","school_sig_gpf","school_sig_sleep","school_sig_bmi","school_sig_asthma")], na.rm =T) > 0, ]

### All CpGs
correlations_all <- corr.test(meta_results_all.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_all.mat <- correlations_all$r
rownames(correlations_all.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_all.mat) <- rownames(correlations_all.mat)

cairo_pdf(file = "figures/correlations_all.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_all.mat, title = "β correlations between time points and outcomes", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20), tl.srt = 45)
dev.off()

### Absolute effect sizes
correlations_abs <- corr.test(abs(meta_results_all.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")]), method = "spearman")
correlations_abs.mat <- correlations_abs$r
rownames(correlations_abs.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_abs.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_abs.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_abs.mat, title = "Absolute β correlations between time points and outcomes (All CpGs)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

### Z scores signed
correlations_z <- corr.test(abs(meta_results_all.data[c("z_birth_cat_adhd","z_school_cat_rev_adhd","z_birth_cat_gpf","z_school_cat_rev_gpf","z_birth_cat_sleep","z_school_cat_rev_sleep","z_birth_cat_bmi","z_school_cat_rev_sleep","z_birth_cat_asthma","z_school_cat_rev_asthma")]), method = "spearman")
correlations_z_abs.mat <- correlations_z_abs$r
rownames(correlations_z_abs.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_z_abs.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_z.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_z_abs.mat, title = "z correlations between time points and outcomes (All CpGs)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

### Z scores absolute
correlations_z_abs <- corr.test(abs(meta_results_all.data[c("z_birth_cat_abs_adhd","z_school_cat_rev_abs_adhd","z_birth_cat_abs_gpf","z_school_cat_rev_abs_gpf","z_birth_cat_abs_sleep","z_school_cat_rev_abs_sleep","z_birth_cat_abs_bmi","z_school_cat_rev_abs_sleep","z_birth_cat_abs_asthma","z_school_cat_rev_abs_asthma")]), method = "spearman")
correlations_z_abs.mat <- correlations_z_abs$r
rownames(correlations_z_abs.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_z_abs.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_z_abs.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_z_abs.mat, title = "Absolute z correlations between time points and outcomes (All CpGs)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()


correlations_birth_sig <- corr.test(meta_results_birth_sig.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_birth_sig.mat <- correlations_birth_sig$r
rownames(correlations_birth_sig.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_birth_sig.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_birth_sig.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_birth_sig.mat, title = "β correlations between time points and outcomes (CpGs nominal at birth)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

correlations_school_sig <- corr.test(meta_results_school_sig.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_school_sig.mat <- correlations_school_sig$r
rownames(correlations_school_sig.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_school_sig.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_school_sig.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_school_sig.mat, title = "β correlations between time points and outcomes (CpGs nominal at school age)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

### Only independent probes
# Load list of correlated CpGs
load("datasets/methylation_birth_mcb_r03.Rdata")

# Correlated CpGs
cor_cpgs <-methylation_birth_mcb_r03$MCBsites 

# cor_cpgs <-methylation_birth_mcb_r03$MCBsites 
mcb.list <- methylation_birth_mcb_r03$MCBinformation[,4]

# Pick per block all CpGs except randomly 1
set.seed(20231211)
mcb_random.list <- lapply(mcb.list, function(mcb) {
  cpgs <- unlist(strsplit(mcb, " "))
  sample(cpgs, length(cpgs)-1)
})

mcb_random <- unlist(mcb_random.list)

meta_results_all_ind.data <- meta_results_all.data[meta_results_all.data$cpg %in% mcb_random, ]
meta_results_birth_sig_ind.data <- meta_results_birth_sig.data[meta_results_birth_sig.data$cpg %in% mcb_random, ]
meta_results_school_sig_ind.data <- meta_results_school_sig.data[meta_results_school_sig.data$cpg %in% mcb_random, ]

correlations_all_ind <- corr.test(meta_results_all_ind.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_all_ind.mat <- correlations_all_ind$r
rownames(correlations_all_ind.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_all_ind.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_all_ind.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_all_ind.mat, title = "β correlations between time points and outcomes (All CpGs)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

correlations_birth_sig_ind <- corr.test(meta_results_birth_sig_ind.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_birth_sig_ind.mat <- correlations_birth_sig_ind$r
rownames(correlations_birth_sig_ind.mat) <- c("ADHD (Birth)","ADHD (School Age)","GPF (Birth)","GPF (School Age)","Sleep (Birth)","Sleep (School Age)","BMI (Birth)","BMI (School Age)","Asthma (Birth)","Asthma (School Age)")
colnames(correlations_birth_sig_ind.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_birth_sig_ind.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_birth_sig_ind.mat, title = "β correlations between time points and outcomes (CpGs nominal at birth)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()

correlations_school_sig_ind <- corr.test(meta_results_school_sig_ind.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi","birth_asthma","school_asthma")], method = "spearman")
correlations_school_sig_ind.mat <- correlations_school_sig_ind$r
rownames(correlations_school_sig_ind.mat) <- c("ADHD (Birth)","ADHD (Childhood)","GPF (Birth)","GPF (Childhood)","Sleep (Birth)","Sleep (Childhood)","BMI (Birth)","BMI (Childhood)","Asthma (Birth)","Asthma (Childhood)")
colnames(correlations_school_sig_ind.mat) <- c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma")

cairo_pdf(file = "figures/correlations_school_sig_ind.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_school_sig_ind.mat, title = "β correlations between time points and outcomes (CpGs nominal at school age)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20))
dev.off()


