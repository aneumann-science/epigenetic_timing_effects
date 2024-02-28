library(psych)
library(corrplot)

### Load meta-analytic results
load("results/meta_results_adhd_alspac.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
meta_results_adhd.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig")]
meta_results.data <- NULL; gc()
names(meta_results_adhd.data) <- c("cpg","birth_adhd","school_adhd","birth_sig_adhd","school_sig_adhd")

load("results/meta_results_gpf_alspac.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
meta_results_gpf.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig")]
meta_results.data <- NULL; gc()
names(meta_results_gpf.data) <- c("cpg","birth_gpf","school_gpf","birth_sig_gpf","school_sig_gpf")

load("results/meta_results_sleep_alspac.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
meta_results_sleep.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig")]
meta_results.data <- NULL; gc()
names(meta_results_sleep.data) <- c("cpg","birth_sleep","school_sleep","birth_sig_sleep","school_sig_sleep")

load("results/meta_results_bmi_alspac.Rdata")
meta_results.data[meta_results.data$pval_birth < 1E-04, "birth_sig"] <- 1
meta_results.data$birth_sig[is.na(meta_results.data$birth_sig)] <- 0
meta_results.data[meta_results.data$pval_school_cat_rev < 1E-04, "school_sig"] <- 1
meta_results.data$school_sig[is.na(meta_results.data$school_sig)] <- 0
meta_results_bmi.data <- meta_results.data[, c("cpg","estimate_birth","estimate_school_cat_rev","birth_sig","school_sig")]
meta_results.data <- NULL; gc()
names(meta_results_bmi.data) <- c("cpg","birth_bmi","school_bmi","birth_sig_bmi","school_sig_bmi")

meta_results_all.data <- Reduce(function(x, y) merge(x, y, by="cpg", all=T), list(meta_results_adhd.data,meta_results_gpf.data,meta_results_sleep.data,meta_results_bmi.data))
meta_results_birth_sig.data <- meta_results_all.data[rowMeans(meta_results_all.data[c("birth_sig_adhd","birth_sig_gpf","birth_sig_sleep","birth_sig_bmi")], na.rm =T) > 0, ]
meta_results_school_sig.data <- meta_results_all.data[rowMeans(meta_results_all.data[c("school_sig_adhd","school_sig_gpf","school_sig_sleep","school_sig_bmi")], na.rm =T) > 0, ]

correlations_all <- corr.test(meta_results_all.data[c("birth_adhd","school_adhd","birth_gpf","school_gpf","birth_sleep","school_sleep","birth_bmi","school_bmi")], method = "spearman")
correlations_all.mat <- correlations_all$r
rownames(correlations_all.mat) <- c("ADHD (Birth)","ADHD (Childhood)","GPF (Birth)","GPF (Childhood)","Sleep (Birth)","Sleep (Childhood)","BMI (Birth)","BMI (Childhood)")
colnames(correlations_all.mat) <- rownames(correlations_all.mat)

cairo_pdf(file = "figures/correlations_alspac.pdf", width = 40, height = 30, pointsize = 50)
corrplot(correlations_all.mat, title = "β correlations between time points and outcomes (ALSPAC only)", method = "color", addCoef.col = "black", mar=c(0,0,2,0),
         col = colorRampPalette(c("blue","white", "red"))(20), tl.srt = 45)
dev.off()
