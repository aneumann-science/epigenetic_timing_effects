library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)
library(ggplot2)

# Load summary data from all cohorts
adhd_genr_birth.data <- fread("sumstats/ADHD/adhd_ewas_genr_20180221.txt")
adhd_alspac_birth.data <- fread("sumstats/ADHD/adhd_alspac_7ct_plate_noPCs_trimmed_20161010.txt")
adhd_inma_birth.data <- fread("sumstats/ADHD/adhd_ewas_inma_20161021.txt")
adhd_predo_birth.data <- fread("sumstats/ADHD/ChildADHD_PREDO_10012017_trimmed_df.txt")

# Rename columns
setnames(adhd_genr_birth.data, c("beta", "se", "t", "n"), c("beta_genr_birth", "se_genr_birth", "t_genr_birth", "n_genr_birth"))
setnames(adhd_alspac_birth.data, c("beta", "se", "t", "n"), c("beta_alspac_birth", "se_alspac_birth", "t_alspac_birth", "n_alspac_birth"))
setnames(adhd_inma_birth.data, c("beta", "se", "t", "n"), c("beta_inma_birth", "se_inma_birth", "t_inma_birth", "n_inma_birth"))
setnames(adhd_predo_birth.data, c("BETA", "SE", "T_VAL", "N"), c("beta_predo_birth", "se_predo_birth", "t_predo_birth", "n_predo_birth"))

# Calculate approximate p-values
adhd_genr_birth.data$P_VAL_genr_birth <- 2 * (1 - pnorm(abs(adhd_genr_birth.data$t_genr_birth)))
adhd_alspac_birth.data$P_VAL_alspac_birth <- 2 * (1 - pnorm(abs(adhd_alspac_birth.data$t_alspac_birth)))
adhd_inma_birth.data$P_VAL_inma_birth <- 2 * (1 - pnorm(abs(adhd_inma_birth.data$t_inma_birth)))
adhd_predo_birth.data$P_VAL_predo_birth <- 2 * (1 - pnorm(abs(adhd_predo_birth.data$t_predo_birth)))

# Merge all cohorts
adhd.data <- Reduce(function(x, y) merge(x, y, by="probeID", all=TRUE), list(adhd_genr_birth.data, adhd_alspac_birth.data, adhd_inma_birth.data, adhd_predo_birth.data))

# Load summary data from all cohorts
adhd_genr_9y.data <- fread("sumstats/ADHD/adhd_ewas_genr_9y_20182102.txt")
adhd_alspac_9y.data <- fread("sumstats/ADHD/adhd_alspac_6ct_plate_noPCs_tr_20180117.txt")
adhd_helix_pak_9y.data <- fread("sumstats/ADHD/adhd_ewas_ALL_HELIX_SLIDE_PAKI_20180504.txt")
adhd_helix_all_eu_9y.data <- fread("sumstats/ADHD/adhd_ewas_ALL_HELIX_SLIDE_EU_20180427.txt")
adhd_glaku_9y.data <- fread("sumstats/ADHD/ADHD_EWAS_GLAKU_12y_20180516.txt")

# Rename columns
setnames(adhd_genr_9y.data, c("beta", "se", "t", "n"), c("beta_genr_9y", "se_genr_9y", "t_genr_9y", "n_genr_9y"))
setnames(adhd_alspac_9y.data, c("beta", "se", "t", "n"), c("beta_alspac_9y", "se_alspac_9y", "t_alspac_9y", "n_alspac_9y"))
setnames(adhd_helix_pak_9y.data, c("beta", "se", "t", "n"), c("beta_helix_pak_9y", "se_helix_pak_9y", "t_helix_pak_9y", "n_helix_pak_9y"))
setnames(adhd_helix_all_eu_9y.data, c("beta", "se", "t", "n"), c("beta_helix_all_eu_9y", "se_helix_all_eu_9y", "t_helix_all_eu_9y", "n_helix_all_eu_9y"))
setnames(adhd_glaku_9y.data, c("beta", "se", "t", "n"), c("beta_glaku_9y", "se_glaku_9y", "t_glaku_9y", "n_glaku_9y"))

# Calculate approximate p-values
adhd_genr_9y.data$P_VAL_genr_9y <- 2 * (1 - pnorm(abs(adhd_genr_9y.data$t_genr_9y)))
adhd_alspac_9y.data$P_VAL_alspac_9y <- 2 * (1 - pnorm(abs(adhd_alspac_9y.data$t_alspac_9y)))
adhd_helix_all_eu_9y.data$P_VAL_helix_all_eu_9y <- 2 * (1 - pnorm(abs(adhd_helix_all_eu_9y.data$t_helix_all_eu_9y)))
adhd_helix_pak_9y.data$P_VAL_helix_pak_9y <- 2 * (1 - pnorm(abs(adhd_helix_pak_9y.data$t_helix_pak_9y)))
adhd_glaku_9y.data$P_VAL_glaku_9y <- 2 * (1 - pnorm(abs(adhd_glaku_9y.data$t_glaku_9y)))

# Merge all cohorts birth
cohorts.list <- list(adhd_genr_birth.data, adhd_alspac_birth.data, adhd_inma_birth.data, adhd_predo_birth.data, adhd_genr_9y.data, adhd_alspac_9y.data, adhd_helix_pak_9y.data, adhd_helix_all_eu_9y.data, adhd_glaku_9y.data)

adhd.data <- Reduce(function(x, y) merge(x, y, by="probeID", all=TRUE), cohorts.list)

# Function to run meta-regression
meta_regression <- function(cpg) {
  tryCatch({
    probe.data <- adhd.data[adhd.data$probeID == cpg, ]
    
    beta <- probe.data[,.(beta_genr_birth,beta_alspac_birth,beta_inma_birth,beta_predo_birth,beta_genr_9y,beta_alspac_9y,beta_helix_all_eu_9y,beta_helix_pak_9y,beta_glaku_9y)]
    beta <- as.numeric(beta)
    se <- probe.data[,.(se_genr_birth,se_alspac_birth,se_inma_birth,se_predo_birth,se_genr_9y,se_alspac_9y,se_helix_all_eu_9y,se_helix_pak_9y,se_glaku_9y)]
    se <- as.numeric(se)
    cohort <- c("GenR","ALSPAC","HELIX_eu","PREDO","GenR","ALSPAC","HELIX_eu","HELIX_pak","GLAKU")
    age <- c(0,0,0,0,10,7,8,7,12)
    age_school <- c(0,0,0,0,1,1,1,1,1)
    age_birth <- c(1,1,1,1,0,0,0,0,0)
    probe_long.data <- data.frame(beta,se,cohort,age,age_school,age_birth)
    
    cont.fit <- rma.mv(beta, se^2, mods = ~ age,
                       random = ~ 1 | cohort,
                       data=probe_long.data)
    
    school.fit <- rma.mv(beta, se^2, mods = ~ age_school,
                         random = ~ 1 | cohort,
                         data=probe_long.data)
    
    birth.fit <- rma.mv(beta, se^2, mods = ~ age_birth,
                        random = ~ 1 | cohort,
                        data=probe_long.data)
    
    birth <- coef(summary(cont.fit))[1,]
    school_age <- coef(summary(cont.fit))[2,]
    birth_cat <- coef(summary(school.fit))[1,]
    school_age_cat <- coef(summary(school.fit))[2,]
    school_age_cat_rev <- coef(summary(birth.fit))[1,]
    birth_cat_rev <- coef(summary(birth.fit))[2,]
    results.data <- data.frame(cpg, birth[c("estimate","se","pval")], school_age[c("estimate","se","pval")], birth_cat[c("estimate","se","pval")], school_age_cat[c("estimate","se","pval")], school_age_cat_rev[c("estimate","se","pval")], birth_cat_rev[c("estimate","se","pval")])
    names(results.data) <- c("cpg","estimate_birth","se_birth","pval_birth","estimate_year_delta","se_year_delta","pval_year_delta","estimate_birth_cat","se_birth_cat","pval_birth_cat","estimate_school_delta_cat","se_school_delta_cat","pval_school_delta_cat","estimate_school_cat_rev","se_school_cat_rev","pval_school_cat_rev","estimate_birth_delta_rev","se_birth_delta_rev","pval_birth_delta_rev")
    return(results.data)
  }, error = function(err) {
    estimate_birth <- NA
    se_birth <- NA
    pval_birth <- NA
    estimate_year_delta <- NA
    se_year_delta <- NA
    pval_year_delta <- NA
    results.data <- data.frame(cpg,estimate_birth,se_birth,pval_birth,estimate_year_delta,se_year_delta,pval_year_delta)
  })
}

### Create a selection of CpG to test
# Filter by n and study number
n_birth <- rowSums(adhd.data[,.(n_genr_birth,n_alspac_birth,n_inma_birth,n_predo_birth)], na.rm = T) 
n_studies_birth <- rowSums(!is.na(adhd.data[,.(n_genr_birth,n_alspac_birth,n_inma_birth,n_predo_birth)]), na.rm = T) 
cpg_birth <- adhd.data$probeID[n_birth >= 1000 & n_studies_birth > 3]

n_9y <- rowSums(adhd.data[,.(n_genr_9y,n_alspac_9y,n_helix_all_eu_9y,n_helix_pak_9y,n_glaku_9y)], na.rm = T) 
n_studies_9y <- rowSums(!is.na(adhd.data[,.(n_genr_9y,n_alspac_9y,n_helix_all_eu_9y,n_helix_pak_9y,n_glaku_9y)]), na.rm = T) 
cpg_9y <- adhd.data$probeID[n_9y >= 1000 & n_studies_9y > 3]

cpg_qc <- cpg_birth[cpg_birth %in% cpg_9y]

# Autosomal probes
annotation.data <- read_csv("datasets/humanmethylation450_ANNOTATION.csv")
aut <- annotation.data$IlmnID[annotation.data$CHR != "X" & annotation.data$CHR != "Y"]

cpg_qc_aut <- cpg_qc[cpg_qc %in% aut]

# Remove Cross-reactive probes
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "450K")
cpg_qc_aut_no_cross <- cpg_qc_aut[!(cpg_qc_aut %in% xloci)]

meta_results.list <- pbmclapply(cpg_qc_aut_no_cross, meta_regression)

meta_results.data <- bind_rows(meta_results.list)
meta_results.data <- meta_results.data[complete.cases(meta_results.data), ]
save(meta_results.data, file = "results/meta_results_adhd_equal_n.Rdata")

# Absolute effect sizes
meta_results.data$predicted_10y <- meta_results.data$estimate_birth + 10*meta_results.data$estimate_year_delta

meta_results.data$estimate_birth_abs <- abs(meta_results.data$estimate_birth)
meta_results.data$predicted_10y_abs <- abs(meta_results.data$predicted_10y)
meta_results.data$dif_abs <- meta_results.data$predicted_10y_abs - meta_results.data$estimate_birth_abs

meta_results.data$estimate_birth_cat_abs <- abs(meta_results.data$estimate_birth_cat)
meta_results.data$estimate_school_cat_rev_abs <- abs(meta_results.data$estimate_school_cat_rev)
meta_results.data$dif_abs_cat <- meta_results.data$estimate_school_cat_rev_abs - meta_results.data$estimate_birth_cat_abs

# Z scores
meta_results.data$z_birth_cat <- meta_results.data$estimate_birth_cat/meta_results.data$se_birth_cat
meta_results.data$z_school_cat_rev <- meta_results.data$estimate_school_cat_rev/meta_results.data$se_school_cat_rev

meta_results.data$z_birth_cat_abs <- abs(meta_results.data$z_birth_cat)
meta_results.data$z_school_cat_rev_abs <- abs(meta_results.data$z_school_cat_rev)
meta_results.data$dif_abs_cat_z <- meta_results.data$z_school_cat_rev_abs - meta_results.data$z_birth_cat_abs 

##### Statistics report
# n cpg
dim(meta_results.data)
# average effect size birth
mean(meta_results.data[,"estimate_birth_cat_abs"])
# average SE birth
mean(meta_results.data[,"se_birth_cat"])
# average z birth
mean(meta_results.data[,"z_birth_cat_abs"])
# nominal significant cpg birth
dim(meta_results.data[meta_results.data$pval_birth_cat < 0.05,])
# FDR significant
dim(meta_results.data[p.adjust(meta_results.data$pval_birth_cat, method = "fdr") < 0.05,])
# Bonferroni significant
dim(meta_results.data[p.adjust(meta_results.data$pval_birth_cat) < 0.05,])

# average effect size school
mean(meta_results.data[,"estimate_school_cat_rev_abs"])
# average SE school
mean(meta_results.data[,"se_school_cat_rev"])
# average z school
mean(meta_results.data[,"z_school_cat_rev_abs"])
# nominal significant cpg school
dim(meta_results.data[meta_results.data$pval_school_cat_rev < 0.05,])
# FDR significant
dim(meta_results.data[p.adjust(meta_results.data$pval_school_cat_rev, method = "fdr") < 0.05,])
# Bonferroni significant
dim(meta_results.data[p.adjust(meta_results.data$pval_school_cat_rev) < 0.05,])

### Change continuous
mean(meta_results.data$dif_abs)
dim(meta_results.data[meta_results.data$pval_year_delta < 0.05 & meta_results.data$dif_abs > 0,])
dim(meta_results.data[meta_results.data$pval_year_delta < 0.05 & meta_results.data$dif_abs < 0,])

# Change categorical
mean(meta_results.data$dif_abs_cat)
dim(meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & meta_results.data$dif_abs_cat > 0,])
dim(meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & meta_results.data$dif_abs_cat < 0,])
dim(meta_results.data[p.adjust(meta_results.data$pval_school_delta_cat, method = "fdr") < 0.05 & meta_results.data$dif_abs_cat > 0,])
dim(meta_results.data[p.adjust(meta_results.data$pval_school_delta_cat, method = "fdr") < 0.05 & meta_results.data$dif_abs_cat < 0,])
dim(meta_results.data[p.adjust(meta_results.data$pval_school_delta_cat) < 0.05 & meta_results.data$dif_abs_cat > 0,])
dim(meta_results.data[p.adjust(meta_results.data$pval_school_delta_cat) < 0.05 & meta_results.data$dif_abs_cat < 0,])

### Violin Plot
violin_adhd_equal_n.data <- data.frame(c(rep("adhd",795488)),c(rep("birth",397744),rep("school-age",397744)),c(meta_results.data$cpg,meta_results.data$cpg),c(meta_results.data$estimate_birth_abs,meta_results.data$estimate_school_cat_rev_abs),c(meta_results.data$z_birth_cat_abs,meta_results.data$z_school_cat_rev_abs),c(meta_results.data$pval_school_delta_cat,meta_results.data$pval_school_delta_cat))
names(violin_adhd_equal_n.data) <- c("outcome","age","cpg","estimate","z","pval_school_delta_cat")
violin_adhd_equal_n.data$age <- as.factor(violin_adhd_equal_n.data$age)
save(violin_adhd_equal_n.data, file = "results/violin_adhd_equal_n.Rdata")

