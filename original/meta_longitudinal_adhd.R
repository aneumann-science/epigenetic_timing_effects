library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)

# Load summary data from all cohorts
adhd_genr_birth.data <- fread("sumstats/ADHD/adhd_ewas_genr_20180221.txt")
adhd_alspac_birth.data <- fread("sumstats/ADHD/adhd_alspac_7ct_plate_noPCs_trimmed_20161010.txt")
adhd_inma_birth.data <- fread("sumstats/ADHD/adhd_ewas_inma_20161021.txt")
adhd_nest_white_birth.data <- fread("sumstats/ADHD/ADHD_cont_House_NEST_white_11142016.txt")
adhd_nest_black_birth.data <- fread("sumstats/ADHD/ADHD_cont_House_NEST_black_11142016.txt")
adhd_predo_birth.data <- fread("sumstats/ADHD/ChildADHD_PREDO_10012017_trimmed_df.txt")

# Rename columns
setnames(adhd_genr_birth.data, c("beta", "se", "t", "n"), c("beta_genr_birth", "se_genr_birth", "t_genr_birth", "n_genr_birth"))
setnames(adhd_alspac_birth.data, c("beta", "se", "t", "n"), c("beta_alspac_birth", "se_alspac_birth", "t_alspac_birth", "n_alspac_birth"))
setnames(adhd_inma_birth.data, c("beta", "se", "t", "n"), c("beta_inma_birth", "se_inma_birth", "t_inma_birth", "n_inma_birth"))
setnames(adhd_nest_white_birth.data, c("BETA", "SE", "P_VAL", "N_for_probe", "N"), c("beta_nest_white_birth", "se_nest_white_birth", "P_VAL_nest_white_birth", "N_for_probe_nest_white_birth", "n_nest_white_birth"))
setnames(adhd_nest_black_birth.data, c("BETA", "SE", "P_VAL", "N_for_probe", "N"), c("beta_nest_black_birth", "se_nest_black_birth", "P_VAL_nest_black_birth", "N_for_probe_nest_black_birth", "n_nest_black_birth"))
setnames(adhd_predo_birth.data, c("BETA", "SE", "T_VAL", "N"), c("beta_predo_birth", "se_predo_birth", "t_predo_birth", "n_predo_birth"))

# Calculate approximate p-values
adhd_genr_birth.data$P_VAL_genr_birth <- 2 * (1 - pnorm(abs(adhd_genr_birth.data$t_genr_birth)))
adhd_alspac_birth.data$P_VAL_alspac_birth <- 2 * (1 - pnorm(abs(adhd_alspac_birth.data$t_alspac_birth)))
adhd_inma_birth.data$P_VAL_inma_birth <- 2 * (1 - pnorm(abs(adhd_inma_birth.data$t_inma_birth)))
adhd_predo_birth.data$P_VAL_predo_birth <- 2 * (1 - pnorm(abs(adhd_predo_birth.data$t_predo_birth)))

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

# Merge all cohorts
cohorts.list <- list(adhd_genr_birth.data, adhd_alspac_birth.data, adhd_inma_birth.data, adhd_nest_white_birth.data, adhd_nest_black_birth.data, adhd_predo_birth.data, adhd_genr_9y.data, adhd_alspac_9y.data, adhd_helix_pak_9y.data, adhd_helix_all_eu_9y.data, adhd_glaku_9y.data)

adhd.data <- Reduce(function(x, y) merge(x, y, by="probeID", all=TRUE), cohorts.list)

# Function to run meta-regression
meta_regression <- function(cpg) {
  tryCatch({
    # Extract data for a specific DNAm site
    probe.data <- adhd.data[adhd.data$probeID == cpg, ]
    
    # Re-arrange data into a long format with column: Beta, SE, cohort, age (continuous), time point (0=birth, 1=childhood), time point (0=birth, 1=childhood), time point (0=childhood, 1=birth) 
    beta <- probe.data[,.(beta_genr_birth,beta_alspac_birth,beta_inma_birth,beta_nest_white_birth,beta_nest_black_birth,beta_predo_birth,beta_genr_9y,beta_alspac_9y,beta_helix_all_eu_9y,beta_helix_pak_9y,beta_glaku_9y)]
    beta <- as.numeric(beta)
    se <- probe.data[,.(se_genr_birth,se_alspac_birth,se_inma_birth,se_nest_white_birth,se_nest_black_birth,se_predo_birth,se_genr_9y,se_alspac_9y,se_helix_all_eu_9y,se_helix_pak_9y,se_glaku_9y)]
    se <- as.numeric(se)
    cohort <- c("GenR","ALSPAC","HELIX_eu","NEST_white","NEXT_black","PREDO","GenR","ALSPAC","HELIX_eu","HELIX_pak","GLAKU")
    age <- c(0,0,0,0,0,0,10,7,8,7,12)
    age_school <- c(0,0,0,0,0,0,1,1,1,1,1)
    age_birth <- c(1,1,1,1,1,1,0,0,0,0,0)
    probe_long.data <- data.frame(beta,se,cohort,age,age_school,age_birth)
    
    # Meta-regression on continuous age
    cont.fit <- rma.mv(beta, se^2, mods = ~ age,
                       random = ~ 1 | cohort,
                       data=probe_long.data)
    
    # Meta-regression on time point (childhood effect)
    school.fit <- rma.mv(beta, se^2, mods = ~ age_school,
                         random = ~ 1 | cohort,
                         data=probe_long.data)
    
    # Meta-regression on time point (birth effect)
    birth.fit <- rma.mv(beta, se^2, mods = ~ age_birth,
                        random = ~ 1 | cohort,
                        data=probe_long.data)
    
    ### extract all results
    # Continuous age model
    birth <- coef(summary(cont.fit))[1,]
    school_age <- coef(summary(cont.fit))[2,]
    # Childhood effect (birth reference)
    birth_cat <- coef(summary(school.fit))[1,]
    school_age_cat <- coef(summary(school.fit))[2,]
    # Birth effect (school reference)
    school_age_cat_rev <- coef(summary(birth.fit))[1,]
    birth_cat_rev <- coef(summary(birth.fit))[2,]
    # Format results
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
n_birth <- rowSums(adhd.data[,.(n_genr_birth,n_alspac_birth,n_inma_birth,n_nest_white_birth,n_nest_black_birth,n_predo_birth)], na.rm = T) 
n_studies_birth <- rowSums(!is.na(adhd.data[,.(n_genr_birth,n_alspac_birth,n_inma_birth,n_nest_white_birth,n_nest_black_birth,n_predo_birth)]), na.rm = T) 
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

###### Perform meta-regression for all probes
meta_results.list <- pbmclapply(cpg_qc_aut_no_cross, meta_regression)

meta_results.data <- bind_rows(meta_results.list)
meta_results.data <- meta_results.data[complete.cases(meta_results.data), ]
save(meta_results.data, file = "results/meta_results_adhd.Rdata")

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
## average standard error
#sqrt(sum(meta_results.data$se_birth_cat^2))/length(meta_results.data$se_birth_cat)
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

# Distribution absolute effect size at either time point
effect_dist <- quantile(c(meta_results.data$estimate_birth_cat_abs,meta_results.data$estimate_school_cat_rev_abs), probs = c(0:10/10,0.95,0.975,0.99))
round(effect_dist, 2)
effect_dist_birth <- quantile(meta_results.data$estimate_birth_cat_abs, probs = c(0:10/10,0.95,0.975,0.99))
round(effect_dist_birth, 2)
effect_childhood <- quantile(meta_results.data$estimate_school_cat_rev_abs, probs = c(0:10/10,0.95,0.975,0.99))
round(effect_childhood, 2)

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

### Change patterns
## Effect size change (n=17383)
changing.data <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05,]
# Positive constant (n=189)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Positive to null (n=748)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] & 
                    abs(changing.data$estimate_school_cat_rev) < effect_dist[9],])
# Null to positive (n=4254)
dim(changing.data[abs(changing.data$estimate_birth_cat) < effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Negative constant (n=45)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
# Negative to null (n=1748)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    abs(changing.data$estimate_school_cat_rev) < effect_dist[9],])
# Null to negative (n=2269)
dim(changing.data[abs(changing.data$estimate_birth_cat) < effect_dist[9] & 
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
# Positive to negative (n=5535)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] &
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
# Negative to positive (n=1265)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Null to Null (n=5069)
dim(changing.data[abs(changing.data$estimate_birth_cat) < effect_dist[9] & 
                    abs(changing.data$estimate_school_cat_rev) < effect_dist[9],])

# Ratios of increasing vs decreasing DNAm sites at different p-value thresholds
alphas <- c(1,0.5,0.05,0.01,1E-03,1E-04,1E-05,1E-06,1E-07,5E-08,1E-08)
increase_decrease_ratio.data <- sapply(alphas, function(alpha) {
  n_increase <- dim(meta_results.data[meta_results.data$pval_school_delta_cat < alpha & meta_results.data$dif_abs_cat > 0,])[1]
  n_decrease <- dim(meta_results.data[meta_results.data$pval_school_delta_cat < alpha & meta_results.data$dif_abs_cat < 0,])[1]
  increase_descrease_ratio <- n_increase/n_decrease
})


###### Uncorrelated CpGs 
###### (not reported in manuscript due to space/complexity
###### but generally similar results to original analysis)
# Load list of correlated CpGs
load("datasets/methylation_birth_mcb_r03.Rdata")

# Correlated CpGs
cor_cpgs <-methylation_birth_mcb_r03$MCBsites 

# cor_cpgs <-methylation_birth_mcb_r03$MCBsites 
mcb.list <- methylation_birth_mcb_r03$MCBinformation[,4]

# Pick per block all CpGs except randomly 1
mcb_random.list <- lapply(mcb.list, function(mcb) {
  cpgs <- unlist(strsplit(mcb, " "))
  sample(cpgs, length(cpgs)-1)
})

mcb_random <- unlist(mcb_random.list)

meta_results_ind.data <- meta_results.data[meta_results.data$cpg %in% mcb_random, ]

##### Statistics report
# n cpg
dim(meta_results_ind.data)
# average effect size birth
mean(meta_results_ind.data[,"estimate_birth_cat_abs"])
# average z birth
mean(meta_results_ind.data[,"z_birth_cat_abs"])
# nominal significant cpg birth
dim(meta_results_ind.data[meta_results_ind.data$pval_birth_cat < 0.05,])
# FDR significant
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_birth_cat, method = "fdr") < 0.05,])
# Bonferroni significant
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_birth_cat) < 0.05,])

# average effect size school
mean(meta_results_ind.data[,"estimate_school_cat_rev_abs"])
# average z school
mean(meta_results_ind.data[,"z_school_cat_rev_abs"])
# nominal significant cpg school
dim(meta_results_ind.data[meta_results_ind.data$pval_school_cat_rev < 0.05,])
# FDR significant
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_cat_rev, method = "fdr") < 0.05,])
# Bonferroni significant
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_cat_rev) < 0.05,])

# Change categorical
mean(meta_results_ind.data$dif_abs_cat)
dim(meta_results_ind.data[meta_results_ind.data$pval_school_delta_cat < 0.05 & meta_results_ind.data$dif_abs_cat > 0,])
dim(meta_results_ind.data[meta_results_ind.data$pval_school_delta_cat < 0.05 & meta_results_ind.data$dif_abs_cat < 0,])
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_delta_cat, method = "fdr") < 0.05 & meta_results_ind.data$dif_abs_cat > 0,])
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_delta_cat, method = "fdr") < 0.05 & meta_results_ind.data$dif_abs_cat < 0,])
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_delta_cat) < 0.05 & meta_results_ind.data$dif_abs_cat > 0,])
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_school_delta_cat) < 0.05 & meta_results_ind.data$dif_abs_cat < 0,])
 
### Violin Plot
# All
violin_adhd.data <- data.frame(c(rep("adhd",860654)),c(rep("birth",430327),rep("school-age",430327)),c(meta_results.data$cpg,meta_results.data$cpg),c(meta_results.data$estimate_birth_cat_abs,meta_results.data$estimate_school_cat_rev_abs),c(meta_results.data$z_birth_cat_abs,meta_results.data$z_school_cat_rev_abs),c(meta_results.data$pval_school_delta_cat,meta_results.data$pval_school_delta_cat))
names(violin_adhd.data) <- c("outcome","age","cpg","estimate","z","pval_school_delta_cat")
violin_adhd.data$age <- as.factor(violin_adhd.data$age)
save(violin_adhd.data, file = "results/violin_adhd.Rdata")

###### Heterogeneity analyses
# Function to run meta-regression and estimate heterogeneity
meta_regression_het <- function(cpg) {
  tryCatch({
    probe.data <- adhd.data[adhd.data$probeID == cpg, ]
    
    beta <- probe.data[,.(beta_genr_birth,beta_alspac_birth,beta_inma_birth,beta_nest_white_birth,beta_nest_black_birth,beta_predo_birth,beta_genr_9y,beta_alspac_9y,beta_helix_all_eu_9y,beta_helix_pak_9y,beta_glaku_9y)]
    beta <- as.numeric(beta)
    se <- probe.data[,.(se_genr_birth,se_alspac_birth,se_inma_birth,se_nest_white_birth,se_nest_black_birth,se_predo_birth,se_genr_9y,se_alspac_9y,se_helix_all_eu_9y,se_helix_pak_9y,se_glaku_9y)]
    se <- as.numeric(se)
    cohort <- c("GenR","ALSPAC","HELIX_eu","NEST_white","NEXT_black","PREDO","GenR","ALSPAC","HELIX_eu","HELIX_pak","GLAKU")
    age_school <- factor(c(0,0,0,0,0,0,1,1,1,1,1))
    probe_long.data <- data.frame(beta,se,cohort,age_school)
    
    school.fit <- rma.mv(beta, se^2, mods = ~ 0 + age_school,
                         random = ~ age_school | cohort,
                         struct="DIAG", cvvc=T,
                         data=probe_long.data)
    
    results.data <- data.frame(cpg, t(coef(school.fit)), t(sqrt(diag(vcov(school.fit)))),abs(t(school.fit$zval)),t(school.fit$pval),t(sqrt(school.fit$tau2)))
    names(results.data) <- c("cpg","beta_birth","beta_school","se_birth","se_school","z_birth_abs","z_school_abs","p_birth","p_school","tau_birth","tau_school")
    results.data$beta_birth_abs <- abs(results.data$beta_birth)
    results.data$beta_school_abs <- abs(results.data$beta_school)
    results.data <- results.data[c("cpg","beta_birth","beta_birth_abs","tau_birth","se_birth","z_birth_abs","p_birth","beta_school","beta_school_abs","tau_school","se_school","z_school_abs","p_school")]
    return(results.data)
  }, error = function(err) {
    beta_birth <- NA
    beta_birth_abs <- NA
    tau_birth <- NA
    se_birth <- NA
    z_birth_abs <- NA
    p_birth <- NA
    beta_school <- NA
    beta_school_abs <- NA
    tau_school <- NA
    se_school <- NA
    z_school_abs <- NA
    p_school <- NA
    results.data <- data.frame(cpg,beta_birth,beta_birth_abs,tau_birth,se_birth,z_birth_abs,p_birth,beta_school,beta_school_abs,tau_school,se_school,z_school_abs,p_school)
  })
}

# Obtain heterogeneity statistics for each probe
meta_results_het.list <- pbmclapply(cpg_qc_aut_no_cross, meta_regression_het, mc.cores = 6)

meta_results_het.data <- bind_rows(meta_results_het.list)
meta_results_het.data <- meta_results_het.data[complete.cases(meta_results_het.data), ]
save(meta_results_het.data, file = "results/meta_results_het_adhd.Rdata")

# Statistics all
length(meta_results_het.data$beta_birth)
mean(meta_results_het.data$beta_birth_abs)
mean(meta_results_het.data$tau_birth)
mean(meta_results_het.data$beta_school_abs)
mean(meta_results_het.data$tau_school)

# Statistics nominal
length(meta_results_het.data$beta_birth[meta_results_het.data$p_birth < 0.05])
mean(meta_results_het.data[meta_results_het.data$p_birth < 0.05,"beta_birth_abs"])
mean(meta_results_het.data[meta_results_het.data$p_birth < 0.05,"tau_birth"])
length(meta_results_het.data$beta_birth[meta_results_het.data$p_school < 0.05])
mean(meta_results_het.data[meta_results_het.data$p_school < 0.05,"beta_school_abs"])
mean(meta_results_het.data[meta_results_het.data$p_school < 0.05,"tau_school"])

### Export example data
# DNAm sites with nominally significant change and association at birth or school
change_cpg_adhd <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                       (meta_results.data$pval_birth_cat < 0.05 | 
                                          meta_results.data$pval_school_cat_rev < 0.05), "cpg"]

# Save example data for these data
adhd_example.data <- adhd.data[adhd.data$probeID %in% change_cpg_adhd, ]
save(adhd_example.data, file = "example_data/adhd_example.data")
