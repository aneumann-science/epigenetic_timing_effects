### Example script demonstrating longitudinal meta-regression
# Load all required libraries
library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)

# Load example data for DNAm associations at birth and childhood with adhd symptoms
load("example_data/adhd_example.data")
adhd.data <- adhd_example.data

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
meta_results.list <- pbmclapply(cpg_qc_aut_no_cross, meta_regression, mc.cores = 4)

meta_results.data <- bind_rows(meta_results.list)
meta_results.data <- meta_results.data[complete.cases(meta_results.data), ]

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
meta_results_het.list <- pbmclapply(cpg_qc_aut_no_cross, meta_regression_het, mc.cores = 4)

meta_results_het.data <- bind_rows(meta_results_het.list)
meta_results_het.data <- meta_results_het.data[complete.cases(meta_results_het.data), ]

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
