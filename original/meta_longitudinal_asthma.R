library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)
library(ggplot2)

# Load summary data from all cohorts
asthma_alspac_birth.data <- fread("sumstats/asthma/ALSPAC_newborns_model2b.csv")
asthma_medall_school.data <- fread("sumstats/asthma/MeDALL.BAMSE_olderkids_model2b.csv")
asthma_epigene_school.data <- fread("sumstats/asthma/EpiGene.BAMSE_olderkids_model2b.csv")
asthma_chs_birth.data <- fread("sumstats/asthma/CHS_CBCT_newborns_model2b.csv")
asthma_chop_school.data <- fread("sumstats/asthma/CHOP_olderkids_model2b.csv")
asthma_eden_birth.data <- fread("sumstats/asthma/EDEN_CBCT_newborns_model2b.csv")
asthma_gala2_school.data <- fread("sumstats/asthma/GALA2_olderkids_model2b.csv")
asthma_genr_birth.data <- fread("sumstats/asthma/GenR_CBCT_newborns_model2b.csv")
asthma_icac_school.data <- fread("sumstats/asthma/ICAC_olderkids_model2b.csv")
asthma_moba1_birth.data <- fread("sumstats/asthma/MoBa1_CBCT_newborns_model2b.csv")
asthma_moba2_birth.data <- fread("sumstats/asthma/MoBa2_CBCT_newborns_model2b.csv")
asthma_nest_birth.data <- fread("sumstats/asthma/NEST_CBCT_newborns_model2b.csv")
asthma_nfbc_school.data <- fread("sumstats/asthma/NFBC_olderkids_model2b.csv")
asthma_piama_school.data <- fread("sumstats/asthma/PIAMA_olderkids_model2b.csv")
asthma_raine_school.data <- fread("sumstats/asthma/Raine_olderkids_model2b.csv")
asthma_stoppa_school.data <- fread("sumstats/asthma/STOPPA_olderkids_model2b.csv")

# CHS is completely duplicated for some reason, keep first observation
asthma_chs_birth.data <- asthma_chs_birth.data[!duplicated(asthma_chs_birth.data)]

# Raine has empty rows, remove them
asthma_raine_school.data <- asthma_raine_school.data[!is.na(asthma_raine_school.data$CpG), ]

# Add missing n
asthma_medall_school.data$N <- 214
asthma_epigene_school.data$N <- 307
asthma_chs_birth.data$N <- 229
asthma_chop_school.data$N <- 382
asthma_eden_birth.data$N <- 150
asthma_gala2_school.data$N <- 193
asthma_genr_birth.data$N <- 661
asthma_moba1_birth.data$N <- 666
asthma_moba2_birth.data$N <- 458
asthma_piama_school.data$N <- 197
asthma_raine_school.data$N <- 509
asthma_stoppa_school.data$N <- 460

# Number of cases
# asthma_alspac_birth.data$N_cases <- 88
# asthma_chs_birth.data$N_cases <- 39
# asthma_eden_birth.data$N_cases <- 34
# asthma_genr_birth.data$N_cases <- 37
# asthma_moba1_birth.data$N_cases <- 149
# asthma_moba2_birth.data$N_cases <- 239
# asthma_nest_birth.data$N_cases <- 45
# sum(88+39+34+37+149+239+45) #by coincidence same amount of cases as in children
# 
# asthma_epigene_school.data$N_cases <- 93
# asthma_medall_school.data$N_cases <- 47
# asthma_chop_school.data$N_cases <- 19
# asthma_gala2_school.data$N_cases <- 106
# asthma_icac_school.data$N_cases <- 92
# asthma_nfbc_school.data$N_cases <- 17
# asthma_piama_school.data$N_cases <- 15
# asthma_raine_school.data$N_cases <- 105
# asthma_stoppa_school.data$N_cases <- 137


# Rename columns
setnames(asthma_alspac_birth.data, c("coef", "se", "pvalue", "N_for_probe"), c("beta_alspac_birth", "se_alspac_birth", "p_alspac_birth", "n_alspac_birth"))
setnames(asthma_medall_school.data, c("coef", "se", "pvalue", "N"), c("beta_medall_school", "se_medall_school", "p_medall_school", "n_medall_school"))
setnames(asthma_epigene_school.data, c("coef", "se", "pvalue", "N"), c("beta_epigene_school", "se_epigene_school", "p_epigene_school", "n_epigene_school"))
setnames(asthma_chs_birth.data, c("coef", "se", "pvalue", "N"), c("beta_chs_birth", "se_chs_birth", "p_chs_birth", "n_chs_birth"))
setnames(asthma_chop_school.data, c("coef", "se", "pvalue", "N"), c("beta_chop_school", "se_chop_school", "p_chop_school", "n_chop_school"))
setnames(asthma_eden_birth.data, c("coef", "se", "pvalue", "N"), c("beta_eden_birth", "se_eden_birth", "p_eden_birth", "n_eden_birth"))
setnames(asthma_gala2_school.data, c("coef", "se", "pvalue", "N"), c("beta_gala2_school", "se_gala2_birth", "p_gala2_school", "n_gala2_school"))
setnames(asthma_genr_birth.data, c("coef", "se", "pvalue", "N"), c("beta_genr_birth", "se_genr_birth", "p_genr_birth", "n_genr_birth"))
setnames(asthma_icac_school.data, c("coef", "se", "pvalue", "N"), c("beta_icac_school", "se_icac_school", "p_icac_school", "n_icac_school"))
setnames(asthma_moba1_birth.data, c("coef", "se", "pvalue", "N"), c("beta_moba1_birth", "se_moba1_birth", "p_moba1_birth", "n_moba1_birth"))
setnames(asthma_moba2_birth.data, c("coef", "se", "pvalue", "N"), c("beta_moba2_birth", "se_moba2_birth", "p_moba2_birth", "n_moba2_birth"))
setnames(asthma_nest_birth.data, c("coef", "se", "pvalue", "N_for_probe"), c("beta_nest_birth", "se_nest_birth", "p_nest_birth", "n_nest_birth"))
setnames(asthma_nfbc_school.data, c("coef", "se", "pvalue", "N"), c("beta_nfbc_school", "se_nfbc_school", "p_nfbc_school", "n_nfbc_school"))
setnames(asthma_piama_school.data, c("coef", "se", "pvalue", "N"), c("beta_piama_school", "se_piama_school", "p_piama_school", "n_piama_school"))
setnames(asthma_raine_school.data, c("coef", "se", "pvalue", "N"), c("beta_raine_school", "se_raine_school", "p_raine_school", "n_raine_school"))
setnames(asthma_stoppa_school.data, c("coef", "se", "pvalue", "N"), c("beta_stoppa_school", "se_stoppa_school", "p_stoppa_school", "n_stoppa_school"))

# Merge all cohorts
cohorts.list <- list(asthma_alspac_birth.data, asthma_medall_school.data, asthma_epigene_school.data, asthma_chs_birth.data,
                     asthma_chop_school.data, asthma_eden_birth.data, asthma_gala2_school.data, asthma_genr_birth.data, asthma_icac_school.data,
                     asthma_moba1_birth.data, asthma_moba2_birth.data, asthma_nest_birth.data, asthma_nfbc_school.data,
                     asthma_piama_school.data, asthma_raine_school.data, asthma_stoppa_school.data)

asthma.data <- Reduce(function(x, y) merge(x, y, by="CpG", all=TRUE), cohorts.list)

# Function to run meta-regression
meta_regression <- function(cpg) {
  tryCatch({
    # Select cpg to analyze
    probe.data <- asthma.data[asthma.data$CpG == cpg, ]
    
    # Define betas
    beta <- probe.data[,.(beta_alspac_birth, beta_medall_school, beta_epigene_school, beta_chs_birth,
                          beta_chop_school, beta_eden_birth, beta_gala2_school, beta_genr_birth,
                          beta_moba1_birth, beta_moba2_birth, beta_nest_birth, beta_nfbc_school,
                          beta_piama_school, beta_raine_school, beta_stoppa_school, beta_icac_school)]
    beta <- as.numeric(beta)
    
    # Define SE
    se <- probe.data[,.(se_alspac_birth, se_medall_school, se_epigene_school, se_chs_birth,
                        se_chop_school, se_eden_birth, se_gala2_birth, se_genr_birth,
                        se_moba1_birth, se_moba2_birth, se_nest_birth, se_nfbc_school,
                        se_piama_school, se_raine_school, se_stoppa_school, se_icac_school)]
    se <- as.numeric(se)
    
    
    cohort <- c("ALSPAC","MEDALL","EPIGENE","CHS","CHOP","EDEN","GALA2","GenR","MOBA1","MOBA2","NEST","NFBC","PIAMA","RAINE","STOPPA","ICAC")
    
    age <- c(0,8.37,8.33,0,7.10,0,9.31,0,0,0,0,16.00,8.06,17.01,12.54,9)
    age_school <- c(0,1,1,0,1,0,1,0,0,0,0,1,1,1,1,1)
    age_birth <- c(1,0,0,1,0,1,0,1,1,1,1,0,0,0,0,0)
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
n_birth <- rowSums(asthma.data[,.(n_alspac_birth, n_chs_birth,
                                  n_eden_birth, n_genr_birth,
                                  n_moba1_birth, n_moba2_birth, n_nest_birth)], na.rm = T) 
n_studies_birth <- rowSums(!is.na(asthma.data[,.(n_alspac_birth, n_chs_birth,
                                                 n_eden_birth, n_genr_birth,
                                                 n_moba1_birth, n_moba2_birth, n_nest_birth)]), na.rm = T) 
cpg_birth <- asthma.data$CpG[n_birth >= 1000 & n_studies_birth > 3]

n_school <- rowSums(asthma.data[,.(n_medall_school, n_epigene_school,
                                   n_chop_school, n_gala2_school, n_nfbc_school,
                                   n_piama_school, n_raine_school, n_stoppa_school, n_icac_school)], na.rm = T) 
n_studies_school <- rowSums(!is.na(asthma.data[,.(n_medall_school, n_epigene_school,
                                                  n_chop_school, n_gala2_school, n_nfbc_school,
                                                  n_piama_school, n_raine_school, n_stoppa_school, n_icac_school)]), na.rm = T) 
cpg_school <- asthma.data$CpG[n_school >= 1000 & n_studies_school > 3]

cpg_qc <- cpg_birth[cpg_birth %in% cpg_school]

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
save(meta_results.data, file = "results/meta_results_asthma.Rdata")

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
## Effect size change (n=27600)
changing.data <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05,]
# Positive constant (n=5)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Positive to null (n=1961)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] & 
                    abs(changing.data$estimate_school_cat_rev) < effect_dist[9],])
# Null to positive (n=9249)
dim(changing.data[abs(changing.data$estimate_birth_cat) < effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Negative constant (n=1)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
# Negative to null (n=1108)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    abs(changing.data$estimate_school_cat_rev) < effect_dist[9],])
# Null to negative (n=3267)
dim(changing.data[abs(changing.data$estimate_birth_cat) < effect_dist[9] & 
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
 # Positive to negative (n=5535)
dim(changing.data[changing.data$estimate_birth_cat >= effect_dist[9] &
                    changing.data$estimate_school_cat_rev <= -effect_dist[9],])
# Negative to positive (n=3441)
dim(changing.data[changing.data$estimate_birth_cat <= -effect_dist[9] & 
                    changing.data$estimate_school_cat_rev >= effect_dist[9],])
# Null to Null (n=9070)
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
violin_asthma.data <- data.frame(c(rep("asthma",865456)),c(rep("birth",432728),rep("school-age",432728)),c(meta_results.data$cpg,meta_results.data$cpg),c(meta_results.data$estimate_birth_cat_abs,meta_results.data$estimate_school_cat_rev_abs),c(meta_results.data$z_birth_cat_abs,meta_results.data$z_school_cat_rev_abs),c(meta_results.data$pval_school_delta_cat,meta_results.data$pval_school_delta_cat))
names(violin_asthma.data) <- c("outcome","age","cpg","estimate","z","pval_school_delta_cat")
violin_asthma.data$age <- as.factor(violin_asthma.data$age)
save(violin_asthma.data, file = "results/violin_asthma.Rdata")

###### Heterogeneity analyses
# Function to run meta-regression and estimate heterogeneity
meta_regression_het <- function(cpg) {
  tryCatch({
    # Select cpg to analyze
    probe.data <- asthma.data[asthma.data$CpG == cpg, ]
    
    # Define betas
    beta <- probe.data[,.(beta_alspac_birth, beta_medall_school, beta_epigene_school, beta_chs_birth,
                          beta_chop_school, beta_eden_birth, beta_gala2_school, beta_genr_birth,
                          beta_moba1_birth, beta_moba2_birth, beta_nest_birth, beta_nfbc_school,
                          beta_piama_school, beta_raine_school, beta_stoppa_school, beta_icac_school)]
    beta <- as.numeric(beta)
    
    # Define SE
    se <- probe.data[,.(se_alspac_birth, se_medall_school, se_epigene_school, se_chs_birth,
                        se_chop_school, se_eden_birth, se_gala2_birth, se_genr_birth,
                        se_moba1_birth, se_moba2_birth, se_nest_birth, se_nfbc_school,
                        se_piama_school, se_raine_school, se_stoppa_school, se_icac_school)]
    se <- as.numeric(se)
    
    
    cohort <- c("ALSPAC","MEDALL","EPIGENE","CHS","CHOP","EDEN","GALA2","GenR","MOBA1","MOBA2","NEST","NFBC","PIAMA","RAINE","STOPPA","ICAC")
    
    age_school <- factor(c(0,1,1,0,1,0,1,0,0,0,0,1,1,1,1,1))
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
save(meta_results_het.data, file = "results/meta_results_het_asthma.Rdata")

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
