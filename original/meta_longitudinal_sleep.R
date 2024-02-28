library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)
library(ggplot2)

# Load summary data from all cohorts
sleep_alspac_birth.data <- fread("sumstats/sleep/ALSPAC_20190726/ALSPAC.duration.modela.ewasresults.birth.txt")
sleep_alspac_school.data <- fread("sumstats/sleep/ALSPAC_20190726/ALSPAC.duration.modela.ewasresults.F7.txt")
sleep_eden_birth.data <- fread("sumstats/sleep/EDEN_20190427_cordupdate20201217/EDENewas_res.duration.modelabirth.csv", drop = "V1")
sleep_genr_birth.data <- fread("sumstats/sleep/GenR_20191220_wasoupdate20201206_initupdate20210129/GenR_duration_model1a_birth.csv", drop = "V1")
sleep_genr_school.data <- fread("sumstats/sleep/GenR_20191220_wasoupdate20201206_initupdate20210129/GenR_duration_model4a_postnatal.csv", drop = "V1")
sleep_glaku_school.data <- fread("sumstats/sleep/GLAKU_20190627/GLAKU.ewas_res.duration.modela.postnatal.csv")
sleep_healthystart_birth.data <- fread("sumstats/sleep/HealthyStart_20190701/HealthyStart.ewas_res.duration.modela.birth.csv", drop = "V1")
sleep_helix_school.data <- fread("sumstats/sleep/HELIX_20191021/helix_duration_model4a_postnatal_20191022.csv", drop = "V1")
sleep_inma_birth.data <- fread("sumstats/sleep/INMA_20190401_initupdate20210215/inma_duration_model1a_birth_20190401.csv", drop = "V1")
sleep_lina_birth.data <- fread("sumstats/sleep/LINA_20191214/LINAewas_res.duration.modelabirth_10LJ.csv", drop = "V1")
sleep_moba1_birth.data <- fread("sumstats/sleep/MOBA1_20210113/MoBa1_duration_Model1a_birth_20210113.csv")
sleep_moba2_birth.data <- fread("sumstats/sleep/MOBA2_20210115/MoBa2_duration_Model1a_birth_20210115.csv")
sleep_predo_birth.data <- fread("sumstats/sleep/PREDO_20190531/PREDO_duration_model1a_birth_20190531.csv", drop = "V1")
sleep_viva_birth.data <- fread("sumstats/sleep/VIVA_20190924/Viva_duration_Model1a_birth_2019-09-24.csv", drop = "V1")
sleep_viva_school.data <- fread("sumstats/sleep/VIVA_20190924/Viva_duration_Model4a_7y_2019-09-24.csv", drop = "V1")

# Rename columns
setnames(sleep_alspac_birth.data, c("BETA", "SE", "P_VAL", "N_for_probe"), c("beta_alspac_birth", "se_alspac_birth", "p_alspac_birth", "n_alspac_birth"))
setnames(sleep_alspac_school.data, c("BETA", "SE", "P_VAL", "N_for_probe"), c("beta_alspac_school", "se_alspac_school", "p_alspac_school", "n_alspac_school"))
setnames(sleep_eden_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_eden_birth", "se_eden_birth", "t_eden_birth", "p_eden_birth", "n_eden_birth"))
setnames(sleep_genr_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_genr_birth", "se_genr_birth", "t_genr_birth", "p_genr_birth", "n_genr_birth"))
setnames(sleep_genr_school.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_genr_school", "se_genr_school", "t_genr_school", "p_genr_school", "n_genr_school"))
setnames(sleep_glaku_school.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_glaku_school", "se_glaku_school", "t_glaku_school", "p_glaku_school", "n_glaku_school"))
setnames(sleep_healthystart_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_healthystart_birth", "se_healthystart_birth", "t_healthystart_birth", "p_healthystart_birth", "n_healthystart_birth"))
setnames(sleep_helix_school.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_helix_school", "se_helix_school", "t_helix_school", "p_helix_school", "n_helix_school"))
setnames(sleep_inma_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_inma_birth", "se_inma_birth", "t_inma_birth", "p_inma_birth", "n_inma_birth"))
setnames(sleep_lina_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_lina_birth", "se_lina_birth", "t_lina_birth", "p_lina_birth", "n_lina_birth"))
setnames(sleep_moba1_birth.data, c("beta", "se", "t", "P_value", "N_for_probe"), c("beta_moba1_birth", "se_moba1_birth", "t_moba1_birth", "p_moba1_birth", "n_moba1_birth"))
setnames(sleep_moba2_birth.data, c("beta", "se", "t", "P_value", "N_for_probe"), c("beta_moba2_birth", "se_moba2_birth", "t_moba2_birth", "p_moba2_birth", "n_moba2_birth"))
setnames(sleep_predo_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_predo_birth", "se_predo_birth", "t_predo_birth", "p_predo_birth", "n_predo_birth"))
setnames(sleep_viva_birth.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_viva_birth", "se_viva_birth", "t_viva_birth", "p_viva_birth", "n_viva_birth"))
setnames(sleep_viva_school.data, c("BETA", "SE", "t", "P_VAL", "N_for_probe"), c("beta_viva_school", "se_viva_school", "t_viva_school", "p_viva_school", "n_viva_school"))

# Merge all cohorts
cohorts.list <- list(sleep_alspac_birth.data, sleep_alspac_school.data, sleep_eden_birth.data, sleep_genr_birth.data,
                     sleep_genr_school.data, sleep_glaku_school.data, sleep_healthystart_birth.data, sleep_helix_school.data,
                     sleep_inma_birth.data, sleep_lina_birth.data, sleep_moba1_birth.data, sleep_moba2_birth.data,
                     sleep_predo_birth.data, sleep_viva_birth.data, sleep_viva_school.data)

sleep.data <- Reduce(function(x, y) merge(x, y, by="probeID", all=TRUE), cohorts.list)

# Function to run meta-regression
meta_regression <- function(cpg) {
  tryCatch({
    # Select cpg to analyze
    probe.data <- sleep.data[sleep.data$probeID == cpg, ]
    
    # Define betas
    beta <- probe.data[,.(beta_alspac_birth, beta_alspac_school, beta_eden_birth, beta_genr_birth,
                          beta_genr_school, beta_glaku_school, beta_healthystart_birth, beta_helix_school,
                          beta_inma_birth, beta_lina_birth, beta_moba1_birth, beta_moba2_birth,
                          beta_predo_birth, beta_viva_birth, beta_viva_school)]
    beta <- as.numeric(beta)
    
    # Define SE
    se <- probe.data[,.(se_alspac_birth, se_alspac_school, se_eden_birth, se_genr_birth,
                        se_genr_school, se_glaku_school, se_healthystart_birth, se_helix_school,
                        se_inma_birth, se_lina_birth, se_moba1_birth, se_moba2_birth,
                        se_predo_birth, se_viva_birth, se_viva_school)]
    se <- as.numeric(se)
    
    
    cohort <- c("ALSPAC","ALSPAC","HELIX","GenR","GenR","GLAKU","HealthyStart","HELIX","HELIX","LINA","MOBA1","MOBA2","PREDO","VIVA","VIVA")
    age <- c(0,7.4,0,0,9.8,12.4,0,7.7,0,0,0,0,0,0,7.8)
    age_school <- c(0,1,0,0,1,1,0,1,0,0,0,0,0,0,1)
    age_birth <- c(1,0,1,1,0,0,1,0,1,1,1,1,1,1,0)
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
n_birth <- rowSums(sleep.data[,.(n_alspac_birth, n_eden_birth, n_genr_birth,
                                 n_healthystart_birth, n_inma_birth, n_lina_birth,
                                 n_moba1_birth, n_moba2_birth,
                                 n_predo_birth, n_viva_birth)], na.rm = T) 
n_studies_birth <- rowSums(!is.na(sleep.data[,.(n_alspac_birth, n_eden_birth, n_genr_birth,
                                                n_healthystart_birth, n_inma_birth, n_lina_birth,
                                                n_moba1_birth, n_moba2_birth,
                                                n_predo_birth, n_viva_birth)]), na.rm = T) 
cpg_birth <- sleep.data$probeID[n_birth >= 1000 & n_studies_birth > 3]

n_school <- rowSums(sleep.data[,.(n_alspac_school, n_genr_school, n_glaku_school,
                                  n_helix_school, n_viva_school)], na.rm = T) 
n_studies_school <- rowSums(!is.na(sleep.data[,.(n_alspac_school, n_genr_school, n_glaku_school,
                                                 n_helix_school, n_viva_school)]), na.rm = T) 
cpg_school <- sleep.data$probeID[n_school >= 1000 & n_studies_school > 3]

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
save(meta_results.data, file = "results/meta_results_sleep.Rdata")

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

### Change patterns
## Effect size change (n=27600)
changing.data <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05,]
# Positive constant (n=17)
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
## Effect size change (n=19618)
changing.data <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05,]
# Positive constant (n=0)
dim(changing.data[(changing.data$estimate_birth_cat > 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$estimate_school_cat_rev > 0 & changing.data$pval_school_cat_rev < 0.05),])
# Positive to null (n=1556)
dim(changing.data[(changing.data$estimate_birth_cat > 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$pval_school_cat_rev >= 0.05),])
# Null to positive (n=4224)
dim(changing.data[(changing.data$pval_birth_cat >= 0.05) & 
                    (changing.data$estimate_school_cat_rev > 0 & changing.data$pval_school_cat_rev < 0.05),])
# Negative constant (n=2)
dim(changing.data[(changing.data$estimate_birth_cat < 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$estimate_school_cat_rev < 0 & changing.data$pval_school_cat_rev < 0.05),])
# Negative to null (n=2181)
dim(changing.data[(changing.data$estimate_birth_cat < 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$pval_school_cat_rev >= 0.05),])
# Null to negative (n=3315)
dim(changing.data[(changing.data$pval_birth_cat >= 0.05) & 
                    (changing.data$estimate_school_cat_rev < 0 & changing.data$pval_school_cat_rev < 0.05),])
# Positive to negative (n=126)
dim(changing.data[(changing.data$estimate_birth_cat > 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$estimate_school_cat_rev < 0 & changing.data$pval_school_cat_rev < 0.05),])
# Negative to positive (n=187)
dim(changing.data[(changing.data$estimate_birth_cat < 0 & changing.data$pval_birth_cat < 0.05) & 
                    (changing.data$estimate_school_cat_rev > 0 & changing.data$pval_school_cat_rev < 0.05),])

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
dim(meta_results_ind.data[p.adjust(meta_results_ind.data$pval_birth_Cat) < 0.05,])

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
violin_sleep.data <- data.frame(c(rep("sleep",862318)),c(rep("birth",431159),rep("school-age",431159)),c(meta_results.data$cpg,meta_results.data$cpg),c(meta_results.data$estimate_birth_cat_abs,meta_results.data$estimate_school_cat_rev_abs),c(meta_results.data$z_birth_cat_abs,meta_results.data$z_school_cat_rev_abs),c(meta_results.data$pval_school_delta_cat,meta_results.data$pval_school_delta_cat))
names(violin_sleep.data) <- c("outcome","age","cpg","estimate","z","pval_school_delta_cat")
violin_sleep.data$age <- as.factor(violin_sleep.data$age)
save(violin_sleep.data, file = "results/violin_sleep.Rdata")

###### Heterogeneity analyses
# Function to run meta-regression and estimate heterogeneity
meta_regression_het <- function(cpg) {
  tryCatch({
    # Select cpg to analyze
    probe.data <- sleep.data[sleep.data$probeID == cpg, ]
    
    # Define betas
    beta <- probe.data[,.(beta_alspac_birth, beta_alspac_school, beta_eden_birth, beta_genr_birth,
                          beta_genr_school, beta_glaku_school, beta_healthystart_birth, beta_helix_school,
                          beta_inma_birth, beta_lina_birth, beta_moba1_birth, beta_moba2_birth,
                          beta_predo_birth, beta_viva_birth, beta_viva_school)]
    beta <- as.numeric(beta)
    
    # Define SE
    se <- probe.data[,.(se_alspac_birth, se_alspac_school, se_eden_birth, se_genr_birth,
                        se_genr_school, se_glaku_school, se_healthystart_birth, se_helix_school,
                        se_inma_birth, se_lina_birth, se_moba1_birth, se_moba2_birth,
                        se_predo_birth, se_viva_birth, se_viva_school)]
    se <- as.numeric(se)
    
    
    cohort <- c("ALSPAC","ALSPAC","HELIX","GenR","GenR","GLAKU","HealthyStart","HELIX","HELIX","LINA","MOBA1","MOBA2","PREDO","VIVA","VIVA")
    age_school <- factor(c(0,1,0,0,1,1,0,1,0,0,0,0,0,0,1))
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
save(meta_results_het.data, file = "results/meta_results_het_sleep.Rdata")

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
