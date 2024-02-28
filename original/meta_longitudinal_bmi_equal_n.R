library(data.table)
library(metafor)
library(pbmcapply)
library(readr)
library(dplyr)
library(ggplot2)

# Load summary data from all cohorts
bmi_alspac_birth.data <- fread("sumstats/BMI/childBMI_cont_61-120_adj_House_noBW_ALSPAC_2015-07-02_MA1801.txt")
bmi_alspac_school.data <- fread("sumstats/BMI/ALSPAC_childBMI_cross_cont_full_model5.csv")
bmi_bamse_school.data <- fread("sumstats/BMI/BAMSE_childBMI_cross_cont_full_model5.csv")
bmi_chs_birth.data <- fread("sumstats/BMI/ChildBMI_cont_60-120_adj_House_noBW_CHS__20150713_MA1801.txt")
bmi_chamacos_birth.data <- fread("sumstats/BMI/ChildBMI_cont_61-120_adj_House_noBW_CHAMACOS_20160302_MA1801.txt")
bmi_chamacos_school.data <- fread("sumstats/BMI/CHAMACOS_childBMI_cross_cont_full_model5.csv")
bmi_chop_school.data <- fread("sumstats/BMI/CHOP_childBMI_cross_cont_full_model5.csv")
bmi_domino_school.data <- fread("sumstats/BMI/DOMInO_cross_childBMI_cont_full_model5_23oct.csv")
bmi_genr_birth.data <- fread("sumstats/BMI/GenR_ChildBMI_61-120_model5_21-11-2016.txt")
bmi_genr_school.data <- fread("sumstats/BMI/GENR_childBMI_cross_model5_naexclusieBFmissings_14062019.csv")
bmi_gecko_birth.data <- fread("sumstats/BMI/GECKO_childBMIZ_60mo_House_noBW_MA1801.txt")
bmi_healthystart_birth.data <- fread("sumstats/BMI/ChildBMI_cont_24-60_adj_House_noBW_HStart_20160803_1212MA.txt")
bmi_helix_school.data <- fread("sumstats/BMI/HELIXeur_cross_childBMI_cont_full_model5.csv")
bmi_inma_birth.data <- fread("sumstats/BMI/Model_5ChildBMI_cont_61120_adj_House_noBW_INMA_20151116_MA1801.txt")
bmi_inma_school.data <- fread("sumstats/BMI/INMA_childBMI_cross_cont_full_model5_okt2018.csv")
bmi_iow_birth.data <- fread("sumstats/BMI/ChildBMI_cont_24-60_adj_House_noBW_IOW_20160920_MA0912.txt")
bmi_iow_school.data <- fread("sumstats/BMI/IOW_childBMI_cross_cont_full_model5.csv")
bmi_moba2_birth.data <- fread("sumstats/BMI/ChildBMI_cont_61-120_adj_House_noBW_MoBa2_20160129.txt")
bmi_nest_black_birth.data <- fread("sumstats/BMI/ChildBMI_adj_61-120_COV_House_noBW_NEST_black_20161104_MA1801.txt")
bmi_nest_white_birth.data <- fread("sumstats/BMI/ChildBMI_adj_61-120_COV_House_noBW_NEST_white_20161104_MA1801.txt")
bmi_predo_birth.data <- fread("sumstats/BMI/ChildBMI_cont_24-60_adj_House_noBW_PREDO_20161118.txt")
bmi_piama_school.data <- fread("sumstats/BMI/PIAMA_childBMI_cross_cont_full_model5.csv")
bmi_viva_birth.data <- fread("sumstats/BMI/ChildBMI_Cont_61to120_Adj_House_noBW_Viva_20160122.txt")
bmi_viva_school.data <- fread("sumstats/BMI/Viva_childBMI_cross_cont_full_model5.csv")

# Add missing n for IOW
bmi_iow_birth.data$N <- 58

# Rename columns
setnames(bmi_alspac_birth.data, c("BETA", "SE", "P_VAL", "N_for_probe"), c("beta_alspac_birth", "se_alspac_birth", "p_alspac_birth", "n_alspac_birth"))
setnames(bmi_alspac_school.data, c("BETA", "SE", "P_VAL", "N_for_probe"), c("beta_alspac_school", "se_alspac_school", "p_alspac_school", "n_alspac_school"))
setnames(bmi_bamse_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_bamse_school", "se_bamse_school", "p_bamse_school", "n_bamse_school"))
setnames(bmi_chamacos_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_chamacos_birth", "se_chamacos_birth", "p_chamacos_birth", "n_chamacos_birth"))
setnames(bmi_chamacos_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_chamacos_school", "se_chamacos_school", "p_chamacos_school", "n_chamacos_school"))
setnames(bmi_chop_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_chop_school", "se_chop_school", "p_chop_school", "n_chop_school"))
setnames(bmi_chs_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_chs_birth", "se_chs_birth", "p_chs_birth", "n_chs_birth"))
setnames(bmi_domino_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_domino_school", "se_domino_school", "p_domino_school", "n_domino_school"))
setnames(bmi_gecko_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_gecko_birth", "se_gecko_birth", "p_gecko_birth", "n_gecko_birth"))
setnames(bmi_genr_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_genr_birth", "se_genr_birth", "p_genr_birth", "n_genr_birth"))
setnames(bmi_genr_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_genr_school", "se_genr_school", "p_genr_school", "n_genr_school"))
setnames(bmi_healthystart_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_healthystart_birth", "se_healthystart_birth", "p_healthystart_birth", "n_healthystart_birth"))
setnames(bmi_helix_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_helix_school", "se_helix_school", "p_helix_school", "n_helix_school"))
setnames(bmi_inma_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_inma_birth", "se_inma_birth", "p_inma_birth", "n_inma_birth"))
setnames(bmi_inma_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_inma_school", "se_inma_school", "p_inma_school", "n_inma_school"))
setnames(bmi_iow_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_iow_birth", "se_iow_birth", "p_iow_birth", "n_iow_birth"))
setnames(bmi_iow_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_iow_school", "se_iow_school", "p_iow_school", "n_iow_school"))
setnames(bmi_moba2_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_moba2_birth", "se_moba2_birth", "p_moba2_birth", "n_moba2_birth"))
setnames(bmi_nest_black_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_nest_black_birth", "se_nest_black_birth", "p_nest_black_birth", "n_nest_black_birth"))
setnames(bmi_nest_white_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_nest_white_birth", "se_nest_white_birth", "p_nest_white_birth", "n_nest_white_birth"))
setnames(bmi_piama_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_piama_school", "se_piama_school", "p_piama_school", "n_piama_school"))
setnames(bmi_predo_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_predo_birth", "se_predo_birth", "p_predo_birth", "n_predo_birth"))
setnames(bmi_viva_birth.data, c("BETA", "SE", "P_VAL", "N"), c("beta_viva_birth", "se_viva_birth", "p_viva_birth", "n_viva_birth"))
setnames(bmi_viva_school.data, c("BETA", "SE", "P_VAL", "N"), c("beta_viva_school", "se_viva_school", "p_viva_school", "n_viva_school"))

# Merge all cohorts
cohorts.list <- list(bmi_alspac_birth.data, bmi_alspac_school.data, bmi_bamse_school.data, bmi_chamacos_birth.data,
                     bmi_chamacos_school.data, bmi_chop_school.data, bmi_chs_birth.data, bmi_domino_school.data,
                     bmi_gecko_birth.data, bmi_genr_birth.data, bmi_genr_school.data, bmi_healthystart_birth.data,
                     bmi_helix_school.data, bmi_inma_birth.data, bmi_inma_school.data, bmi_iow_birth.data, bmi_iow_school.data,
                     bmi_moba2_birth.data, bmi_nest_black_birth.data, bmi_nest_white_birth.data,
                     bmi_piama_school.data, bmi_predo_birth.data, bmi_viva_birth.data, bmi_viva_school.data)

bmi.data <- Reduce(function(x, y) merge(x, y, by="probeID", all=TRUE), cohorts.list)

# Function to run meta-regression
meta_regression <- function(cpg) {
  tryCatch({
    # Select cpg to analyze
    probe.data <- bmi.data[bmi.data$probeID == cpg, ]
    
    # Define betas
    beta <- probe.data[,.(beta_alspac_birth, beta_alspac_school, beta_bamse_school, beta_chamacos_birth,
                          beta_chamacos_school, beta_chop_school, beta_chs_birth, beta_domino_school,
                          beta_gecko_birth, beta_genr_birth, beta_genr_school, beta_healthystart_birth,
                          beta_helix_school, beta_inma_birth, beta_inma_school, beta_iow_birth, beta_iow_school,beta_moba2_birth, beta_nest_black_birth, beta_nest_white_birth,
                          beta_piama_school, beta_predo_birth, beta_viva_birth, beta_viva_school)]
    beta <- as.numeric(beta)
    
    # Define SE
    se <- probe.data[,.(se_alspac_birth, se_alspac_school, se_bamse_school, se_chamacos_birth,
                        se_chamacos_school, se_chop_school, se_chs_birth, se_domino_school,
                        se_gecko_birth, se_genr_birth, se_genr_school, se_healthystart_birth,
                        se_helix_school, se_inma_birth, se_inma_school, se_iow_birth, se_iow_school,
                        se_moba2_birth, se_nest_black_birth, se_nest_white_birth,
                        se_piama_school, se_predo_birth, se_viva_birth, se_viva_school)]
    se <- as.numeric(se)
    
    
    cohort <- c("ALSPAC","ALSPAC","BAMSE","CHAMACOS","CHAMACOS","CHOP","CHS","DOMINO","GECKO","GenR","GenR","HEALTHYSTART","HELIX","INMA","INMA","IOW","IOW","MOBA2","NEST_B","NEST_W","PIAMA","PREDO","VIVA","VIVA")
    age <- c(0,9.7,8.3,0,9.1,5.5,0,5.0,0,0,9.8,0,7.4,0,4.5,0,10,0,0,0,8.1,0,0,7.8)
    age_school <- c(0,1,1,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,1)
    age_birth <- c(1,0,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,1,1,1,0,1,1,0)
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
n_birth <- rowSums(bmi.data[,.(n_alspac_birth, n_chamacos_birth,
                               n_chs_birth,
                               n_gecko_birth, n_genr_birth, n_healthystart_birth, 
                               n_inma_birth, n_iow_birth,
                               n_moba2_birth, n_nest_black_birth, n_nest_white_birth,
                               n_predo_birth, n_viva_birth)], na.rm = T) 
n_studies_birth <- rowSums(!is.na(bmi.data[,.(n_alspac_birth, n_chamacos_birth,
                                              n_chs_birth,
                                              n_gecko_birth, n_genr_birth, n_healthystart_birth, 
                                              n_inma_birth, n_iow_birth,
                                              n_moba2_birth, n_nest_black_birth, n_nest_white_birth,
                                              n_predo_birth, n_viva_birth)]), na.rm = T) 
cpg_birth <- bmi.data$probeID[n_birth >= 1000 & n_studies_birth > 3]

n_school <- rowSums(bmi.data[,.(n_alspac_school, n_bamse_school,
                                n_chamacos_school, n_chop_school, n_domino_school,
                                n_genr_school,
                                n_helix_school, n_inma_school, n_iow_school,
                                n_piama_school, n_viva_school)], na.rm = T) 
n_studies_school <- rowSums(!is.na(bmi.data[,.(n_alspac_school, n_bamse_school,
                                               n_chamacos_school, n_chop_school, n_domino_school,
                                               n_genr_school,
                                               n_helix_school, n_inma_school, n_iow_school,
                                               n_piama_school, n_viva_school)]), na.rm = T) 
cpg_school <- bmi.data$probeID[n_school >= 1000 & n_studies_school > 3]

cpg_qc <- cpg_birth[cpg_birth %in% cpg_school]

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
save(meta_results.data, file = "results/meta_results_bmi_equal_n.Rdata")

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
violin_bmi_equal_n.data <- data.frame(c(rep("bmi",871296)),c(rep("birth",435648),rep("school-age",435648)),c(meta_results.data$cpg,meta_results.data$cpg),c(meta_results.data$estimate_birth_abs,meta_results.data$estimate_school_cat_rev_abs),c(meta_results.data$z_birth_cat_abs,meta_results.data$z_school_cat_rev_abs),c(meta_results.data$pval_school_delta_cat,meta_results.data$pval_school_delta_cat))
names(violin_bmi_equal_n.data) <- c("outcome","age","cpg","estimate","z","pval_school_delta_cat")
violin_bmi_equal_n.data$age <- as.factor(violin_bmi_equal_n.data$age)
save(violin_bmi_equal_n.data, file = "results/violin_bmi_equal_n.Rdata")

