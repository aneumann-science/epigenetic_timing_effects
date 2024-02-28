library(missMethyl)

### Load meta-analytic results and filter for nominally significant change
load("results/meta_results_adhd.Rdata")
change_cpg_adhd <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                       (meta_results.data$pval_birth_cat < 0.05 | 
                                          meta_results.data$pval_school_cat_rev < 0.05), "cpg"]
cpg_adhd <- meta_results.data$cpg

load("results/meta_results_gpf.Rdata")
change_cpg_gpf <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                      (meta_results.data$pval_birth_cat < 0.05 | 
                                         meta_results.data$pval_school_cat_rev < 0.05), "cpg"]
cpg_gpf <- meta_results.data$cpg

load("results/meta_results_sleep.Rdata")
change_cpg_sleep <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                        (meta_results.data$pval_birth_cat < 0.05 | 
                                           meta_results.data$pval_school_cat_rev < 0.05), "cpg"]
cpg_sleep <- meta_results.data$cpg

load("results/meta_results_bmi.Rdata")
change_cpg_bmi <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                      (meta_results.data$pval_birth_cat < 0.05 | 
                                         meta_results.data$pval_school_cat_rev < 0.05), "cpg"]
cpg_bmi <- meta_results.data$cpg

load("results/meta_results_asthma.Rdata")
change_cpg_asthma <- meta_results.data[meta_results.data$pval_school_delta_cat < 0.05 & 
                                         (meta_results.data$pval_birth_cat < 0.05 | 
                                            meta_results.data$pval_school_cat_rev < 0.05), "cpg"]
cpg_asthma <- meta_results.data$cpg

change_cpg <- c(change_cpg_adhd,change_cpg_gpf,change_cpg_sleep,change_cpg_bmi,change_cpg_asthma)
change_cpg <- change_cpg[!duplicated(change_cpg)]

cpg <- c(cpg_adhd,cpg_gpf,cpg_sleep,cpg_bmi,cpg_asthma)
cpg <- cpg[!duplicated(cpg)]

enrichment.data <- gometh(change_cpg, all.cpg = cpg)
write.csv(enrichment.data, "results/enrichment.csv")

enrichment_adhd.data <- gometh(change_cpg_adhd, all.cpg = cpg_adhd)
enrichment_adhd.data <- enrichment_adhd.data[order(enrichment_adhd.data$P.DE), ]
write.csv(enrichment_adhd.data, "results/enrichment_adhd.csv")

enrichment_gpf.data <- gometh(change_cpg_gpf, all.cpg = cpg_gpf)
enrichment_gpf.data <- enrichment_gpf.data[order(enrichment_gpf.data$P.DE), ]
write.csv(enrichment_gpf.data, "results/enrichment_gpf.csv")

enrichment_sleep.data <- gometh(change_cpg_sleep, all.cpg = cpg_sleep)
enrichment_sleep.data <- enrichment_sleep.data[order(enrichment_sleep.data$P.DE), ]
write.csv(enrichment_sleep.data, "results/enrichment_sleep.csv")

enrichment_bmi.data <- gometh(change_cpg_bmi, all.cpg = cpg_bmi)
enrichment_bmi.data <- enrichment_bmi.data[order(enrichment_bmi.data$P.DE), ]
write.csv(enrichment_bmi.data, "results/enrichment_bmi.csv")

enrichment_asthma.data <- gometh(change_cpg_asthma, all.cpg = cpg_adhd)
enrichment_asthma.data <- enrichment_asthma.data[order(enrichment_asthma.data$P.DE), ]
write.csv(enrichment_asthma.data, "results/enrichment_asthma.csv")
