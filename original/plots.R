library("ggplot2")
library("introdataviz")
library(cowplot)

load("results/violin_adhd.Rdata")
load("results/violin_gpf.Rdata")
load("results/violin_sleep.Rdata")
load("results/violin_bmi.Rdata")
load("results/violin_asthma.Rdata")

violin.data <- rbind(violin_adhd.data,violin_gpf.data,violin_sleep.data,violin_bmi.data,violin_asthma.data)
violin.data$outcome <- factor(violin.data$outcome, levels = c("adhd","gpf","sleep","bmi","asthma"))
levels(violin.data$outcome) <- c("ADHD","GPF","Sleep","BMI", "Asthma")
levels(violin.data$age) <- c("Birth","Childhood")

violin_b <- ggplot(violin.data, aes(x=outcome, y=estimate, fill=age)) + 
  geom_split_violin() + scale_y_continuous(trans='sqrt') + 
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.5), alpha = 0.3) +
  ggtitle("Effect size", subtitle = "Full sample size") + theme_gray(base_size = 24) + ylab("Effect Size (|β|)") + xlab("Outcome") + guides(fill=guide_legend(title="Methylation \nAssessment"))
violin_b

violin_z <- ggplot(violin.data, aes(x=outcome, y=z, fill=age)) + 
  geom_split_violin() + scale_y_continuous(trans='sqrt') + 
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.5), alpha = 0.3) +
  ggtitle("Statistical significance", subtitle = "Full sample size") + theme_gray(base_size = 24) + ylab("Statistical significance (|Z|)") + xlab("Outcome") + guides(fill=guide_legend(title="Methylation \nAssessment"))
violin_z

### Equal n
load("results/violin_adhd_equal_n.Rdata")
load("results/violin_sleep_equal_n.Rdata")
load("results/violin_bmi_equal_n.Rdata")
load("results/violin_asthma_equal_n.Rdata")

violin_equal_n.data <- rbind(violin_adhd_equal_n.data,violin_gpf.data,violin_sleep_equal_n.data,violin_bmi_equal_n.data,violin_asthma_equal_n.data)
violin_equal_n.data$outcome <- factor(violin_equal_n.data$outcome, levels = c("adhd","gpf","sleep","bmi","asthma"))
violin_equal_n_sig.data <- violin_equal_n.data[violin_equal_n.data$pval_school_delta_cat < 0.05, ]

violin_b_equal_n <- ggplot(violin_equal_n.data, aes(x=outcome, y=estimate, fill=age)) + 
  geom_split_violin() + scale_y_continuous(trans='sqrt') + 
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.5), alpha = 0.3) +
  ggtitle("Effect size", subtitle = "Equal sample sizes birth and childhood") + theme_gray(base_size = 24) + ylab("Effect Size (|β|)") + xlab("Outcome") + guides(fill=guide_legend(title="Methylation \nAssessment"))
violin_b_equal_n

violin_z_equal_n <- ggplot(violin_equal_n.data, aes(x=outcome, y=z, fill=age)) + 
  geom_split_violin() + scale_y_continuous(trans='sqrt') + 
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.5), alpha = 0.3) +
  ggtitle("Statistical significance", subtitle = "Equal sample sizes birth and childhood") + theme_gray(base_size = 24) + ylab("Statistical significance (|Z|)") + xlab("Outcome") + guides(fill=guide_legend(title="Methylation \nAssessment"))
violin_z_equal_n

violin_b_z <- plot_grid(violin_b + theme(legend.position="none"),
                        violin_z + theme(legend.position="none"),
                        violin_b_equal_n + theme(legend.position="none"),
                        violin_z_equal_n + theme(legend.position="none"))

legend <- get_legend(
  # create some space to the left of the legend
  violin_b
)

violin_b_z_legend <- plot_grid(violin_b_z, legend, rel_widths = c(3, .4))
ggsave(violin_b_z_legend, file = "figures/violin.pdf", width = 6000, height = 4500, units = "px", device = cairo_pdf)

### Change mean coefficients by outcome
change_all.data <- data.frame(Outcome=c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma"),age=c("Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood"),estimate=c(1.03,1.39,1.23,1.50,0.97,1.59,0.77,1.10,2.70,2.94))
change_all.data$Outcome <- factor(change_all.data$Outcome, levels = c("ADHD","GPF","Sleep","BMI", "Asthma"))

change_all <- ggplot(change_all.data, aes(x=age, y=estimate, group=Outcome, color=Outcome)) + 
  geom_line(size=1.5) + geom_point(size=4, shape=17) + ylim(0.7,3.7) +
  ggtitle("Effect size", subtitle = "Full sample size") + theme_gray(base_size = 24) + ylab("Mean Effect Size (|β|)") + xlab("Methylation Assessment") + guides(fill=guide_legend(title="Methylation \nAssessment"))
change_all

change_z_all.data <- data.frame(Outcome=c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma"),age=c("Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood"),estimate=c(1.02,0.78,0.78,0.78,0.76,0.77,0.75,0.86,0.82,0.77))
change_z_all.data$Outcome <- factor(change_z_all.data$Outcome, levels = c("ADHD","GPF","Sleep","BMI", "Asthma"))

change_z_all <- ggplot(change_z_all.data, aes(x=age, y=estimate, group=Outcome, color=Outcome)) + 
  geom_line(size=1.5) + geom_point(size=4, shape=17) + ylim(0.7,1.2) +
  ggtitle("Statistical significance", subtitle = "Full sample size") + theme_gray(base_size = 24) + ylab("Mean Statistical Significance (|Z|)") + xlab("Methylation Assessment") + guides(fill=guide_legend(title="Methylation \nAssessment"))
change_z_all

change_equal_n.data <- data.frame(Outcome=c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma"),age=c("Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood"),estimate=c(1.17,1.42,1.23,1.50,1.11,1.66,1.12,1.10,2.89,3.62))
change_equal_n.data$Outcome <- factor(change_equal_n.data$Outcome, levels = c("ADHD","GPF","Sleep","BMI", "Asthma"))

change_equal_n <- ggplot(change_equal_n.data, aes(x=age, y=estimate, group=Outcome, color=Outcome)) + 
  geom_line(size=1.5) + geom_point(size=4, shape=17) + ylim(0.7,3.7) +
  ggtitle("Effect size", subtitle = "Equal sample sizes birth and childhood") + theme_gray(base_size = 24) + ylab("Mean Effect Size (|β|)") + xlab("Methylation Assessment") + guides(fill=guide_legend(title="Methylation \nAssessment"))
change_equal_n

change_z_equal_n.data <- data.frame(Outcome=c("ADHD","ADHD","GPF","GPF","Sleep","Sleep","BMI","BMI","Asthma","Asthma"),age=c("Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood","Birth", "Childhood"),estimate=c(1.05,0.79,0.78,0.78,0.76,0.76,0.77,0.87,0.84,0.77))
change_z_equal_n.data$Outcome <- factor(change_equal_n.data$Outcome, levels = c("ADHD","GPF","Sleep","BMI", "Asthma"))

change_z_equal_n <- ggplot(change_z_equal_n.data, aes(x=age, y=estimate, group=Outcome, color=Outcome)) + 
  geom_line(size=1.5) + geom_point(size=4, shape=17) + ylim(0.7,1.2) +
  ggtitle("Statistical significance", subtitle = "Equal sample sizes birth and childhood") + theme_gray(base_size = 24) + ylab("Mean Statistical Significance (|Z|)") + xlab("Methylation Assessment") + guides(fill=guide_legend(title="Methylation \nAssessment"))
change_z_equal_n

change_p_z <- plot_grid(change_all + theme(legend.position="none"),
                        change_z_all + theme(legend.position="none"),
                        change_equal_n + theme(legend.position="none"),
                        change_z_equal_n + theme(legend.position="none"))

legend <- get_legend(
  # create some space to the left of the legend
  change_all
)

change_p_z_legend <- plot_grid(change_p_z, legend, rel_widths = c(3, .4))
ggsave(change_p_z_legend, file = "figures/change.pdf", width = 6000, height = 4500, units = "px", device = cairo_pdf)

### Increasing/decreasing ratio by p-value threshold
ratio.data <- read.csv("results/ratio.csv")
ratio.data <- ratio.data[,1:7]
ratio_long.data <- data.frame(Outcome = c(rep("ADHD",6),rep("GPF",6),rep("Sleep",6),rep("BMI",6),rep("Asthma",6)),
                              ratio = c(t(ratio.data[1,2:7]),t(ratio.data[2,2:7]),t(ratio.data[3,2:7]),t(ratio.data[4,2:7]),t(ratio.data[5,2:7])),
                              logp = -log10(rep(c(1,0.5,0.05,0.01,0.001,0.0001),5)))
 

ratio <- ggplot(ratio_long.data,aes(y = ratio, x = logp, colour = Outcome)) +
  geom_line(size=1.5) + geom_point(size=4, shape=17) +
  ylab("Number of DNAm sites with increasing effect size/decreasing effect size") + xlab("-log10 p-value threshold for change in association from birth to childhood") +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) + theme_gray(base_size = 24) + 
  ggtitle("Effect size change ratio across different p-value thresholds")
ggsave(ratio, file = "figures/ratio.pdf", width = 6000, height = 4500, units = "px", device = cairo_pdf)

                     