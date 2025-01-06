library(meta)
library(grid)
library(readxl)

#---- AUC VL in placebo---#
master_file = "~/RSV_HIC_meta_analysis/database_for_meta_analysis.xlsx"
auc_VL_pbo <- read_excel(path = master_file, sheet = "auc_VL_pbo_after_first dose")

auc_VL_pbo.model <- metamean(
  n = n,
  mean = mean,
  sd = sd,
  studlab = Study,
  data = auc_VL_pbo,
  sm = "MRAW",
  title = "auc_VL_pbo",
  control = list(stepadj=0.5, maxiter=1000),
  common = TRUE,
  verbose = TRUE,
  test = "knha", method.tau = "REML", method.random.ci = "HK", prediction = T
)

auc_VL_pbo.plot <- forest.meta(auc_VL_pbo.model,
                               leftcols = c("studlab", "n"),
                               leftlabs = c("Study", "N"),
                               rightcols = c("effect", "sd", "ci", "w.common", "w.random"),
                               rightlabs = c("Mean (in number)", "SD", "95%-CI", "Weight \n (common)", "Weight \n (random)"),
                               print.tau2 = FALSE,
                               digits.sd = 2,
                               smlab = "Mean (visual)",
                               header.line  = "", fontsize = 10, spacing = 2, just = "center", print.Q = T, print.pval.Q = T, just.addcols = "center", colgap.forest.left = "2cm", colgap.forest.right = "2cm")

grid.text("Mean viral load AUC (log10 PFUe.hr/mL) after first dose for placebo - human viral challenge model", .5, .9, gp=gpar(cex=1.5, fontface = "bold"), just = "centre", vjust = -1)
#funnel(auc_VL_pbo.model, random = TRUE, studlab = T, xlab = "Viral load AUC (log10 PFUe.hr/mL) after first dose for placebo")


#-- Delta method for SE calculation---#
treatment_diff = read_excel(path = master_file, sheet = "auc_VL_high_dose")
A = 1/treatment_diff$mean_pbo
B = -treatment_diff$mean_trt/treatment_diff$mean_pbo^2
sigma_T_square = treatment_diff$SD_trt^2
sigma_P_square = treatment_diff$SD_pbo^2

n_trt = treatment_diff$n_trt
n_pbo = treatment_diff$n_pbo
n = n_trt + n_pbo

SD_trt_diff_RR_perc = sqrt((A^2*(n/n_trt)*sigma_T_square + B^2* (n/n_pbo)*sigma_P_square)*10000) 
SD_trt_diff_RR = sqrt((A^2*(n/n_trt)*sigma_T_square + B^2* (n/n_pbo) *sigma_P_square)) 

#---- Treatment difference meta-analysis, where SE are calculated using delta method ---#
treatment_diff$SD_delta_method = round(treatment_diff$SD_delta_method, 2)
#---- RMR of AUC VL---#
RR_percent.model <- metamean(
  n = n,
  mean = trt_diffpc,
  sd = SD_delta_method,
  studlab = Study,
  data = treatment_diff,
  sm = "MRAW",
  test = "knha", method.tau = "REML", method.random.ci = "HK"
)

RR_percent.plot <- forest.meta(RR_percent.model,
                               leftcols = c("studlab", "n_pbo", "n_trt", "n", "mean_pbo", "SD_pbo", "mean_trt", "SD_trt"),
                               leftlabs = c("Study", "N \n (Placebo)", "N \n (Drug)", "N \n (Total)", "Mean \n (Placebo)", "SD \n (Placebo)",  "Mean \n (Drug)", "SD \n (Drug)"),
                               rightcols = c("effect", "sd", "ci", "w.common", "w.random"),
                               rightlabs = c("Relative Mean Reduction (in number)", "SD", "95%-CI", "Weight \n (common)", "Weight \n (random)"),
                               print.tau2 = FALSE,
                               digits.sd = 2,
                               smlab = "Relative Mean Reduction (Visual)",
                               header.line  = "", fontsize = 10, spacing = 2, just = "center", print.Q = T, print.pval.Q = T, just.addcols = "center", colgap.forest.left = "1cm", colgap.forest.right = "1cm")

grid.text("Relative reduction in mean viral load AUC (log10 PFUe.hr/mL)  - Human Viral Challenge Model", x = .5, y = .9, gp=gpar(cex=1.5, fontface = "bold"), just = "centre", vjust = -1)

#funnel(RR_percent.model, random = TRUE, studlab = T, xlab = "Relative reduction in mean viral load AUC (log10 PFUe.hr/mL)")


#--- mean VL at peak for placebo ---#
mean_peak_VL_pbo <- read_excel(path = master_file, sheet = "mean_peak_VL")
FreqFit(means = mean_peak_VL_pbo$mean, sds = mean_peak_VL_pbo$sd, n = mean_peak_VL_pbo$n)
EmpFit(means = mean_peak_VL_pbo$mean, sds = mean_peak_VL_pbo$sd, n = mean_peak_VL_pbo$n)  

mean_peak_VL_pbo.model <- metamean(
  n = n,
  mean = mean,
  sd = sd,
  studlab = Study,
  data = mean_peak_VL_pbo,
  sm = "MRAW",
  title = "mean_peak_VL_pbo",
  control = list(stepadj=0.5, maxiter=1000),
  common = TRUE,
  verbose = TRUE,
  test = "knha", method.tau = "REML", method.random.ci = "HK", prediction = T
)

mean_peak_VL_pbo.plot <- forest.meta(mean_peak_VL_pbo.model,
                                     leftcols = c("studlab", "n"),
                                     leftlabs = c("Study", "N"),
                                     rightcols = c("effect", "sd", "ci", "w.common", "w.random"),
                                     rightlabs = c("Mean (in number)", "SD", "95%-CI", "Weight \n (common)", "Weight \n (random)"),
                                     print.tau2 = FALSE,
                                     digits.sd = 2,
                                     smlab = "Mean (visual)",
                                     header.line  = "", fontsize = 10, spacing = 2, just = "center", print.Q = T, print.pval.Q = T, just.addcols = "center", colgap.forest.left = "2cm", colgap.forest.right = "2cm")

grid.text("Mean viral load at peak (log10 PFUe/mL) for placebo - human viral challenge model", .5, .9, gp=gpar(cex=1.5, fontface = "bold"), just = "centre", vjust = -1)
#funnel(auc_VL_pbo.model, random = TRUE, studlab = T, xlab = "Viral load AUC (log10 PFUe/mL) after first dose for placebo")


#---- AUC Total Symptom Score in placebo---#
auc_symsco_rr <- read_excel(master_file, sheet = "auc_total_symptom_score_pb")
auc_symsco_rr = auc_symsco_rr[c(1,3,7,8), ]

AUC_pbo.model <- metamean(
  n = n_pbo,
  mean = mean_pbo,
  sd = SD_pbo,
  studlab = Study,
  data = auc_symsco_rr,
  sm = "MRAW",
  test = "knha",
  method.tau = "REML",
  method.random.ci = "HK",
  title = "AUC Placebo",
  prediction = T
)

AUC_pbo.plot <- forest.meta(AUC_pbo.model,
                            leftcols = c("studlab", "n"),
                            leftlabs = c("Study", "N"),
                            rightcols = c("effect", "sd", "ci", "w.common", "w.random"),
                            rightlabs = c("Mean (in number)", "SD", "95%-CI", "Weight \n (common)", "Weight \n (random)"),
                            print.tau2 = FALSE,
                            digits.sd = 2,
                            smlab = "Mean (visual)",
                            header.line  = "", fontsize = 10, spacing = 2, just = "center", print.Q = T, print.pval.Q = T, just.addcols = "center", colgap.forest.left = "2cm", colgap.forest.right = "2cm")

grid.text("Mean Total Symptom Score AUC for Placebo - Human Viral Challenge Model", x = .5, y = .9, gp=gpar(cex=2, fontface = "bold"), just = "centre", vjust = -1)


#---- AUC Total Symptom Score RMR---#
auc_symsco_rr$Trt_diff_rr = (auc_symsco_rr$mean_pbo - auc_symsco_rr$mean_trt)/auc_symsco_rr$mean_pbo

A = 1/auc_symsco_rr$mean_pbo
B = -auc_symsco_rr$mean_trt/auc_symsco_rr$mean_pbo^2
sigma_T_square = auc_symsco_rr$SD_trt^2
sigma_P_square = auc_symsco_rr$SD_pbo^2

n_trt = auc_symsco_rr$n_trt
n_pbo = auc_symsco_rr$n_pbo
n = n_trt + n_pbo

SD_trt_diff_RR_perc = sqrt((A^2*(n/n_trt)*sigma_T_square + B^2* (n/n_pbo)*sigma_P_square)*10000) 
SD_trt_diff_RR = sqrt((A^2*(n/n_trt)*sigma_T_square + B^2* (n/n_pbo) *sigma_P_square)) 

auc_symsco_rr$n = n
RR_percent.model <- metamean(
  n = n,
  mean = Trt_diff_rr,
  sd = SD_trt_diff_RR,
  studlab = Study,
  data = auc_symsco_rr,
  sm = "MRAW",
  test = "knha",
  method.tau = "REML",
  method.random.ci = "HK",
  title = "Relative Reduction"
)

RR_percent.plot <- forest.meta(RR_percent.model,
                               leftcols = c("studlab", "n_pbo", "n_trt", "n", "mean_pbo", "SD_pbo", "mean_trt", "SD_trt"),
                               leftlabs = c("Study", "N \n (Placebo)", "N \n (Drug)", "N \n (Total)", "Mean \n (Placebo)", "SD \n (Placebo)",  "Mean \n (Drug)", "SD \n (Drug)"),
                               rightcols = c("effect", "sd", "ci", "w.common", "w.random"),
                               rightlabs = c("Relative Mean Reduction (in number)", "SD", "95%-CI", "Weight \n (common)", "Weight \n (random)"),
                               print.tau2 = FALSE,
                               digits.sd = 2,
                               smlab = "Relative Mean Reduction (Visual)",
                               header.line  = "", fontsize = 10, spacing = 2, just = "center", print.Q = T, print.pval.Q = T, just.addcols = "center", colgap.forest.left = "1cm", colgap.forest.right = "1cm")

grid.text("Relative Reduction in Mean Total Symptom Score AUC - Human Viral Challenge Model", x = .5, y = .9, gp=gpar(cex=2, fontface = "bold"), just = "centre", vjust = -1)







#---- DESCRIPTIVE STATISTICS -----#
#--- Time to mean peak viral load ---#
time_to_mean_peak_VL = read_excel(path = master_file, sheet = "Mean time to peak VL")
summary(na.omit(time_to_mean_peak_VL$`Time to mean VL peak`))
sd(na.omit(time_to_mean_peak_VL$`Time to mean VL peak`))/sqrt(length(na.omit( time_to_mean_peak_VL$`Time to mean VL peak`)))
t.test(time_to_mean_peak_VL$`Time to mean VL peak`)

#--- Time to mean peak total symptom score ---#
time_to_mean_peak_tss = read_excel(path = master_file, sheet = "time_to_symptompeak")
summary(na.omit(time_to_mean_peak_tss$`time to mean total symptom score peak`))
sd(na.omit(time_to_mean_peak_tss$`time to mean total symptom score peak`))/sqrt(length(na.omit(time_to_mean_peak_tss$`time to mean total symptom score peak`)))

t.test(time_to_mean_peak_tss$`time to mean total symptom score peak`)


