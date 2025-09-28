# Rebecca Satterwhite
# plate 3 was dropped bc the machine failed



# plot GCs
########

library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(growthcurver)
library(purrr)



# # # # # plate map set-up - decode randomized variables
map = read.csv("...platemap.csv")  
data1 = read.csv("...parsed_p1.csv")
data2 = read.csv("...parsed_p2.csv")
data4 = read.csv("...parsed_p4.csv")
data5 = read.csv("...parsed_p5.csv")

# Convert the "time" column from seconds to hours
data1$Time <- data1$Time / 60 / 60
data2$Time <- data2$Time / 60 / 60
data4$Time <- data4$Time / 60 / 60
data5$Time <- data5$Time / 60 / 60

# # # # # pair raw data to platemap & combine to 1 dataframe
data1$Temp = NULL
data1$Cycle = NULL
data2$Temp = NULL
data2$Cycle = NULL
data4$Temp = NULL
data4$Cycle = NULL
data5$Temp = NULL
data5$Cycle = NULL

reshaped1 <- reshape2::melt(data1, id=c("Time"), variable.name="Well", value.name="OD600")
reshaped2 <- reshape2::melt(data2, id=c("Time"), variable.name="Well", value.name="OD600")
reshaped4 <- reshape2::melt(data4, id=c("Time"), variable.name="Well", value.name="OD600")
reshaped5 <- reshape2::melt(data5, id=c("Time"), variable.name="Well", value.name="OD600")

annotated1 <- inner_join(reshaped1, map, by="Well")
annotated2 <- inner_join(reshaped2, map, by="Well")
annotated4 <- inner_join(reshaped4, map, by="Well")
annotated5 <- inner_join(reshaped5, map, by="Well")

annotated1$plate=1
annotated2$plate=2
annotated4$plate=4
annotated5$plate=5

annotated2$Time=annotated1$Time # make time points per plate identical for nice plotting
annotated4$Time=annotated1$Time 
annotated5$Time=annotated1$Time 

annotated = rbind(annotated1, annotated2,annotated4,annotated5)
annotated = subset(annotated, subset = Time <172801) # # # # # make the end clean at 48 hr 
dim(annotated) # 74112     8

# # # # ready for plotting
plotting_data = annotated
dim(plotting_data)  # 74112     8

# # # # # normalize by average NC
temp = subset(plotting_data, subset = Time == 0.000)
norm = mean(temp$OD600) # 0.09656276
var(temp$OD600)         # 2.83639e-06 confirm this is very small
plotting_data$OD600 = plotting_data$OD600 - norm

# # # # # remove negative controls
plotting_data = subset(plotting_data, subset = Strain != "blank")
dim(plotting_data) # 46320     8

# # # # # calculate stats! averaged across plates
stats <- plotting_data %>% dplyr::group_by(Strain, Time, Env) %>%
  dplyr::summarise(N=length(OD600), Average=mean(OD600), StDev=sd(OD600))
names(stats) =c("Pathogen","Time","Env","N","Average","StDev")
dim(stats) # 1930    6

# # # # # set up for plotting 
stats <- stats %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
stats <- stats %>% dplyr::mutate(Env = dplyr::recode(Env, "kb"="KB","lb"="LB"))
stats$Pathogen <- factor(stats$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))

# # # # # Fig: mean and error bars
pdf(file = ".../f_curves_10panels.pdf",width=5,height=4)
ggplot(data = stats, aes(x = Time, y = Average, color = Pathogen, group = Pathogen)) + 
  geom_ribbon(aes(ymin = Average - StDev, ymax = Average + StDev, fill = Pathogen), alpha = 0.1) + 
  geom_line(size = 0.5) + facet_grid(Env ~ Pathogen) + 
  scale_color_manual(name = "Pathogen",
    values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",
    values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  scale_x_continuous(limits = c(0, 48), expand = c(0, 0)) +
  theme_classic() + #theme(legend.position = "none") + 
  labs(title = " ", x = "Time (Hours)", y = "Absorbance at 600 nm")
dev.off()


# # # # # Fig: only media panels 
pdf(file = ".../f_curves_2panels.pdf",width=5,height=4)
ggplot(data = stats, aes(x = Time, y = Average, color = Pathogen, group = Pathogen)) + 
  geom_ribbon(aes(ymin = Average - StDev, ymax = Average + StDev, fill = Pathogen), alpha = 0.1) + 
  geom_line(size = 0.5) + facet_grid(~Env) + 
  scale_color_manual(name = "Pathogen",
                     values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",
                    values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  scale_x_continuous(limits = c(0, 48), expand = c(0, 0)) +
  theme_classic() + #theme(legend.position = "none") + 
  labs(title = " ", x = "Time (Hours)", y = "Absorbance at 600 nm")
dev.off()

########

# plot K
########

dim(annotated) # 74112     8

# drop blanks & rep 2
annotated = subset(annotated, subset = Strain != "blank")
annotated = subset(annotated, subset = Bio_rep==1 & plate !=2)
dim(annotated) # 17370     8

# run Growthcurver
results <- annotated %>%
  group_by(Strain, plate, Env) %>% group_split() %>%
  map_dfr(function(df_group) {sg <- SummarizeGrowth(df_group$Time, df_group$OD600)
  vals <- sg$vals
  tibble(Strain=unique(df_group$Strain), plate = unique(df_group$plate), 
         Env=unique(df_group$Env), k=vals$k, r=vals$r,
         gen_time=vals$t_gen, auc_l=vals$auc_l, note=vals$note, DT=vals$DT)})
dim(results) # 30  8
unique(results$note) # confirm no notes/fits worked fine
head(results) # Strain plate Env       k     r gen_time auc_l note 
names(results) = c("Pathogen", "Plate", "Env", "K", "r", "gen", "auc", "note")
results <- results %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
results <- results %>% dplyr::mutate(Env = dplyr::recode(Env, "kb"="KB","lb"="LB"))
results$Pathogen <- factor(results$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))

# # # # # # # # # # t-tests
kb_dc = subset(results, subset = Strain =="DC" & Env =="kb")
kb_14 = subset(results, subset = Strain =="A14" & Env =="kb")
t.test(kb_dc$k, kb_14$k) # p-value = 0.3207 true difference in means is equal to 0

lb_dc = subset(results, subset = Strain =="DC" & Env =="lb")
lb_14 = subset(results, subset = Strain =="A14" & Env =="lb")
t.test(lb_14$k, lb_dc$k, alternative = "less") # p-value = 0.0403 14 is sig less than DC in LB

lb_es = subset(results, subset = Strain =="ES" & Env =="lb")
lb_14 = subset(results, subset = Strain =="A14" & Env =="lb")
t.test(lb_es$k, lb_14$k) # p-value = 0.2131 true difference in means is equal to 0


# average across bioreps and get error
stats <- results %>% group_by(Pathogen, Env) %>%
  summarise(N=length(K), k_mean=mean(K), k_sd=sd(K), r_mean =mean(r), r_sd = sd(r))
dim(stats) # 10  7

names(stats) = c("Pathogen", "Env", "N","k_mean", "k_sd", "r_mean", "r_sd")
stats <- stats %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
stats <- stats %>% dplyr::mutate(Env = dplyr::recode(Env, "kb"="KB","lb"="LB"))
stats$Pathogen <- factor(stats$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))

# # # # # set up for plotting 
# # # points
pdf(file = ".../f_invitro_rk.pdf",width=5,height=4)
pk = ggplot(data=stats, aes(x=Env, y=k_mean, color=Pathogen)) + 
  facet_grid(~Pathogen) + geom_point(size=2) + 
  scale_color_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  guides(alpha = "none") + theme_classic() + geom_errorbar(width=0.8, aes(ymin=(k_mean - 0.5*k_sd), ymax=(k_mean + 0.5*k_sd)))+
  geom_jitter(data = results, aes(x = Env, y = K, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  labs(x="Pathogen",y="Carrying Capacity") +
  theme_minimal()+ 
   theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_text(size=10), 
        strip.text.x=element_text(size=11), strip.text.y=element_blank()) + 
  theme(legend.position = "none")

pr = ggplot(data=stats, aes(x=Env, y=r_mean, color=Pathogen)) + 
  facet_grid(~Pathogen) + geom_point(size=2) + 
  scale_color_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  guides(alpha = "none") + theme_classic() + geom_errorbar(width=0.8, aes(ymin=(r_mean - (0.5*r_sd)), ymax=(r_mean + (0.5*r_sd))))+
  geom_jitter(data = results, aes(x = Env, y = r, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
   labs(x="Environment",y="Growth Rate") +
  theme_minimal()+ 
  
  theme(axis.text.x = element_text(size=9, angle = 45, hjust = 1),
        axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), 
        strip.text.x=element_blank(), strip.text.y=element_blank()) + 
  theme(legend.position = "none")

combined_plot <- (pk / pr) + plot_annotation(tag_levels = 'A')
combined_plot
dev.off()




# # # box plot
names(results) =c("Pathogen","plate", "Env","k", "r", "gen_time", "auc_l", "note")
results <- results %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
results <- results %>% dplyr::mutate(Env = dplyr::recode(Env, "kb"="KB","lb"="LB"))
results$Pathogen <- factor(results$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))
results$Env <- as.factor(results$Env)

pdf(file = ".../f_k-box-2panels.pdf",width=4,height=3.5)
ggplot(results, aes(x = Pathogen, y = k, fill = Pathogen)) +
  geom_boxplot(size = 0.5,outlier.size = 0.5, width = 0.8)+#, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Env) +
  labs(y = "Carrying Capacity", x = "Environmnet", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_x_discrete(limits = c("A9", "DC3000", "ES4326", "1448A", "NP29")) +
  scale_y_continuous(labels = function(x) signif(x, 1)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size=9, angle = 45, hjust = 1),
        axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), 
        strip.text.x=element_text(size=8), strip.text.y=element_blank()) + 
  theme(legend.position = "none")
dev.off()


# # # box plot - 1 plot with envs next to each other
pdf(file = "...f_k-r-boxes-1panel.pdf",width=4,height=5)
kplot=ggplot(results, aes(x = Env, y = k, color=Pathogen, fill=Pathogen)) +
  geom_boxplot(alpha=0.3, size = 1, outlier.size = 0.5, width = 0.9) +
  facet_grid(~ Pathogen) + labs(y = "Carrying capacity", x = "Environment") +
  scale_color_manual(values = c("A9" = "sienna1", "DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(values = c("A9" = "sienna1", "DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  theme_minimal()+ theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10),legend.position = "none")

rplot=ggplot(results, aes(x = Env, y = r, color=Pathogen, fill=Pathogen)) +
  geom_boxplot(alpha=0.3, size = 1, outlier.size = 0.5, width = 0.9) +
  facet_grid(~ Pathogen) + labs(y = "Growth rate", x = "Environment") +
  scale_color_manual(values = c("A9" = "sienna1", "DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(values = c("A9" = "sienna1", "DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  theme_minimal()+ theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
        axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),legend.position = "none")

combined_plot <- (kplot / rplot) + plot_annotation(tag_levels = 'A')
combined_plot
dev.off()

# # # # # make k & r into one ggplot
results_long <- results %>% pivot_longer(cols = c(k, r), names_to = "Parameter", values_to = "Value")
color_values <- c("A9" = "sienna1", "DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")

pdf(file = ".../f_k-r-boxes-comb.pdf",width=4,height=5)
ggplot(results_long, aes(x = Env, y = Value, color = Pathogen, fill = Pathogen)) +
  geom_boxplot(alpha = 0.3, size = 1, outlier.size = 0.5, width = 0.9) +
  facet_grid(Parameter ~ Pathogen, scales = "free_y") +  # one row per 'k' and 'r'
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  labs(x = "Environment", y = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 10),legend.position = "none")
dev.off()
########


# # # # # # # Formally compare K across strains
library(car)
library(emmeans)


# # # check assumptions 
# for Levene, H0 = variances are equal
# for Shapiro, H0 = data are normally distributed
leveneTest(K ~ Pathogen, data = results)     # 0.001191 ** var not equal
shapiro.test(results$K)                      #  0.07362 data are normally distributed

# # # # # # # # # 2 way ANCOVA K
#######
mod_k <- lm(K ~ Pathogen * Env + Plate , data = results)
anova(mod_k)
# Response: k
#               Df Sum Sq Mean Sq  F value    Pr(>F)    
# Pathogen      4 0.12313 0.03078  12.5335 3.699e-05 ***
# Env           1 0.82028 0.82028 333.9788 1.630e-13 ***
# plate         1 0.00009 0.00009   0.0362  0.851116    
# Pathogen:Env  4 0.06499 0.01625   6.6154  0.001642 ** 
# Residuals    19 0.04667 0.00246      


# # # Tukey post hoc tests
emmeans(mod_k, pairwise ~ Pathogen | Env, adjust = "tukey") 
# $emmeans
# Env = KB:
#   Pathogen emmean     SE df lower.CL upper.CL
# A9        0.790 0.0286 19    0.731    0.850
# DC3000    0.702 0.0286 19    0.642    0.762
# ES4326    0.762 0.0286 19    0.702    0.822
# 1448A     0.721 0.0286 19    0.661    0.781
# NP29      0.864 0.0286 19    0.804    0.923
# 
# Env = LB:
#   Pathogen emmean     SE df lower.CL upper.CL
# A9        0.500 0.0286 19    0.440    0.560
# DC3000    0.518 0.0286 19    0.458    0.578
# ES4326    0.296 0.0286 19    0.236    0.356
# 1448A     0.347 0.0286 19    0.287    0.407
# NP29      0.526 0.0286 19    0.466    0.586
# 
# Confidence level used: 0.95 
# 
# $contrasts
# Env = KB:
#   contrast        estimate     SE df t.ratio p.value
# A9 - DC3000      0.08805 0.0405 19   2.176  0.2308
# A9 - ES4326      0.02822 0.0405 19   0.697  0.9546
# A9 - 1448A       0.06930 0.0405 19   1.713  0.4501
# A9 - NP29       -0.07309 0.0405 19  -1.806  0.3987
# DC3000 - ES4326 -0.05982 0.0405 19  -1.478  0.5878
# DC3000 - 1448A  -0.01875 0.0405 19  -0.463  0.9898
# DC3000 - NP29   -0.16114 0.0405 19  -3.982  0.0063**
# ES4326 - 1448A   0.04108 0.0405 19   1.015  0.8453
# ES4326 - NP29   -0.10132 0.0405 19  -2.504  0.1315
# 1448A - NP29    -0.14239 0.0405 19  -3.519  0.0172**
# 
# Env = LB:
#   contrast        estimate     SE df t.ratio p.value
# A9 - DC3000     -0.01789 0.0405 19  -0.442  0.9914
# A9 - ES4326      0.20410 0.0405 19   5.044  0.0006**
# A9 - 1448A       0.15316 0.0405 19   3.785  0.0097**
# A9 - NP29       -0.02642 0.0405 19  -0.653  0.9640
# DC3000 - ES4326  0.22199 0.0405 19   5.486  0.0002**
# DC3000 - 1448A   0.17105 0.0405 19   4.227  0.0037**
# DC3000 - NP29   -0.00853 0.0405 19  -0.211  0.9995
# ES4326 - 1448A  -0.05094 0.0405 19  -1.259  0.7182
# ES4326 - NP29   -0.23052 0.0405 19  -5.697  0.0002**
# 1448A - NP29    -0.17958 0.0405 19  -4.438  0.0023**



emmeans(mod_k, pairwise ~ Env, adjust = "tukey") # overall main effect of Env
# $emmeans
# Env emmean     SE df lower.CL upper.CL
# KB   0.768 0.0128 19    0.741    0.795
# LB   0.437 0.0128 19    0.410    0.464
# 
# Results are averaged over the levels of: Pathogen 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate     SE df t.ratio p.value
# KB - LB     0.331 0.0181 19  18.275  <.0001
# 
# Results are averaged over the levels of: Pathogen 


emmeans(mod_k, pairwise ~ Pathogen, adjust = "tukey") # overall main effect of pathogen
# $emmeans
# Pathogen emmean     SE df lower.CL upper.CL
# A9        0.645 0.0202 19    0.603    0.687
# DC3000    0.610 0.0202 19    0.568    0.652
# ES4326    0.529 0.0202 19    0.487    0.571
# 1448A     0.534 0.0202 19    0.492    0.576
# NP29      0.695 0.0202 19    0.653    0.737
# 
# Results are averaged over the levels of: Env 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast        estimate     SE df t.ratio p.value
# A9 - DC3000      0.03508 0.0286 19   1.226  0.7369
# A9 - ES4326      0.11616 0.0286 19   4.060  0.0053 **
# A9 - 1448A       0.11123 0.0286 19   3.887  0.0078 **
# A9 - NP29       -0.04976 0.0286 19  -1.739  0.4354
# DC3000 - ES4326  0.08108 0.0286 19   2.834  0.0706
# DC3000 - 1448A   0.07615 0.0286 19   2.661  0.0983
# DC3000 - NP29   -0.08483 0.0286 19  -2.965  0.0544
# ES4326 - 1448A  -0.00493 0.0286 19  -0.172  0.9998
# ES4326 - NP29   -0.16592 0.0286 19  -5.799  0.0001 **
# 1448A - NP29    -0.16099 0.0286 19  -5.626  0.0002 **
# Results are averaged over the levels of: Env 
# P value adjustment: tukey method for comparing a family of 5 estimates 

#######
# # # # # # # # # 2 way ANCOVA r

mod_r <- lm(r ~ Pathogen * Env + Plate , data = results)
anova(mod_r)
# Response: r
#                Df   Sum Sq  Mean Sq  F value    Pr(>F)    
#   Pathogen      4 0.085695 0.021424  18.1202 2.795e-06 ***
#   Env           1 0.138949 0.138949 117.5240 1.412e-09 ***
#   plate         1 0.009957 0.009957   8.4218 0.0091378 ** 
#   Pathogen:Env  4 0.037199 0.009300   7.8658 0.0006459 ***
#   Residuals    19 0.022464 0.001182      

# # # Tukey post hoc tests
emmeans(mod_r, pairwise ~ Pathogen | Env, adjust = "tukey") 
# Env = KB:
#   Pathogen emmean     SE df lower.CL upper.CL
# A9        0.116 0.0199 19   0.0747    0.158
# DC3000    0.169 0.0199 19   0.1272    0.210
# ES4326    0.152 0.0199 19   0.1105    0.194
# 1448A     0.172 0.0199 19   0.1306    0.214
# NP29      0.173 0.0199 19   0.1316    0.215
# 
# Env = LB:
#   Pathogen emmean     SE df lower.CL upper.CL
# A9        0.147 0.0199 19   0.1051    0.188
# DC3000    0.265 0.0199 19   0.2230    0.306
# ES4326    0.307 0.0199 19   0.2651    0.348
# 1448A     0.330 0.0199 19   0.2885    0.372
# NP29      0.415 0.0199 19   0.3734    0.457
# 
# Confidence level used: 0.95 
# 
# $contrasts
# Env = KB:
#   contrast        estimate     SE df t.ratio p.value
# A9 - DC3000     -0.05244 0.0281 19  -1.868  0.3667
# A9 - ES4326     -0.03573 0.0281 19  -1.273  0.7102
# A9 - 1448A      -0.05591 0.0281 19  -1.991  0.3072
# A9 - NP29       -0.05683 0.0281 19  -2.024  0.2925
# DC3000 - ES4326  0.01671 0.0281 19   0.595  0.9741
# DC3000 - 1448A  -0.00347 0.0281 19  -0.123  0.9999
# DC3000 - NP29   -0.00438 0.0281 19  -0.156  0.9999
# ES4326 - 1448A  -0.02018 0.0281 19  -0.719  0.9496
# ES4326 - NP29   -0.02110 0.0281 19  -0.751  0.9412
# 1448A - NP29    -0.00092 0.0281 19  -0.033  1.0000
# 
# Env = LB:
#   contrast        estimate     SE df t.ratio p.value
# A9 - DC3000     -0.11787 0.0281 19  -4.198  0.0039 ******A9 low
# A9 - ES4326     -0.16003 0.0281 19  -5.700  0.0001 ******A9 low
# A9 - 1448A      -0.18336 0.0281 19  -6.531  <.0001 ******A9 low
# A9 - NP29       -0.26829 0.0281 19  -9.556  <.0001 ****** NP high
# DC3000 - ES4326 -0.04216 0.0281 19  -1.502  0.5738
# DC3000 - 1448A  -0.06549 0.0281 19  -2.333  0.1777
# DC3000 - NP29   -0.15042 0.0281 19  -5.358  0.0003 ******NP high
# ES4326 - 1448A  -0.02334 0.0281 19  -0.831  0.9176
# ES4326 - NP29   -0.10827 0.0281 19  -3.856  0.0083 ******NP high
# 1448A - NP29    -0.08493 0.0281 19  -3.025  0.0482 ******NP high
# 
# P value adjustment: tukey method for comparing a family of 5 estimates 

emmeans(mod_r, pairwise ~ Env , adjust = "tukey") 
#  KB - LB    -0.136 0.0126 19 -10.841  <.0001

