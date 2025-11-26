# 20251125 revised script for all figures generated in R
# For Microbial Ecology - WIDTH 84 mm/8.4 cm, or 174 mm/17.4 cm - HEIGHT up to 234 mm/23.4 cm.



# # # # # Set up paths for data import/output ### UPDATES REQUIRED HERE
#########

project_root <- ".../figures-revision1"    # <--- update this once # INPUT
output_dir <- ".../output"                 # <--- update this once # OUTPUT


# Paths for each figure
folder_fig2 <- file.path(project_root, "data_Fig2")
folder_fig3 <- file.path(project_root, "data_Fig3")
folder_fig4 <- file.path(project_root, "data_Fig4")
folder_fig6 <- file.path(project_root, "data_Fig6")
folder_fig7 <- file.path(project_root, "data_Fig7")
folder_figS1 <- file.path(project_root, "data_FigS1")
#########

# # # # # FIG 1 - schematic pdf, biorender

# # # # # FIG 2 - IN VITRO width = 17.4, height = 10
#########

# # # # # # import R packages and set a folder for analysis
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(growthcurver)


# # # # # plate map set-up - decode randomized variables
in_file <- "aug_platemap_master.csv"

map14 <- read.csv(file.path(folder_fig2, in_file))
map15 <- read.csv(file.path(folder_fig2, in_file))
map16 <- read.csv(file.path(folder_fig2, in_file))
map14$Date="aug14"
map15$Date="aug15"
map16$Date="aug16"

map14$Env = recode(map14$Env,"media1"="AB","media2"="KB","media3"="LB")
map15$Env = recode(map15$Env,"media1"="LB","media2"="AB","media3"="KB")
map16$Env = recode(map16$Env,"media1"="KB","media2"="LB","media3"="AB")

map14$ID = recode(map14$ID,"1"="lib_14","2"="wt_np","3"="lib_es","4"="wt_14","5"="lib_a9","6"="wt_a9","7"="wt_es","8"="lib_np","9"="lib_dc","10"= "wt_dc")
map15$ID = recode(map15$ID,"1"="lib_14","2"="wt_np","3"="lib_es","4"="wt_14","5"="lib_a9","6"="wt_a9","7"="wt_es","8"="lib_np","9"="lib_dc","10"= "wt_dc")
map16$ID = recode(map16$ID,"1"="lib_es","2"="lib_dc","3"="lib_np","4"="wt_es","5"="wt_np","6"="wt_dc","7"="lib_14","8"="lib_a9","9"="wt_a9","10"= "wt_14")

# # # # # import raw OD by time data 
data14 <- read.csv(file.path(folder_fig2, "aug14data.csv"), header = TRUE)
data15 <- read.csv(file.path(folder_fig2, "aug15data.csv"), header = TRUE)
data16 <- read.csv(file.path(folder_fig2, "aug16data.csv"), header = TRUE)

# # # # # confirm that time points are identical
# if they are not, values will not plot correctly
data15$Time=data14$Time
data16$Time=data14$Time

# # # # # pair raw data to platemaps & combine to 1 dataframe
reshaped14 <- melt(data14, id=c("Time", "Temperature"), variable.name="Well", value.name="OD600")
annotated14 <- inner_join(reshaped14, map14, by="Well")
dim(annotated14) # 10464     8

reshaped15 <- melt(data15, id=c("Time", "Temperature"), variable.name="Well", value.name="OD600")
annotated15 <- inner_join(reshaped15, map15, by="Well")
dim(annotated15) # 10464     8

reshaped16 <- melt(data16, id=c("Time", "Temperature"), variable.name="Well", value.name="OD600")
annotated16 <- inner_join(reshaped16, map16, by="Well")
dim(annotated16) # 10464     8

full = rbind(annotated14,annotated15,annotated16)
dim(full) #  31392     8

# Convert the "time" column from seconds to hours
full$Time <- full$Time / 60 / 60
tail(full)

# # # # # remove water negative controls
wet = subset(full, subset = Env == "water")
full <- anti_join(full,wet)
dim(full) # 23544     8

# # # # # remove media negative controls
ncs = subset(full, subset = ID == "nc")
full <- anti_join(full,ncs)
dim(full) # 9810    8

wts = full[grep("wt_", full$ID),]
dim(wts) # 4905    8




# # # # # get stats
wts <- wts %>% dplyr::mutate(ID = dplyr::recode(ID, "wt_es"="ES4326","wt_np"="NP29","wt_14"="1448A","wt_a9"="A9", "wt_dc"= "DC3000"))

# run Growthcurver (here Date = plate)
results <- wts %>%
  group_by(ID, Env, Date) %>% group_split() %>%
  map_dfr(function(df_group) {sg <- SummarizeGrowth(df_group$Time, df_group$OD600)
  vals <- sg$vals
  tibble(Pathogen=unique(df_group$ID), Plate = unique(df_group$Date), 
         Env=unique(df_group$Env), K=vals$k, r=vals$r,
         gen_time=vals$t_gen, AUC=vals$auc_l, note=vals$note, DT=vals$DT)})
dim(results) # 45  8
#print(results, n=45)
unique(results$note) # confirm no notes/fits worked fine
results$Pathogen <- factor(results$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))
dim(results) # 45  8

# average across bioreps and get error
stats <- results %>% group_by(Pathogen, Env) %>%
  summarise(N=length(K), k_mean=mean(K), k_sd=sd(K), r_mean =mean(r), r_sd = sd(r))
dim(stats) # 15  7


# # # PLOT FIG 1 
pk = ggplot(data=stats, aes(x=Env, y=k_mean, color=Pathogen)) + 
  facet_grid(~Pathogen) + geom_point(size=3.5) + 
  scale_color_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  guides(alpha = "none") + theme_classic() + geom_errorbar(width=0.9, aes(ymin=(k_mean - 0.5*k_sd), ymax=(k_mean + 0.5*k_sd)))+
  geom_jitter(data = results, aes(x = Env, y = K, color = Pathogen), 
              size = 2, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  labs(x="",y="Carrying Capacity") +
  theme_minimal() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), axis.text.y = element_text(size=9),
      axis.title.x=element_blank(), axis.title.y=element_text(size=12),
      strip.text.x=element_blank(), strip.text.y=element_blank())

pr = ggplot(data=stats, aes(x=Env, y=r_mean, color=Pathogen)) + 
  facet_grid(~Pathogen) + geom_point(size=3.5) + 
  scale_color_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  scale_fill_manual(name = "Pathogen",values = c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet", "1448A" = "mediumaquamarine", "NP29" = "olivedrab3"), guide = "none")+  
  guides(alpha = "none") + theme_classic() + geom_errorbar(width=0.9, aes(ymin=(r_mean - (0.5*r_sd)), ymax=(r_mean + (0.5*r_sd))))+
  geom_jitter(data = results, aes(x = Env, y = r, color = Pathogen), 
              size = 2, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  labs(x="", y="Growth Rate") +
  theme_minimal() +  theme(legend.position = "none") + 
  theme(axis.text.x = element_blank(),axis.text.y=element_text(size=9),
        axis.title.x=element_blank(), axis.title.y=element_text(size=12), 
        strip.text.x=element_text(size=12), strip.text.y=element_blank())
  
combined_plot <- (pr / pk) + plot_annotation(tag_levels = 'a') 
 
ggsave(
  filename = file.path(output_dir, "Fig2-rev.eps"),
  plot = combined_plot,
  device = cairo_ps,
  width = 17.4,
  height = 10,
  units = "cm"
)

ggsave(
  filename = file.path(output_dir, "Fig2-rev.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 17.4,
  height = 10,
  units = "cm",
  dpi = 600
)

#########

# # # # # FIG 3 - IN PLANTA width = 17.4, height = 10
#########

#PIPELINE
# 1) import data
# 2) drop libs, drop controls
# 3) convert CT to Q
# 4) account for dilution factor 
# 5) change Q=NA to Q=0
# 6) find k (highest Q)
# 7) stats to get mean and error around k
# 8) ANOVA on k data frame (not stats data frame)
# 9) find r with linear regression

library(ggplot2)
library(tidyr)
library(patchwork)
library(emmeans)
library(car)
library(dplyr)  
library(broom)
library(MESS)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(emmeans)
library(car)
library(dplyr)  


# # # import raw data & put into 1 dataframe
data1 <- read.table(file.path(folder_fig3, "2_1_6.txt"),  header = TRUE)
data2 <- read.table(file.path(folder_fig3, "2_2_5a.txt"), header = TRUE)
data3 <- read.table(file.path(folder_fig3, "2_4_2.txt"),  header = TRUE)
data4 <- read.table(file.path(folder_fig3, "2_4_3.txt"),  header = TRUE)
data5 <- read.table(file.path(folder_fig3, "2_11_6.txt"), header = TRUE)
data6 <- read.table(file.path(folder_fig3, "2_13_4.txt"), header = TRUE)
data7 <- read.table(file.path(folder_fig3, "2_14_6.txt"), header = TRUE)
data8 <- read.table(file.path(folder_fig3, "2_15_2.txt"), header = TRUE)
data = rbind(data1,data2,data3,data4,data5,data6,data7,data8, stringsAsFactors = FALSE) 
dim(data) # 1720    9

# # # drop plate controls, plant cell counts, old treatments, & extra columns
data = subset(data, subset = target != "At_gene") # remove At gene datapoints
data = subset(data, subset = task != "NTC")
data = subset(data, subset = cv != "na")
data = subset(data, subset = pv != "na")
data = subset(data, subset = cv != "DERBY")
data = subset(data, subset = pv != "B7")
data = subset(data, subset = task != "STANDARD")
data = subset(data, select =-c(well,target,task))
dim(data) # 465   6


# # # account for dilution factor = ct*1/total dilution factor
# add 1- because ct is inversely proportional to quantity
# diluted DNA (1/2*2/8) = 0.125 = total dilution factor
data$ct = as.numeric(data$ct) # throws a warning as Undetermined go to NA
data$ct = data$ct*1-(1/0.1) 

# # # use standard curve correlation coefficients to calculate quantity as x = log10(titer)
b = 37.2
m = -3.5
data$x = NA
data$x = as.numeric(data$x)
data$x = (data$ct-b)/m

# convert Q=NA to Q=0
data_na <- data[is.na(data$ct), ]
data = anti_join(data, data_na)
data_na$x = 0
data_na$ct = 40
data = rbind(data, data_na)
dim(data) # 465   7
data <- data %>% dplyr::mutate(pv = dplyr::recode(pv, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
data <- data %>% dplyr::mutate(cv = dplyr::recode(cv, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))

# # # # subset Mocks at t=24 
mocks = subset(data, subset = pv=="MOCK" & time == 24 & x < 5)

# # # # find K (highest mean Q)
statsk <- data %>% group_by(cv, pv, time) %>% summarise(N = length(x),
                                                        K = mean(x),Ksd = sd(x), .groups = "drop_last") %>% slice_max(K, n = 1, with_ties = FALSE) 
statsk <- statsk %>% dplyr::mutate(pv = dplyr::recode(pv, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
statsk <- statsk %>% dplyr::mutate(cv = dplyr::recode(cv, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))
#write.csv(statsk, ".../statsk.csv")

dim(statsk) # 30  6
stats=statsk
stats$SE = stats$Ksd / sqrt(stats$N) # add SE
stats = stats[order(stats$cv),]
dim(stats) # 30  7

# # # # make full dataframe of K datapoints
k_data <- dplyr::semi_join(data, stats, by = c("cv", "pv", "time"))
dim(k_data) # 129   7
k_data$ct=NULL 
k_data$ran=NULL 
#head(k_data)
names(k_data)= c("Host", "Pathogen", "Time", "Plate", "Quantity")
k_data <- k_data %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
k_data <- k_data %>% dplyr::mutate(Host = dplyr::recode(Host, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))
k_data$Pathogen <- factor(k_data$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29", "MOCK")) # make them plot in order
k_data$Host <- factor(k_data$Host, levels = c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA"))  

# # # # combine K (final) with MOCKS t24 (initial)
# take average of the mocks
head(mocks) # ct     cv   pv time plate ran        x
mocks = subset(mocks, subset = plate !=3)
statsm <- mocks %>% group_by(cv) %>% summarise(Time = 0, Quantity = mean(x), .groups = "drop_last") 
names(statsm) = c("Host", "Time",  "Quantity")
statsm <- statsm %>% dplyr::mutate(Host = dplyr::recode(Host, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))
statsm$Plate = 9
head(statsm) # Host        Time     N Quantity   Qsd Plate

# add pathogens to avg mocks 
k_data = subset(k_data, subset = Pathogen != "MOCK")
head(k_data) # Host Pathogen Time Plate Quantity
host_pathogens <- k_data %>% select(Host, Pathogen) %>% distinct()
expanded <- host_pathogens %>% dplyr::left_join(statsm, by = "Host") %>% mutate(Time = 0) # throws a warning
head(expanded) # Host Pathogen Plate Time Quantity     
unique(expanded$Time) # for just time0

# add avg mocks to K data
data_all = rbind(expanded,k_data)
head(data_all)
tail(data_all)
unique(data_all$Pathogen)

# # # # TABLE S1
# # # # Fit linear model to log(CFU) ~ Time, extract slope (rate) & AUC
# The slope (coef[2]) is the estimated growth rate per day.library(MESS)
growth_summary <- data_all %>% group_by(Host, Pathogen) %>% arrange(Time) %>%
  summarise(model = list(lm(Quantity ~ Time)),  # store the model
            auc = auc(Time, Quantity, type = "linear"),
            .groups = "drop") %>%
  mutate(slope = map_dbl(model, ~ coef(.x)[2]),
         r_squared = map_dbl(model, ~ summary(.x)$r.squared),
         p_value = map_dbl(model, ~ coef(summary(.x))[2, 4]) ) %>% select(-model) 

growth_summary = growth_summary[order(growth_summary$Host, growth_summary$Pathogen),]
print(growth_summary, n=25)

# add K to this table
head(statsk)
statsk$N=NULL
statsk$Ksd=NULL
names(statsk) = c("Host", "Pathogen", "Time", "K")
head(growth_summary)
growth_summary2 <- growth_summary %>% left_join(select(statsk, Host, Pathogen, K) %>% distinct(), by = c("Host", "Pathogen"))
head(growth_summary2)
#write.csv(growth_summary2, ".../table_planta.csv")


# Step 1: Rename columns in mock_df to match k_data
mock_df <- mocks %>%
  rename(Host = cv,
         Pathogen = pv,
         Time = time,
         Plate = plate,
         Quantity = x) %>%
  mutate(Infection = "Mock")

k_data <- k_data %>%
  mutate(Infection = "Real")

# Step 2: Combine data frames
combined_df <- bind_rows(
  mock_df %>% select(Host, Pathogen, Time, Plate, Quantity, Infection),
  k_data %>% select(Host, Pathogen, Time, Plate, Quantity, Infection))

# Step 3: Statistical test: Mann-Whitney U test (Wilcoxon rank-sum test in R)
mock_quant <- combined_df %>% filter(Infection == "Mock") %>% pull(Quantity)
real_quant <- combined_df %>% filter(Infection == "Real") %>% pull(Quantity)
test_result <- wilcox.test(mock_quant, real_quant, alternative = "less") # 7.874e-13

# ttest is same as Mann-Whitney
shapiro.test(mock_quant) # normal
shapiro.test(real_quant)# normal
t.test(mock_quant, real_quant, alternative = "less") # < 2.2e-16

# # # add zero time point (statsm) to all data points (data)

# get stats for mocks and data
head(statsm)
statsm <- mocks %>% group_by(cv) %>% summarise(Time = 0, N=length(x), Quantity = mean(x), Qsd=sd(x), .groups = "drop_last") 
names(statsm) = c("Host", "Time", "N", "Quantity", "Qsd")
statsm <- statsm %>% dplyr::mutate(Host = dplyr::recode(Host, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))
head(statsm) # Host        Time     N Quantity   Qsd 

head(data)
statsall <- data %>% group_by(cv, pv, time) %>% summarise(N = length(x), Quantity = mean(x), Qsd = sd(x), .groups = "drop_last") 
head(statsall)
names(statsall) = c("Host","Pathogen", "Time", "N", "Quantity", "Qsd")

# # # combine mocks with data
# add pathogens to mock
expanded <- host_pathogens %>% dplyr::left_join(statsm, by = "Host") %>% mutate(Time = 0)
head(statsall)
head(expanded)
mat_dat = rbind(expanded, statsall)

mat_dat <- mat_dat %>% dplyr::mutate(Pathogen = dplyr::recode(Pathogen, "A14"="1448A","ES"="ES4326","DC"="DC3000","NP"="NP29"))
mat_dat <- mat_dat %>% dplyr::mutate(Host = dplyr::recode(Host, "AFRICAN"="BEAN","BRANDY"="TOMATO-DC","GOLDEN"="TOMATO-A9","RRS10"="A.THALIANA"))
mat_dat$Pathogen <- factor(mat_dat$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29", "MOCK")) # make them plot in order
mat_dat$Host <- factor(mat_dat$Host, levels = c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA"))  
mat_dat = subset(mat_dat, subset = Pathogen != "MOCK")

k_data= subset(k_data, subset = Pathogen != "MOCK")
k_data <- k_data %>% group_by(Host, Pathogen) %>%
  mutate(mean_estimate = mean(Quantity, na.rm = TRUE),
         diff_from_mean = abs(Quantity - mean_estimate)) %>%
  arrange(diff_from_mean) %>% slice_head(n = 4) %>% ungroup() %>%
  dplyr::select(-mean_estimate, -diff_from_mean) 

data_all=data_all[order(data_all$Host, data_all$Pathogen),]


# Extract rows where Time == 0 and Plate == 9
time0_plate9 <- subset(data_all, Time == 0 & Plate == 9)

# We want to replicate these rows for all plates except 9
plates <- unique(data_all$Plate)
plates_to_copy <- setdiff(plates, 9)

# Create a copy of these rows for each plate in plates_to_copy
new_rows <- do.call(rbind, lapply(plates_to_copy, function(p) {
  new_row <- time0_plate9
  new_row$Plate <- p
  new_row}))

data_all_new <- rbind(data_all, new_rows)
data_all_new = subset(data_all_new, subset = Plate !=9)
data_all_new=data_all_new[order(data_all_new$Host, data_all_new$Pathogen),]

auc_trap <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}


data_auc <- data_all_new %>% group_by(Host, Pathogen, Plate) %>% arrange(Time) %>% summarise(
  model = list(lm(Quantity ~ Time)), auc = auc_trap(Time, Quantity), .groups = "drop" ) %>%
  mutate( slope = map_dbl(model, ~ coef(.x)[["Time"]] %||% NA_real_),
          r_squared = map_dbl(model, ~ summary(.x)$r.squared %||% NA_real_),
          r = sqrt(r_squared) * sign(slope),
          p_value = map_dbl(model, ~ { coefs <- coef(summary(.x))
          if ("Time" %in% rownames(coefs)) coefs["Time", "Pr(>|t|)"] else NA_real_}) ) %>%
  select(Host, Pathogen, Plate, slope, auc, r_squared, r, p_value)


growth_summary <- data_all_new %>% group_by(Host, Pathogen, Plate) %>% arrange(Time) %>% summarise(
  model = list(lm(Quantity ~ Time)), auc = auc_trap(Time, Quantity), .groups = "drop" ) %>%
  mutate(slope = map_dbl(model, ~ coef(.x)[["Time"]] %||% NA_real_),
         r_squared = map_dbl(model, ~ summary(.x)$r.squared %||% NA_real_),
         r = sqrt(r_squared) * sign(slope),
         p_value = map_dbl(model, ~ { coefs <- coef(summary(.x))
         if ("Time" %in% rownames(coefs)) coefs["Time", "Pr(>|t|)"] else NA_real_}) ) %>%
  select(Host, Pathogen, Plate, slope, auc, r_squared, r, p_value)

growth_summary <- growth_summary[!is.na(growth_summary$slope), ]
dim(growth_summary) # 79  8 has raw estimates of rate and auc for stats

growth_summary$Pathogen <- factor(growth_summary$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29", "MOCK")) # make them plot in order
growth_summary$Host <- factor(growth_summary$Host, levels = c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA"))  
dim(growth_summary) # 79  8



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # statistics on K

# for Levene, H0 = variances are equal # for Shapiro, H0 = data are normally distributed
leveneTest(Quantity ~ Pathogen, data = k_data)     #  0.9878 var are equal
shapiro.test(k_data$Quantity)                      #  0.1751 data are normally distributed

# 2 way ANCOVA 
mod <- lm(Quantity ~ Pathogen * Host + Plate , data = k_data)
anova(mod)
# Response: Quantity
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Pathogen       4  5.8650 1.46626  3.9217 0.005895 **
# Host           4  1.7019 0.42549  1.1380 0.344789   
# Plate          1  1.2950 1.29503  3.4637 0.066450 . 
# Pathogen:Host 16  6.9771 0.43607  1.1663 0.313206  
# Residuals     79 29.5370 0.37389                    


# # # planned contrasts for native host vs all other hosts (K local adaptation)
a9 = subset(k_data, subset = Pathogen == "A9")
dc = subset(k_data, subset = Pathogen == "DC3000")
es = subset(k_data, subset = Pathogen == "ES4326")
a14 = subset(k_data, subset = Pathogen == "1448A")
np = subset(k_data, subset = Pathogen == "NP29")

# A9 
a9.aov <- aov(Quantity ~ Plate + Host, data = a9, contrasts = list(Host = "contr.sum"))
a9.emml <- emmeans(a9.aov, consec ~ Host)
a9$Host=as.factor(a9$Host)
levels(a9$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modela9 <- aov(Quantity~Plate+Host, data = a9)
a9_mod <- emmeans(modela9, specs = ~ Host)
contrast(a9_mod, method = list(c(-4,1,1,1,1))) # 0.1114  

# DC 
dc.aov <- aov(Quantity ~ Plate + Host, data = dc, contrasts = list(Host = "contr.sum"))
dc.emml <- emmeans(dc.aov, consec ~ Host)
dc$Host=as.factor(dc$Host)
levels(dc$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modeldc <- aov(Quantity~Plate+Host, data = dc)
dc_mod <- emmeans(modeldc, specs = ~ Host)
contrast(dc_mod, method = list(c(1,-4,1,1,1))) # 0.7328

# ES 
es.aov <- aov(Quantity ~ Plate + Host, data = es, contrasts = list(Host = "contr.sum"))
es.emml <- emmeans(es.aov, consec ~ Host)
es$Host=as.factor(es$Host)
levels(es$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modeles <- aov(Quantity~Plate+Host, data = es)
es_mod <- emmeans(modeles, specs = ~ Host)
contrast(es_mod, method = list(c(1,1,-4,1,1))) # 0.6580

# 1448A 
a14.aov <- aov(Quantity ~ Plate + Host, data = a14, contrasts = list(Host = "contr.sum"))
a14.emml <- emmeans(a14.aov, consec ~ Host)
a14$Host=as.factor(a14$Host)
levels(a14$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modela14 <- aov(Quantity~Plate+Host, data = a14)
a14_mod <- emmeans(modela14, specs = ~ Host)
contrast(a14_mod, method = list(c(1,1,1,-4,1))) # 0.0447 # SIG

# NP 
np.aov <- aov(Quantity ~ Plate + Host, data = np, contrasts = list(Host = "contr.sum"))
np.emml <- emmeans(np.aov, consec ~ Host)
np$Host=as.factor(np$Host)
levels(np$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modelnp <- aov(Quantity~Plate+Host, data = np)
np_mod <- emmeans(modelnp, specs = ~ Host)
contrast(np_mod, method = list(c(1,1,1,1,-4))) # 0.0425 # SIG



# # # planned contrasts for focal pathogen vs. all other pathogens on each host (K - local dominance)
growth_summary = subset(growth_summary, subset = Host != "MOCK")
tom_a9 = subset(k_data, subset = Host == "TOMATO-A9")
tom_dc = subset(k_data, subset = Host == "TOMATO-DC")
bean = subset(k_data, subset = Host == "BEAN")
radish = subset(k_data, subset = Host == "RADISH")
at = subset(k_data, subset = Host == "A.THALIANA")

# TOM-A9 
tom_a9.aov <- aov(Quantity ~ Plate + Pathogen, data = tom_a9, contrasts = list(Host = "contr.sum"))
tom_a9.emml <- emmeans(tom_a9.aov, consec ~ Pathogen)
tom_a9$Pathogen=as.factor(tom_a9$Pathogen)
levels(tom_a9$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modeltom_a9 <- aov(Quantity~Plate+Pathogen, data = tom_a9)
tom_a9_mod <- emmeans(modeltom_a9, specs = ~ Pathogen)
contrast(tom_a9_mod, method = list(c(-4,1,1,1,1))) # 0.0368 **** A9 clearly does the best

# TOM-DC
tom_dc.aov <- aov(Quantity ~ Plate + Pathogen, data = tom_dc, contrasts = list(Host = "contr.sum"))
tom_dc.emml <- emmeans(tom_dc.aov, consec ~ Pathogen)
tom_dc$Pathogen=as.factor(tom_dc$Pathogen)
levels(tom_dc$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modeltom_dc <- aov(Quantity~Plate+Pathogen, data = tom_dc)
tom_dc_mod <- emmeans(modeltom_dc, specs = ~ Pathogen)
contrast(tom_dc_mod, method = list(c(1,-4,1,1,1))) # 0.0023  **** DC does best 

# RADISH
radish.aov <- aov(Quantity ~ Plate + Pathogen, data = radish, contrasts = list(Host = "contr.sum"))
radish.emml <- emmeans(radish.aov, consec ~ Pathogen)
radish$Pathogen=as.factor(radish$Pathogen)
levels(radish$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelradish <- aov(Quantity~Plate+Pathogen, data = radish)
radish_mod <- emmeans(modelradish, specs = ~ Pathogen)
contrast(radish_mod, method = list(c(1,1,-4,1,1))) # 0.9343 ES not better
contrast(radish_mod, method = list(c(1,-4,1,1,1))) # 0.0674 DC not better
contrast(radish_mod, method = list(c(-4,1,1,1,1))) # 0.9568 A9 not better

# BEAN
bean.aov <- aov(Quantity ~ Plate + Pathogen, data = bean, contrasts = list(Host = "contr.sum"))
bean.emml <- emmeans(bean.aov, consec ~ Pathogen)
bean$Pathogen=as.factor(bean$Pathogen)
levels(bean$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelbean <- aov(Quantity~Plate+Pathogen, data = bean)
bean_mod <- emmeans(modelbean, specs = ~ Pathogen)
contrast(bean_mod, method = list(c(1,1,1,-4,1))) # 0.4858 14 does not win
contrast(bean_mod, method = list(c(1,-4,1,1,1))) # 0.2273 DC not better

# AT
at.aov <- aov(Quantity ~ Plate + Pathogen, data = at, contrasts = list(Host = "contr.sum"))
at.emml <- emmeans(at.aov, consec ~ Pathogen)
at$Pathogen=as.factor(at$Pathogen)
levels(at$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelat <- aov(Quantity~Plate+Pathogen, data = at)
at_mod <- emmeans(modelat, specs = ~ Pathogen)
contrast(at_mod, method = list(c(1,1,1,1,-4))) # 0.3152 NP does not win
contrast(at_mod, method = list(c(1,1,-4,1,1))) # 0.4202 ES does not win
contrast(at_mod, method = list(c(1,-4,1,1,1))) # 0.6758 DC does not win


# # # # # # # pathogen vs. all other pathogens across hosts (K - overall performance)
k_data = subset(k_data, subset = Pathogen != "MOCK")
k_data <- droplevels(subset(k_data, Pathogen != "MOCK"))

model <- lm(Quantity ~ Pathogen + Host, data = k_data)
emm_path <- emmeans(model, ~ Pathogen)
contrasts_vs_all <- contrast(emm_path, method = "eff")  # "eff" = effect contrasts
summary(contrasts_vs_all, infer = TRUE)
# contrast      estimate    SE df lower.CL upper.CL t.ratio p.value
# A9 effect       0.1521 0.131 96  -0.1915   0.4958   1.163  0.3094
# DC3000 effect   0.3899 0.128 96   0.0529   0.7270   3.040  0.0152***
# ES4326 effect  -0.0715 0.123 96  -0.3951   0.2521  -0.581  0.5628
# 1448A effect   -0.2970 0.120 96  -0.6127   0.0187  -2.472  0.0379***
# NP29 effect    -0.1735 0.118 96  -0.4846   0.1375  -1.466  0.2431
# Results are averaged over the levels of: Host 
# Confidence level used: 0.95 
# Conf-level adjustment: bonferroni method for 5 estimates 
# P value adjustment: fdr method for 5 tests (False Discovery Rate)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # statistics on r

# # # 2 way ANCOVA
# for Levene, H0 = variances are equal # for Shapiro, H0 = data are normally distributed
growth_summary$logslope = log10(growth_summary$slope)
leveneTest(logslope ~ Pathogen, data = growth_summary)     #  var are equal
shapiro.test(growth_summary$logslope)                      #  data are normally distributed
mod <- lm(logslope ~ Pathogen * Host + Plate , data = growth_summary)
anova(mod)
# Response: slope
#               Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Pathogen       4 0.8439 0.21097 18.6455 1.280e-09 ***
# Host           4 1.2778 0.31944 28.2318 1.397e-12 ***
# Plate          1 0.0303 0.03035  2.6823    0.1074    
# Pathogen:Host 16 6.1467 0.38417 33.9525 < 2.2e-16 ***
#  Residuals     53 0.5997 0.01131 


# # # # # # # planned contrasts for native host vs all other hosts (r local adaptation)
a9 = subset(growth_summary, subset = Pathogen == "A9")
dc = subset(growth_summary, subset = Pathogen == "DC3000")
es = subset(growth_summary, subset = Pathogen == "ES4326")
a14 = subset(growth_summary, subset = Pathogen == "1448A")
np = subset(growth_summary, subset = Pathogen == "NP29")
mock = subset(growth_summary, subset = Pathogen == "MOCK")

# A9 
a9.aov <- aov(slope ~ Plate + Host, data = a9, contrasts = list(Host = "contr.sum"))
a9.emml <- emmeans(a9.aov, consec ~ Host)
a9$Host=as.factor(a9$Host)
levels(a9$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modela9 <- aov(slope~Plate+Host, data = a9)
a9_mod <- emmeans(modela9, specs = ~ Host)
contrast(a9_mod, method = list(c(-4,1,1,1,1))) # 0.0002****

# # # A9 & DC are so similar
contrast(a9_mod, method = list(c(1,-4,1,1,1))) # 0.0005 this constrast is also sig for A9 in Tom-DC
# is A9 r in the tomatoes equivalent?
ta9 = subset(a9, subset = Host =="TOMATO-A9")
tdc = subset(a9, subset = Host =="TOMATO-DC")
t.test(ta9$slope, tdc$slope) # p-value = 0.5637, true difference in means is not equal to 0


# DC 
dc.aov <- aov(slope ~ Plate + Host, data = dc, contrasts = list(Host = "contr.sum"))
dc.emml <- emmeans(dc.aov, consec ~ Host)
dc$Host=as.factor(dc$Host)
levels(dc$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modeldc <- aov(slope~Plate+Host, data = dc)
dc_mod <- emmeans(modeldc, specs = ~ Host)
contrast(dc_mod, method = list(c(1,-4,1,1,1))) # 0.8225

# ES 
es.aov <- aov(slope ~ Plate + Host, data = es, contrasts = list(Host = "contr.sum"))
es.emml <- emmeans(es.aov, consec ~ Host)
es$Host=as.factor(es$Host)
levels(es$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modeles <- aov(slope~Plate+Host, data = es)
es_mod <- emmeans(modeles, specs = ~ Host)
contrast(es_mod, method = list(c(1,1,-4,1,1))) # 0.0415*******
contrast(es_mod, method = list(c(1,1,1,-4,1))) # 0.0074 on bean*******
contrast(es_mod, method = list(c(1,1,1,1,-4))) # 0.0003 on At*******

# 1448A 
a14.aov <- aov(slope ~ Plate + Host, data = a14, contrasts = list(Host = "contr.sum"))
a14.emml <- emmeans(a14.aov, consec ~ Host)
a14$Host=as.factor(a14$Host)
levels(a14$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modela14 <- aov(slope~Plate+Host, data = a14)
a14_mod <- emmeans(modela14, specs = ~ Host)
contrast(a14_mod, method = list(c(1,1,1,-4,1))) # 0.3651

# NP 
np.aov <- aov(slope ~ Plate + Host, data = np, contrasts = list(Host = "contr.sum"))
summary(np.aov)
np.emml <- emmeans(np.aov, consec ~ Host)
np$Host=as.factor(np$Host)
levels(np$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modelnp <- aov(slope~Plate+Host, data = np)
np_mod <- emmeans(modelnp, specs = ~ Host)
contrast(np_mod, method = list(c(1,1,1,1,-4))) # <.0001*****



# # # # # # # planned contrasts for focal pathogen vs. all other pathogens on each native host (r - local dominance)
growth_summary = subset(growth_summary, subset = Host != "MOCK")
tom_a9 = subset(growth_summary, subset = Host == "TOMATO-A9")
tom_dc = subset(growth_summary, subset = Host == "TOMATO-DC")
bean = subset(growth_summary, subset = Host == "BEAN")
radish = subset(growth_summary, subset = Host == "RADISH")
at = subset(growth_summary, subset = Host == "A.THALIANA")

# TOM-A9 
tom_a9.aov <- aov(slope ~ Plate + Pathogen, data = tom_a9, contrasts = list(Host = "contr.sum"))
tom_a9.emml <- emmeans(tom_a9.aov, consec ~ Pathogen)
tom_a9$Pathogen=as.factor(tom_a9$Pathogen)
levels(tom_a9$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modeltom_a9 <- aov(slope~Plate+Pathogen, data = tom_a9)
tom_a9_mod <- emmeans(modeltom_a9, specs = ~ Pathogen)
contrast(tom_a9_mod, method = list(c(-4,1,1,1,1))) # <.0001 **** A9 clearly does the best

# TOM-DC
tom_dc.aov <- aov(slope ~ Plate + Pathogen, data = tom_dc, contrasts = list(Host = "contr.sum"))
tom_dc.emml <- emmeans(tom_dc.aov, consec ~ Pathogen)
tom_dc$Pathogen=as.factor(tom_dc$Pathogen)
levels(tom_dc$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modeltom_dc <- aov(slope~Plate+Pathogen, data = tom_dc)
tom_dc_mod <- emmeans(modeltom_dc, specs = ~ Pathogen)
contrast(tom_dc_mod, method = list(c(1,-4,1,1,1))) # 0.8465 DC does not do the best
contrast(tom_dc_mod, method = list(c(-4,1,1,1,1))) # <.0001 A9 clearly does the best

# RADISH
radish.aov <- aov(slope ~ Plate + Pathogen, data = radish, contrasts = list(Host = "contr.sum"))
radish.emml <- emmeans(radish.aov, consec ~ Pathogen)
radish$Pathogen=as.factor(radish$Pathogen)
levels(radish$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelradish <- aov(slope~Plate+Pathogen, data = radish)
radish_mod <- emmeans(modelradish, specs = ~ Pathogen)
contrast(radish_mod, method = list(c(1,1,-4,1,1))) # 0.0087 ES these 3 all do better than the rest/equally well
contrast(radish_mod, method = list(c(1,-4,1,1,1))) # 0.0002 DC
contrast(radish_mod, method = list(c(1,1,1,-4,1))) # 0.0201 14

# BEAN
bean.aov <- aov(slope ~ Plate + Pathogen, data = bean, contrasts = list(Host = "contr.sum"))
bean.emml <- emmeans(bean.aov, consec ~ Pathogen)
bean$Pathogen=as.factor(bean$Pathogen)
levels(bean$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelbean <- aov(slope~Plate+Pathogen, data = bean)
bean_mod <- emmeans(modelbean, specs = ~ Pathogen)
contrast(bean_mod, method = list(c(1,1,1,-4,1))) # 0.3129 14 does not win
contrast(bean_mod, method = list(c(1,1,-4,1,1))) # <.0001 ES clearly does best

# AT
at.aov <- aov(slope ~ Plate + Pathogen, data = at, contrasts = list(Host = "contr.sum"))
at.emml <- emmeans(at.aov, consec ~ Pathogen)
at$Pathogen=as.factor(at$Pathogen)
levels(at$Pathogen) # check order of levels "A9"     "DC3000" "ES4326" "1448A"  "NP29"   "MOCK"  
modelat <- aov(slope~Plate+Pathogen, data = at)
at_mod <- emmeans(modelat, specs = ~ Pathogen)
contrast(at_mod, method = list(c(1,1,1,1,-4))) # 0.0013 NP these 3 all do equally well/better than the rest
contrast(at_mod, method = list(c(1,1,-4,1,1))) # 0.0025 ES ES & NP not sig different (t-test)
contrast(at_mod, method = list(c(-4,1,1,1,1))) # 0.0099 A9 



# # # # # # # contrasts for each pathogen vs all others (r - overall growth)
model <- lm(slope ~ Pathogen + Host, data = growth_summary)
emm_path <- emmeans(model, ~ Pathogen)
contrasts_vs_all <- contrast(emm_path, method = "eff")  # "eff" = effect contrasts
summary(contrasts_vs_all, infer = TRUE)

# contrast      estimate     SE df  lower.CL upper.CL t.ratio p.value
# A9 effect       0.0293 0.0108 70  0.000761  0.05789   2.718  0.0413 ***
# DC3000 effect  -0.0105 0.0108 70 -0.039029  0.01810  -0.970  0.3355
# ES4326 effect   0.0174 0.0108 70 -0.011150  0.04598   1.614  0.1631
# 1448A effect   -0.0211 0.0108 70 -0.049697  0.00743  -1.959  0.1353
# NP29 effect    -0.0151 0.0099 70 -0.041363  0.01107  -1.530  0.1631


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # plot Fig. 2
print(k_data, n = 90)

k_data$Time = NULL
k_data$Infection = NULL
names(k_data) = c("Host", "Pathogen", "Plate", "K")
merged_data <- merge(growth_summary, k_data[, c("Host", "Pathogen", "Plate", "K")],
                     by = c("Host", "Pathogen", "Plate"),all.x = TRUE)
merged_data <- subset(merged_data, !is.na(K))

host_levels <- c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA")
pathogen_levels <- c("A9", "DC3000", "ES4326", "1448A", "NP29", "MOCK")

stats <- merged_data %>% group_by(Host, Pathogen) %>% 
  summarise(N = length(auc), K_mean = mean(K),Ksd = sd(K),R_mean = mean(slope), Rsd = sd(slope),  
            AUC_mean = mean(auc),AUCsd = sd(auc), .groups = "drop_last")

# Hosts to shade per pathogen 
shade_info <- data.frame(Pathogen = c("A9", "DC3000", "ES4326", "1448A", "NP29"),
                         Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"))

# Convert Host names to numeric x-axis positions
shade_info$xpos <- match(shade_info$Host, host_levels)
rects <- shade_info |> transform(xmin = xpos - 0.5,xmax = xpos + 0.5,ymin = -Inf, ymax = Inf)
#rects$Pathogen <- factor(rects$Pathogen, levels = levels(stats$Pathogen))

stats$Host <- factor(stats$Host, levels = host_levels)
stats$Pathogen <- factor(stats$Pathogen, levels = pathogen_levels)
rects$Host <- factor(rects$Host, levels = host_levels)
rects$Pathogen <- factor(rects$Pathogen, levels = pathogen_levels)

# add astric labels for the significant contrasts
signif_labels_r <- data.frame(Host = c("TOMATO-A9", "A.THALIANA"),Pathogen = c("A9", "NP29"), y_pos = 0.18, label = "*")
signif_labels_r$Host <- factor(signif_labels_r$Host, levels = host_levels)
signif_labels_r$Pathogen <- factor(signif_labels_r$Pathogen, levels = pathogen_levels)

signif_labels_k <- data.frame(Host = c("BEAN", "A.THALIANA"),Pathogen = c("1448A", "NP29"), y_pos = 8.1, label = "*")
signif_labels_k$Host <- factor(signif_labels_k$Host, levels = host_levels)
signif_labels_k$Pathogen <- factor(signif_labels_k$Pathogen, levels = pathogen_levels)




plot_r <- ggplot(stats, aes(x = Host, y = R_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2.2) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(R_mean - (0.5*Rsd)), ymax=(R_mean + (0.5*Rsd))))+
  geom_jitter(data = merged_data, aes(x = Host, y = slope, color = Pathogen),
              size = 1.5, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1) +
  geom_text(data = signif_labels_r, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(y = "Growth Rate", x = "",  color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 11)) +
  theme(strip.text = element_text(size = 11))


plot_k <- ggplot(stats, aes(x = Host, y = K_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2.2) +
  geom_errorbar(width=0.9, size=0.7, aes(ymin=(K_mean - (0.5*Ksd)), ymax=(K_mean + (0.5*Ksd))))+
  geom_jitter(data = k_data, aes(x = Host, y = K, color = Pathogen), 
              size = 1.5, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  geom_text(data = signif_labels_k, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  
  facet_wrap(~ Pathogen, nrow = 1)+#, strip.position = "right") +
  labs(y = "Log10(K)", color = "Pathogen") +
  coord_cartesian(ylim = c(5, 8.5)) +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), #axis.title.x = element_text(size = 11), 
        axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_blank())

tight_margin <- theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))  # top, right, bottom, left
plot_k <- plot_k + tight_margin
plot_r <- plot_r + tight_margin

combined_plot <- (plot_r / plot_k) + 
  plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = 'a')

ggsave(
  filename = file.path(output_dir, "Fig3-rev.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 17.4,
  height = 10,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "Fig3-rev.eps"),
  plot = combined_plot,
  device = cairo_ps,      # EPS device with embedded fonts
  width = 17.4,
  height = 10,
  units = "cm"
)
#########

# # # # # FIG 4 - pangenome pdf, anvio
# # # # # FIG 5 - tree pdf, itol

# # # # # FIG 6 - EFFECTORS VENN + TABLE width = 17.4, height = 10
#########
#DC     36 families, 54 total alleles, 47 unique alleles - Martel 2022
#ES     26 families, 34 total alleles, 33 unique alleles - Martel 2022
#1448A  22 families, 32 total alleles (PsytTec), 27 functionally verified (Zumaquero 2010)

library(VennDiagram)
library(ggVennDiagram)
library(gt)
library(tidyverse)
library(grid)
library(gridExtra)
library(functional)
library(gggenes)

data <- read.csv(file.path(folder_fig6, "data_Fig6.csv"), header = TRUE)
data = data[order(data$Allele),]

summary_table <- data %>% summarise(across(c(DC3000, ES4326, X1448A, A9, NP29), sum))
#DC3000 ES4326 A1448 A9 NP29 # alleles
#    53     34    32 62   12

# set family names
data <- data %>% mutate(family = case_when(
  grepl("^Hop", Allele) ~ paste0("Hop", substr(Allele, 4, 5)),
  grepl("^Avr", Allele) ~ paste0("Avr", substr(Allele, 4, 4)),
  TRUE ~ NA_character_),
  family = gsub("[0-9]+$", "", family))

data <- data %>% dplyr::mutate(family = dplyr::recode(family, "AvrP"="AvrPto","AvrR"="AvrRpm"))
strain_cols <- c("A9", "DC3000", "ES4326", "X1448A", "NP29")

family_sets <- lapply(strain_cols, function(col) {
  data %>% filter(.data[[col]] != 0) %>%   # filter rows where strain count is non-zero
    pull(family) %>% unique()}) # manually confirm 36 for DC, 26 for ES, 22 for 1448a

names(family_sets) <- c("A9", "DC3000", "ES4326", "1448A","NP29")
count_only_sets <- lapply(family_sets, function(x) rep("", length(x)))
my_colors <- c("A9" = "sienna1", "DC3000" = "skyblue4", "ES4326" = "violet","1448A" = "mediumaquamarine","NP29" = "olivedrab3")

venn.plot <- venn.diagram(x = family_sets, filename = NULL, fill = my_colors, col  = my_colors,     # <- outlines match fills
  alpha = 0.4, lwd = 4,
  fontfamily = "Helvetica", cex = 1.2,  # inner labels (numbers)
  cat.cex = 1.2, cat.fontfamily = "Helvetica", cat.col = my_colors, cat.fontface = "bold",# outer labels (strains)
  cat.pos = c(-48, -7, -115, 136, 20), cat.dist = c(0.21, 0.22, 0.28, 0.23, 0.22)) 

venn_grob <- grobTree(venn.plot)

# Get unique families per strain
unique_families_per_strain <- lapply(names(family_sets), function(strain) {
  other_strains <- setdiff(names(family_sets), strain)       # all other strain names
  other_families <- unlist(family_sets[other_strains])      # combine all other families
  setdiff(family_sets[[strain]], other_families)            # families only in this strain
})


# convert list to a data frame for tableGrob
names(unique_families_per_strain) <- c("A9", "DC3000", "ES4326", "1448A","NP29")
max_len <- max(sapply(unique_families_per_strain, length))

# Convert to a data frame with padded NAs for shorter columns
table_df <- as.data.frame(lapply(unique_families_per_strain, function(x) {
  length(x) <- max_len  # pad with NAs
  x
}), stringsAsFactors = FALSE)

# Reorder columns alphabetically
table_df <- table_df[, order(colnames(table_df))]

# Replace NAs with blanks
table_df[is.na(table_df)] <- ""

# Rename the column "1448A" 
colnames(table_df) <- gsub("^X", "", colnames(table_df))

# Drop NP29 which has no unique effectors
table_df$NP29 = NULL

# Define theme with column-specific background colors
alpha = 0.5
my_cols_trans <- adjustcolor(my_colors, alpha.f = alpha)

my_theme <- ttheme_default(core = list(fg_params = list(cex = 1),        
  bg_params = list(fill = rep(my_cols_trans, each = nrow(table_df))), padding = unit(c(7, 14), "points")),
  colhead = list(fg_params = list(cex = 1.1, fontface = "bold"), bg_params = list(fill = my_cols_trans),
                 padding = unit(c(7, 20), "points")
  ))

table_grob <- tableGrob(table_df, rows = NULL, theme = my_theme)

# make combined plot with a/b labels
venn_plot <- wrap_elements(venn_grob)
table_plot <- wrap_elements(table_grob)

combined_plot <- (venn_plot + table_plot) + plot_layout(widths = c(1.8, 1.2)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  filename = file.path(output_dir, "Fig6-R-rev.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 17.4,
  height = 10,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "Fig6-R-rev.eps"),
  plot = combined_plot,
  device = cairo_ps,   # EPS device with embedded fonts
  width = 17.4,
  height = 10,
  units = "cm"
)#########

# # # # # FIG 7 - VFs presence/absence width = 17.4, height = 10
#########

a9  <- read.delim(file.path(folder_fig7, "ABRicate_A9.tab"))
dc  <- read.delim(file.path(folder_fig7, "ABRicate_DC.tab"))
es  <- read.delim(file.path(folder_fig7, "ABRicate_ES.tab"))
a14 <- read.delim(file.path(folder_fig7, "ABRicate_14.tab"))
np  <- read.delim(file.path(folder_fig7, "ABRicate_NP.tab"))
data = rbind(a9,dc,es,a14,np)
data <- data %>% mutate(X.FILE = dplyr::recode(X.FILE, "A9.fasta"="A9","DC3000.fna"="DC3000","ES4326.fna"="ES4326","A1448A.fna"="1448A","NP29.fasta"="NP29"))

# group by strain, type
table <- data %>% group_by(X.FILE,GENE) %>% summarise(count = n(), .groups = "drop")
table = table[order(table$count),]
tail(table)
dim(table) # 89  3

# add categories to table for plotting
biof = table[grep("alg", table$GENE),]
biof$Category = "biofilm"

sec1 = table[grep("clp", table$GENE),]
sec2 = table[grep("hcp", table$GENE),]
sec3 = table[grep("hsi", table$GENE),]
sec = rbind(sec1,sec2,sec3)
sec$Category = "T6SS"

fla1 = table[grep("fle", table$GENE),]
fla2 = table[grep("flg", table$GENE),]
fla3 = table[grep("fli", table$GENE),]
fla = rbind(fla1,fla2,fla3)
fla$Category = "flagella"

pyo1 = table[grep("mbt", table$GENE),]
pyo2 = table[grep("pvd", table$GENE),]
pyo = rbind(pyo1,pyo2)
pyo$Category = "pyover"

table = rbind(biof,fla,pyo)
names(table) = c("Strain", "Gene", "Count", "Category")
dim(table) #  68  4

# manually add additional toxins
# all have mangotoxin mgoB but not the rest of the operon
# all have ttr - not tabtoxin but 1/2 self-protective genes

tail(table)
table[69,1] = "NP29"
table[69,2] = "syrB1"
table[69,3] = 1
table[69,4] = "s-myc"

table[70,1] = "NP29"
table[70,2] = "syrB2"
table[70,3] = 1
table[70,4] = "s-myc"

table[71,1] = "NP29"
table[71,2] = "sypA"
table[71,3] = 1
table[71,4] = "s-pep" 

table[72,1] = "NP29"
table[72,2] = "sypB"
table[72,3] = 1
table[72,4] = "s-pep" 

table[73,1] = "1448A"
table[73,2] = "amtA"
table[73,3] = 1
table[73,4] = "phaseo" 

table[74,1] = "1448A"
table[74,2] = "argK"
table[74,3] = 1
table[74,4] = "phaseo" # also has 17 gene Pht cluster

table[75,1] = "A9"
table[75,2] = "cfa"
table[75,3] = 1
table[75,4] = "corona" 

table[76,1] = "1448A"
table[76,2] = "cfa"
table[76,3] = 1
table[76,4] = "corona" 

table[77,1] = "DC3000"
table[77,2] = "cfa"
table[77,3] = 1
table[77,4] = "corona" 

table[78,1] = "DC3000"
table[78,2] = "cfa2"
table[78,3] = 1
table[78,4] = "corona" 

table[79,1] = "ES4326"
table[79,2] = "cfa"
table[79,3] = 1
table[79,4] = "corona" 

table[80,1] = "ES4326"
table[80,2] = "cfa2"
table[80,3] = 1
table[80,4] = "corona" 

table[81,1] = "ES4326"
table[81,2] = "cfaB"
table[81,3] = 1
table[81,4] = "corona" 

table[82,1] = "NP29"
table[82,2] = "cfa"
table[82,3] = 1
table[82,4] = "corona" 

table[83,1] = "NP29"
table[83,2] = "cfaB"
table[83,3] = 1
table[83,4] = "corona" 

table[84,1] = "DC3000"
table[84,2] = "Cma cluster"
table[84,3] = 1
table[84,4] = "corona" 

table[85,1] = "ES4326"
table[85,2] = "Cma cluster"
table[85,3] = 1
table[85,4] = "corona" 

table[86,1] = "NP29"
table[86,2] = "sypC"
table[86,3] = 1
table[86,4] = "s-pep" 

table[87,1] = "NP29"
table[87,2] = "syrD"
table[87,3] = 1
table[87,4] = "s-myc" 

table[88,1] = "1448A"
table[88,2] = "Pht cluster"
table[88,3] = 1
table[88,4] = "phaseo" # also has 17 gene Pht cluster

table[89,1] = "1448A"
table[89,2] = "cfaB"
table[89,3] = 1
table[89,4] = "corona" 

strain_order <- c("A9","DC3000","ES4326","1448A","NP29")
table$Strain <- factor(table$Strain, levels = strain_order)
table$Count <- as.factor(table$Count)
table$Gene <- factor(table$Gene, levels = unique(table$Gene))

plot = ggplot(table, aes(x = Gene, y = Strain, fill = Strain)) + geom_tile(color = "white") +
  labs(x = "Gene Name", y = "Pathogen", fill = "Count") +
  scale_fill_manual(values = c("A9" = "sienna1","DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  facet_grid(~ Category, scales = "free_x", space = "free_x") +
  theme_minimal() + theme(legend.position = "none",
                          panel.grid = element_blank(),panel.background = element_rect(fill = "white", color = NA),
                          plot.background = element_rect(fill = "white", color = NA),
                          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
                          axis.text.y = element_text(size = 11),
                          axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
                          strip.background = element_blank(),
                          strip.text = element_text(size = 12))#, face = "bold"))

ggsave(
  filename = file.path(output_dir, "Fig7-R-rev.pdf"),
  plot = plot,
  device = "pdf",
  width = 17.4,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "Fig7-R-rev.eps"),
  plot = plot,
  device = cairo_ps,   # EPS device with embedded fonts
  width = 17.4,
  height = 7,
  units = "cm"
)
#########

# # # # # FIG S1 - OD/CFU width = 8.4, height = 6.4
#########
library(dplyr) 
library(ggplot2)

data <- read.table(file.path(folder_figS1, "data_FigS1.txt"), header = TRUE)
dim(data) # 91  5
head(data)

# cfu/ml = counts/(1/df)*vol plated
data[6] = (data$count/(1/data$df))*data$ml_plated
data[6] = log10(data[6])
head(data)

namevec = c("strain","od","count","df","ml_plated","cfuml")
names(data) =namevec
#head(data)

# remove 2 ES obs = 0 (drop as plating errors)
dim(data) # 75  6
data=subset(data, subset = cfuml > 1)
dim(data) # 73  6
data=subset(data, subset = count <2000)
dim(data) # 71  6

# function to calculate 95% CIs for a vector using a t-distribution
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)}

# calculate stats! median and CIs for all replicates combined
stats <- data %>% group_by(strain, od) %>%
  summarise(N=length(cfuml), Mean=mean(cfuml), CI95=conf_int95(cfuml), StDev=sd(cfuml))
dim(stats) # 15  6
head(stats)
stats$od=as.factor(stats$od)

summary(glm(stats$Mean~stats$od))
mo = summary(lm(stats$Mean~stats$od)) # 0.9613 

stats <- stats %>% dplyr::mutate(strain = dplyr::recode(strain, "A1448"="1448A"))
strain_order <- c("A9","DC3000","ES4326","1448A","NP29")
stats$strain <- factor(stats$strain, levels = strain_order)

dodge =position_dodge(.3) # how much jitter on the x-axis?
plot = ggplot(data=stats, aes(x=od, y=Mean, color=strain)) + 
  geom_point(size=3, position = dodge) + 
  geom_errorbar(aes(x=od, width = 1, ymin=Mean-(1/2*CI95), ymax=Mean+(1/2*CI95)),position = dodge) +
  labs(y = "log10(CFU/mL)", color="Pathogen") +
  scale_x_discrete(name ="OD600", limits=c("0.01","0.1","1")) +
  scale_color_manual(values = c("A9" = "sienna1","DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray80", size = 0.5), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),axis.title.y = element_text(size = 12))

test = (aov(cfuml~strain, data)) # ANOVA
summary(test)
#             Df Sum Sq  Mean Sq  F value  Pr(>F)
# strain       4   3.85  0.9623   0.568    0.687
# Residuals   82 138.92  1.6941  

ggsave(
  filename = file.path(output_dir, "FigS1-R-rev.pdf"),
  plot = plot,
  device = "pdf",
  width = 8.4,
  height = 6.4,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "FigS1-R-rev.eps"),
  plot = plot,
  device = cairo_ps,   # EPS device with embedded fonts
  width = 8.4,
  height = 6.4,
  units = "cm"
)


#########


