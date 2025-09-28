# Rebecca Satterwhite
# 5/30/25  This is the main planta script for matrix fig + stats
# k = highest titer achieved, mean and SE

# PIPELINE
# 1) import data
# 2) drop libs, drop controls
# 3) convert CT to Q
# 4) account for dilution factor 
# 5) change Q=NA to Q=0
# 6) find k (highest Q)
# 7) stats to get mean and error around k
# 8) ANOVA on k data frame (not stats data frame)
# 9) find r with linear regression


# # # # # # # set up 
########
library(ggplot2)
library(tidyr)
library(patchwork)
library(emmeans)
library(car)
library(dplyr)  
library(broom)
library(MESS)
library(purrr)

# # # import raw data & put into 1 dataframe
data1 = read.table(".../2_1_6.txt", header=TRUE) 
data2 = read.table(".../2_2_5a.txt", header=TRUE) 
data3 = read.table(".../2_4_2.txt", header=TRUE)  
data4 = read.table(".../2_4_3.txt", header=TRUE) 
data5 = read.table(".../2_11_6.txt", header=TRUE)  
data6 = read.table(".../2_13_4.txt", header=TRUE)  
data7 = read.table(".../2_14_6.txt", header=TRUE)  
data8 = read.table(".../2_15_2.txt", header=TRUE)  
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

# # # # check that mock infections had sig lower pop sizes

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

p_val <- test_result$p.value
p_label <- paste0("p = ", signif(p_val, 3))  # rounds to 3 significant digits

mock_check <- ggplot(combined_df, aes(x = Infection, y = Quantity)) +
  geom_boxplot() +
  #labs(title = "Comparison of Bacterial Quantities: Mock vs Real Infections") +
  theme_minimal() +
  annotate("text", x = 1, y = max(combined_df$Quantity, na.rm = TRUE)*0.95, label = p_label, size = 4)

ggsave(".../f_mock_comparison.pdf", mock_check, width = 3, height = 3)


########


# # # # # # # initial plots
########
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

mat = ggplot(mat_dat, aes(x = Time, y = Quantity, fill = Pathogen, color = Pathogen)) +
  geom_point(size=1.1)+
  geom_errorbar(width=0.8, size = 0.7, aes(ymin=(Quantity - Qsd), ymax=(Quantity + Qsd)))+
  facet_grid(Host~Pathogen) +
  labs(y = "Log10(Quantity)", x = "Host", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3", "MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+#theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  
#strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))
ggsave(".../f_matrix.pdf", mat, width = 5, height = 5)




# # # # # # # plots of k 

k_data= subset(k_data, subset = Pathogen != "MOCK")

plot_kp <- ggplot(k_data, aes(x = Host, y = Quantity, fill = Pathogen, color = Pathogen)) +
  geom_boxplot(alpha=0.3, size=0.8, outlier.size=1, outlier.alpha=1, width=0.7)+#, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Host, nrow = 1) +
  labs(y = "Log10(Quantity)", x = "Host", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3", "MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+#theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  
#strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))

plot_kh <- ggplot(k_data, aes(x = Pathogen, y = Quantity, fill = Pathogen, color = Pathogen)) +
  geom_boxplot(alpha=0.3, size = 0.8, outlier.size=1, outlier.alpha=1, width = 0.7)+ #, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Host, nrow = 1)+#, strip.position = "right") +
  labs(y = "Log10(Quantity)", x = "Pathogen", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ #theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  # smaller text & less margin 
#strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))

combined_plot <- (plot_kp / plot_kh) + plot_annotation(tag_levels = 'A')
# stacked vertically; use | for side-by-side
ggsave("/.../f_k-box-by-host.pdf", combined_plot, width = 7, height = 6)




# # # the one i'm using (by pathogen)

k_data= subset(k_data, subset = Pathogen != "MOCK")

plot_kp <- ggplot(k_data, aes(x = Host, y = Quantity, fill = Pathogen, color = Pathogen)) +
  geom_boxplot(alpha=0.3, size=0.8, outlier.size=1, outlier.alpha=1, width=0.7)+#, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Pathogen, nrow = 1) +
  labs(y = "Log10(Quantity)", x = "Host", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3", "MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+#theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  
      #strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))

plot_kh <- ggplot(k_data, aes(x = Pathogen, y = Quantity, fill = Pathogen, color = Pathogen)) +
  geom_boxplot(alpha=0.3, size = 0.8, outlier.size=1, outlier.alpha=1, width = 0.7)+ #, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Host, nrow = 1)+#, strip.position = "right") +
  labs(y = "Log10(Quantity)", x = "Pathogen", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ #theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  # smaller text & less margin 
  #strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))

combined_plot <- (plot_kp / plot_kh) + plot_annotation(tag_levels = 'A')
# stacked vertically; use | for side-by-side
ggsave("/.../f_k-box-530.pdf", combined_plot, width = 7, height = 6)

# # # # # # # point plots with mean, EBs, and data points

k_data = subset(k_data, subset = Pathogen != "MOCK")
k_data <- k_data %>% group_by(Host, Pathogen) %>%
  mutate(mean_estimate = mean(Quantity, na.rm = TRUE),
         diff_from_mean = abs(Quantity - mean_estimate)) %>%
  arrange(diff_from_mean) %>% slice_head(n = 4) %>% ungroup() %>%
  dplyr::select(-mean_estimate, -diff_from_mean) 

plot_kp <- ggplot(stats, aes(x = Host, y = K, fill = Pathogen, color = Pathogen)) +
  geom_point(size=3) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(K - SE), ymax=(K + SE)))+
  geom_jitter(data = k_data,aes(x = Host, y = Quantity, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1) +
  labs(y = "Log10(Quantity)", x = "Host", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+#theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  
#strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))
plot_kp +
  geom_text(data = asterisks_df, aes(x = Host, y = y, label = label),
            inherit.aes = FALSE, size = 5, vjust = -0.5)

plot_kh <- ggplot(stats, aes(x = Pathogen, y = K, fill = Pathogen, color = Pathogen)) +
  geom_point(size=3) +
  geom_errorbar(width=0.9, size=0.7, aes(ymin=(K - SE), ymax=(K + SE)))+
  geom_jitter(data = k_data,aes(x = Pathogen, y = Quantity, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Host, nrow = 1)+#, strip.position = "right") +
  labs(y = "Log10(Quantity)", x = "Pathogen", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ #theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  theme(strip.text = element_text(size = 10, margin = margin(t = 1, b = 1)))#,  # smaller text & less margin 
#strip.background = element_rect(size = 0.2), panel.spacing = unit(0.1, "lines"))

combined_plot <- (plot_kp / plot_kh) + plot_annotation(tag_levels = 'A')
# stacked vertically; use | for side-by-side
ggsave(".../f_k-point.pdf", combined_plot, width = 6, height = 5)


########

# # # # # # # plots of r & K
########

head(data_all)
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


# data_all_new <- data_all_new %>%
#   mutate(Time = dplyr::recode(as.character(Time), "24" = "1", "72" = "3", "120" = "5"),Time = as.numeric(Time)) 

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

# # # # # # plot k/r/auc 
head(growth_summary)
head(k_data)
k_data$Time = NULL
names(k_data) = c("Host", "Pathogen", "Plate", "K")
merged_data <- merge(growth_summary, k_data[, c("Host", "Pathogen", "Plate", "K")],
                     by = c("Host", "Pathogen", "Plate"),all.x = TRUE)

host_levels <- c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA")

stats <- merged_data %>% group_by(Host, Pathogen) %>% summarise(N = length(auc), K_mean = mean(K),Ksd = sd(K),
             R_mean = mean(slope),Rsd = sd(slope),  AUC_mean = mean(auc),AUCsd = sd(auc), .groups = "drop_last")

plot_r <- ggplot(stats, aes(x = Host, y = R_mean, fill = Pathogen, color = Pathogen)) +
  geom_point(size=3) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(R_mean - Rsd), ymax=(R_mean + Rsd)))+
  geom_jitter(data = merged_data, aes(x = Host, y = slope, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1) +
  labs(y = "Rate", x = "",  color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+#theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_blank())

plot_auc <- ggplot(stats, aes(x = Host, y = AUC_mean, fill = Pathogen, color = Pathogen)) +
  geom_point(size=3) +
  geom_errorbar(width=0.9, size=0.7, aes(ymin=(AUC_mean - AUCsd), ymax=(AUC_mean + AUCsd)))+
  geom_jitter(data = merged_data, aes(x = Host, y = auc, color = Pathogen),
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1)+#, strip.position = "right") +
  labs(y = "AUC", x = "", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ #theme_light() +
  theme(legend.position = "none", 
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 11, margin = margin(t = 1, b = 1))) 

plot_k <- ggplot(stats, aes(x = Host, y = K_mean, fill = Pathogen, color = Pathogen)) +
  geom_point(size=3) +
  geom_errorbar(width=0.9, size=0.7, aes(ymin=(K_mean - Ksd), ymax=(K_mean + Ksd)))+
  geom_jitter(data = k_data, aes(x = Host, y = K, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1)+#, strip.position = "right") +
  labs(y = "K", x = "", color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ #theme_light() + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 10),
        #axis.title.x = element_text(size = 11), 
        #axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_blank())


combined_plot <- (plot_auc/ plot_k / plot_r) + plot_annotation(tag_levels = 'A')
# stacked vertically; use | for side-by-side
#ggsave(".../f_growth_stats.pdf", combined_plot, width = 6, height = 5)

# # # # plot just K and r
# Hosts to shade per pathogen (example)
shade_info <- data.frame(Pathogen = c("A9", "DC3000", "ES4326", "1448A", "NP29"),
  Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"))

# Convert Host names to numeric x-axis positions
shade_info$xpos <- match(shade_info$Host, host_levels)
rects <- shade_info |> transform(xmin = xpos - 0.5,xmax = xpos + 0.5,ymin = -Inf, ymax = Inf)
rects$Pathogen <- factor(rects$Pathogen, levels = levels(stats$Pathogen))


# add astric labels for the significant contrasts
signif_labels <- data.frame(
  Host = c("TOMATO-A9", "A.THALIANA"),
  Pathogen = c("A9", "NP29"),
  y_pos = 0.18, label = "*")
host_levels <- c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA")
pathogen_levels <- c("A9", "DC3000", "ES4326", "1448A", "NP29", "MOCK")
signif_labels$Host <- factor(signif_labels$Host, levels = host_levels)
signif_labels$Pathogen <- factor(signif_labels$Pathogen, levels = pathogen_levels)


plot_r <- ggplot(stats, aes(x = Host, y = R_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(R_mean - (0.5*Rsd)), ymax=(R_mean + (0.5*Rsd))))+
  geom_jitter(data = merged_data, aes(x = Host, y = slope, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1) +
 geom_text(data = signif_labels, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(y = "Growth Rate", x = "",  color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 11))


# # # # plot just r by pathogen and host
plot_r_path <- ggplot(stats, aes(x = Host, y = R_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(R_mean - (0.5*Rsd)), ymax=(R_mean + (0.5*Rsd))))+
  geom_jitter(data = merged_data, aes(x = Host, y = slope, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Pathogen, nrow = 1) +
  geom_text(data = signif_labels, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  #coord_cartesian(ylim = c(0, 0.2)) +
  labs(y = "Growth Rate", x = "",  color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 11))

plot_r_host <- ggplot(stats, aes(x = Pathogen, y = R_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2) + geom_errorbar(width=0.8, size = 0.7, aes(ymin=(R_mean - (0.5*Rsd)), ymax=(R_mean + (0.5*Rsd))))+
  #geom_jitter(data = merged_data, aes(x = Host, y = slope, color = Pathogen), 
  #            size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Host, nrow = 1) +
  #geom_text(data = signif_labels, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  #coord_cartesian(ylim = c(0, 0.2)) +
  labs(y = "Growth Rate", x = "",  color = "Pathogen") +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3")) +
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_text(size = 11))
combined_plot <- (plot_r_path/ plot_r_host) + plot_annotation(tag_levels = 'A')
#ggsave(".../f_r-byPath-byHost.pdf", combined_plot, width = 6, height = 5)




# # # 
stats$Host <- factor(stats$Host, levels = host_levels)
stats$Pathogen <- factor(stats$Pathogen, levels = pathogen_levels)
rects$Host <- factor(rects$Host, levels = host_levels)
rects$Pathogen <- factor(rects$Pathogen, levels = pathogen_levels)
signif_labels$Host <- factor(signif_labels$Host, levels = host_levels)
signif_labels$Pathogen <- factor(signif_labels$Pathogen, levels = pathogen_levels)
# Apply to all datasets
signif_labels <- data.frame(
  Host = c("BEAN", "A.THALIANA"),
  Pathogen = c( "1448A", "NP29"),
  y_pos = 8.2, label = "*")
signif_labels$Host <- factor(signif_labels$Host, levels = host_levels)
signif_labels$Pathogen <- factor(signif_labels$Pathogen, levels = pathogen_levels)


plot_k <- ggplot(stats, aes(x = Host, y = K_mean, fill = Pathogen, color = Pathogen)) +
  geom_rect(data = rects, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = Host),
            fill = "grey93", color = NA, inherit.aes = FALSE) +
  geom_point(size=2) +
  geom_errorbar(width=0.9, size=0.7, aes(ymin=(K_mean - (0.5*Ksd)), ymax=(K_mean + (0.5*Ksd))))+
  geom_jitter(data = k_data, aes(x = Host, y = K, color = Pathogen), 
              size = 0.9, width = 0.2, alpha = 0.4, inherit.aes = FALSE) +
  geom_text(data = signif_labels, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  facet_wrap(~ Pathogen, nrow = 1)+#, strip.position = "right") +
  labs(y = "Log10(K)", x = "Host", color = "Pathogen") +
  coord_cartesian(ylim = c(5, 8.5)) +
  scale_fill_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  scale_color_manual(values = c("A9"="sienna1","DC3000"="skyblue4","ES4326"="violet","1448A"="mediumaquamarine","NP29"="olivedrab3","MOCK" = "gray75"))+
  theme_minimal()+ 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) + 
  theme(strip.text = element_blank())

tight_margin <- theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))  # top, right, bottom, left
plot_k <- plot_k + tight_margin
plot_r <- plot_r + tight_margin

combined_plot <- (plot_r / plot_k) + 
  plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = 'A')

#ggsave(".../f_growth_stats_rk.pdf", combined_plot, width = 6, height = 4)
#########

# # # # # # # stats - 2 way ANCOVA & contrasts

# # # # # # # ANCOVA r
#######
# for Levene, H0 = variances are equal # for Shapiro, H0 = data are normally distributed
growth_summary$logslope = log10(growth_summary$slope)
leveneTest(logslope ~ Pathogen, data = growth_summary)     #  0.06835 . var are equal
shapiro.test(growth_summary$logslope)                      #  1.151e-08 data are not normally distributed
mod <- lm(logslope ~ Pathogen * Host + Plate , data = growth_summary)
anova(mod)
# Response: slope
#               Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Pathogen       4 0.8439 0.21097 18.6455 1.280e-09 ***
# Host           4 1.2778 0.31944 28.2318 1.397e-12 ***
# Plate          1 0.0303 0.03035  2.6823    0.1074    
# Pathogen:Host 16 6.1467 0.38417 33.9525 < 2.2e-16 ***
#  Residuals     53 0.5997 0.01131 
#########

# # # # # # # contrasts r for local adaptation
#########
# # # # # # # planned contrasts for native host vs all other hosts
# A9, ES, NP are sig

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
########

# # # # # # # contrasts r for local dominance
########
# # # # # # # planned contrasts for focal pathogen vs. all other pathogens on each native host
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

# are these different rs?
es_t = subset(radish, subset = Pathogen =="ES4326")
a14_t = subset(radish, subset = Pathogen =="1448A")
dc_t = subset(radish, subset = Pathogen =="DC3000")
temp = rbind(es_t,a14_t,dc_t)
mod= aov(slope~Plate+Pathogen, data = temp)
summary(mod)
#              Df Sum Sq  Mean Sq    F value  Pr(>F)
# Plate        1 0.000398 0.0003980   0.722  0.434
# Pathogen     2 0.001154 0.0005770   1.046  0.417
# Residuals    5 0.002757 0.0005515           

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

es_t = subset(at, subset = Pathogen =="ES4326")
np_t = subset(at, subset = Pathogen =="NP29")
a9_t = subset(at, subset = Pathogen =="A9")
temp = rbind(es_t,np_t, a9_t)
mod = aov(slope~Plate+Pathogen, data = temp)
summary(mod)

########

# # # # # # # contrasts r for each pathogen vs all others 
########
#growth_summary = subset(growth_summary, subset = Pathogen != "MOCK")
#growth_summary <- droplevels(subset(growth_summary, Pathogen != "MOCK"))

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

########

# # # # # # # ANCOVA K
########
# # # # # # 2 way ANCOVA K
# for Levene, H0 = variances are equal # for Shapiro, H0 = data are normally distributed
leveneTest(Quantity ~ Pathogen, data = k_data)     #  0.9878 var are equal
shapiro.test(k_data$Quantity)                      #  0.1751 data are normally distributed

mod <- lm(Quantity ~ Pathogen * Host + Plate , data = k_data)
anova(mod)
# Response: Quantity
#               Df. Sum Sq Mean Sq  F value  Pr(>F)  
# Pathogen        5  71.466 14.2933 12.5590 1.813e-09 ***
# Host            4   1.470  0.3676  0.3230    0.8620    
# Plate           1   0.467  0.4667  0.4101    0.5234    
# Pathogen:Host  20  21.189  1.0594  0.9309    0.5504    
# Residuals     100 113.809  1.1381                              


# # # # # # # planned contrasts for native host vs all other hosts
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

# 1448A - SIG # # # # 
a14.aov <- aov(Quantity ~ Plate + Host, data = a14, contrasts = list(Host = "contr.sum"))
a14.emml <- emmeans(a14.aov, consec ~ Host)
a14$Host=as.factor(a14$Host)
levels(a14$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modela14 <- aov(Quantity~Plate+Host, data = a14)
a14_mod <- emmeans(modela14, specs = ~ Host)
contrast(a14_mod, method = list(c(1,1,1,-4,1))) # 0.0447 # SIG

# NP - SIG
np.aov <- aov(Quantity ~ Plate + Host, data = np, contrasts = list(Host = "contr.sum"))
summary(np.aov )
np.emml <- emmeans(np.aov, consec ~ Host)
np$Host=as.factor(np$Host)
levels(np$Host) # check order of levels "TOMATO-A9"  "TOMATO-DC"  "RADISH"     "BEAN"       "A.THALIANA"
modelnp <- aov(Quantity~Plate+Host, data = np)
np_mod <- emmeans(modelnp, specs = ~ Host)
contrast(np_mod, method = list(c(1,1,1,1,-4))) # 0.0425 # SIG


# # # plot estimated marginal means & contrasts
contrasts_df <- data.frame(
  Pathogen = c("A9", "DC3000", "ES4326", "1448A", "NP29"),
  Estimate = c(-2.56, -0.449, 0.554, -3.26, -3.45),
  SE = c(1.49, 1.29, 1.23, 1.5, 1.58),
  df = c(12, 13, 15, 17, 18),
  t_ratio = c(-1.718, -0.349, 0.452, -2.167, -2.183),
  p_value = c(0.1114, 0.7328, 0.6580, 0.0447, 0.0425))
contrasts_df$Significance <- ifelse(contrasts_df$p_value < 0.05, "Significant", "Not Significant")
contrasts_df$Pathogen <- factor(contrasts_df$Pathogen, levels = c("NP29","1448A","ES4326","DC3000","A9"))
contrasts_df$p_label <- paste0("p = ", formatC(contrasts_df$p_value, format = "f", digits = 3))

contrasts_plot <- ggplot(contrasts_df, aes(y = Pathogen, x = Estimate, color = Pathogen)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Estimate - SE, xmax = Estimate + SE), height = 0.4) +
  geom_text(aes(x = 2.5 + 0.5, label = p_label), hjust = 0, color = "black", size = 3) +
  labs(y = "Pathogen", x = "Estimated Contrast") +
  scale_color_manual(values = c(
    "A9" = "sienna1","DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine","NP29" = "olivedrab3")) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(contrasts_df$Estimate - contrasts_df$SE), 2.5), clip = "off")+
  theme(legend.position = "none", plot.margin = margin(5, 60, 5, 5),
        axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) 
#ggsave(".../f_contrasts.pdf", contrasts_plot, width = 5, height = 2.5)


# estimated marginal means across host
# these won't show the difference unless you plot for all contrasts
# Host    emmean  SE   df  lower.CL upper.CL
# a9_mod  # 6.85 0.339 12     6.11     7.58
# dc_mod  # 6.64 0.294 13     6.00     7.27
# es_mod  # 6.01 0.264 15     5.44     6.57
# a14_mod # 5.75 0.271 17     5.17     6.32
# np_mod  # 6.71 0.357 18     5.96     7.46


df_A9 <- data.frame(
  Pathogen = "A9", Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"),
  emmean = c(6.85, 6.17, 6.05, 6.41, 6.20),
  SE = c(0.339, 0.339, 0.339, 0.269, 0.298),
  df = 12,
  lower.CL = c(6.11, 5.43, 5.31, 5.83, 5.55),
  upper.CL = c(7.58, 6.90, 6.79, 7.00, 6.85))

df_DC <- data.frame(
  Pathogen = "DC3000", Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"),
  emmean = c(6.44, 6.64, 6.43, 6.74, 6.47),
  SE = c(0.228, 0.294, 0.228, 0.294, 0.294),
  df = 13,
  lower.CL = c(5.95, 6.00, 5.94, 6.11, 5.84),
  upper.CL = c(6.94, 7.27, 6.93, 7.38, 7.11))

df_ES <- data.frame(
  Pathogen = "ES4326",Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"),
  emmean = c(5.85, 5.74, 6.01, 6.32, 6.68),
  SE = c(0.296, 0.264, 0.264, 0.296, 0.336),
  df = 15,
  lower.CL = c(5.22, 5.18, 5.44, 5.69, 5.96),
  upper.CL = c(6.48, 6.30, 6.57, 6.95, 7.39))

df_A14 <- data.frame(
  Pathogen = "1448A",Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"),
  emmean = c(5.65, 6.02, 5.75, 6.60, 5.71),
  SE = c(0.246, 0.348, 0.271, 0.348, 0.246),
  df = 17,
  lower.CL = c(5.14, 5.28, 5.17, 5.86, 5.19),
  upper.CL = c(6.17, 6.75, 6.32, 7.33, 6.23))

df_NP <- data.frame(
  Pathogen = "NP29",Host = c("TOMATO-A9", "TOMATO-DC", "RADISH", "BEAN", "A.THALIANA"),
  emmean = c(5.71, 5.95, 5.97, 5.77, 6.71),
  SE = c(0.412, 0.292, 0.357, 0.269, 0.357),
  df = 18,
  lower.CL = c(4.84, 5.34, 5.22, 5.21, 5.96),
  upper.CL = c(6.57, 6.57, 6.72, 6.34, 7.46))

# Combine all into a single data frame
combined_df <- rbind(df_A9, df_DC, df_ES, df_A14, df_NP)
combined_df$Pathogen <- factor(combined_df$Pathogen, levels = c("A9","DC3000","ES4326","1448A","NP29"))
combined_df$Host <- factor(combined_df$Host, levels = c("TOMATO-A9","TOMATO-DC","RADISH","BEAN","A.THALIANA"))

head(combined_df)

# add astric labels for the significant contrasts
signif_labels <- data.frame(Pathogen = c("1448A", "NP29"),Host = c("BEAN", "A.THALIANA"),
  emmean = c(combined_df$emmean[combined_df$Pathogen == "1448A" & combined_df$Host == "BEAN"], 
    combined_df$emmean[combined_df$Pathogen == "NP29" & combined_df$Host == "A.THALIANA"]))
signif_labels$Pathogen <- factor(signif_labels$Pathogen, levels = c("A9", "DC3000", "ES4326", "1448A", "NP29"))
signif_labels$y_pos <- 7.10
#signif_labels$y_pos <- signif_labels$emmean + 0.5


plot_em <- ggplot(combined_df, aes(x = Host, y = emmean, color = Pathogen)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.3, linewidth=0.8) +
  geom_text(data = signif_labels, aes(x = Host, y = y_pos, label = "*"), inherit.aes = FALSE, size = 10, color = "black") +
  facet_grid(~ Pathogen) +
  scale_color_manual(values = c("A9" = "sienna1","DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine","NP29" = "olivedrab3")) +
  labs(x = "Host", y = "Estimated Marginal Mean", color = "Pathogen") +
  theme_minimal() + theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + 
  theme(strip.text = element_text(size = 12))  

#ggsave(".../f_emms-hosts-531.pdf", plot_em, width = 7, height = 3)


# # # # # # # planned contrasts for pathogen vs. all other pathogens across hosts

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

# default plot
temp = plot(contrasts_vs_all)
#ggsave(".../f_temp.pdf", temp, width = 7, height = 3)

# nicer plot
contrasts_df <- as.data.frame(contrasts_vs_all)
contrasts_df$Pathogen <- sub(" effect", "", contrasts_df$contrast)
contrasts_df$p_label <- ifelse(contrasts_df$p.value < 0.001, "<0.001", sprintf("p = %.3f", contrasts_df$p.value))

contrasts_df$Pathogen <- factor(contrasts_df$Pathogen, levels = c("NP29","1448A","ES4326","DC3000","A9"))
pathogen_colors <- c("A9"="sienna1", "DC3000"="skyblue4","ES4326"="violet", "1448A"="mediumaquamarine","NP29"="olivedrab3")

contrast_plot <- ggplot(contrasts_df, aes(y = Pathogen, x = estimate, color = Pathogen)) +
  geom_point(size = 4) + geom_errorbarh(aes(xmin = estimate - SE, xmax = estimate + SE), height = 0.4, linewidth = 0.8) +
  geom_text(aes(x = 0.55, label = p_label), hjust = 0, color = "black", size = 3) +
  scale_color_manual(values = pathogen_colors) + labs(x="Estimated Contrast", y="Pathogen") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.margin = margin(5, 60, 5, 5),
        axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 11),axis.title.x = element_text(size = 11)) +
  coord_cartesian(xlim = c(min(contrasts_df$estimate - contrasts_df$SE), 0.5), clip = "off")
#ggsave(".../f_contrasts-pathogens-across-hosts.pdf", contrast_plot, width = 5, height = 2.5)

######

# # # # # # # contrasts K for local adaptation
######
# radish
rad.aov <- aov(Quantity ~ Plate + Pathogen, data = rad, contrasts = list(Pathogen = "contr.sum"))
rad.emml <- emmeans(rad.aov, consec ~ Pathogen)
rad$Pathogen=as.factor(rad$Pathogen)
levels(rad$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelrad <- aov(Quantity~Plate+Pathogen, data = rad)
rad_mod <- emmeans(modelrad, specs = ~ Pathogen)
contrast(rad_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.9568
contrast(rad_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.0674
contrast(rad_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.9343
contrast(rad_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.2709
contrast(rad_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.4058


# At
at.aov <- aov(Quantity ~ Plate + Pathogen, data = at, contrasts = list(Pathogen = "contr.sum"))
at.emml <- emmeans(at.aov, consec ~ Pathogen)
at$Pathogen=as.factor(at$Pathogen)
levels(at$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelat <- aov(Quantity~Plate+Pathogen, data = at)
at_mod <- emmeans(modelat, specs = ~ Pathogen)
contrast(at_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.7160
contrast(at_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.6758
contrast(at_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.4202
contrast(at_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.0322 **
contrast(at_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.3152

# bean
bean.aov <- aov(Quantity ~ Plate + Pathogen, data = bean, contrasts = list(Pathogen = "contr.sum"))
bean.emml <- emmeans(bean.aov, consec ~ Pathogen)
bean$Pathogen=as.factor(bean$Pathogen)
levels(bean$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelbean <- aov(Quantity~Plate+Pathogen, data = bean)
bean_mod <- emmeans(modelbean, specs = ~ Pathogen)
contrast(bean_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.7756
contrast(bean_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.2273
contrast(bean_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.9852
contrast(bean_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.4858
contrast(bean_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.0329 **

# Tom A9
ta9.aov <- aov(Quantity ~ Plate + Pathogen, data = ta9, contrasts = list(Pathogen = "contr.sum"))
ta9.emml <- emmeans(ta9.aov, consec ~ Pathogen)
ta9$Pathogen=as.factor(ta9$Pathogen)
levels(ta9$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelta9 <- aov(Quantity~Plate+Pathogen, data = ta9)
ta9_mod <- emmeans(modelta9, specs = ~ Pathogen)
contrast(ta9_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.0368**
contrast(ta9_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.2285
contrast(ta9_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.3810
contrast(ta9_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.1089
contrast(ta9_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.2469

# Tom DC
tdc.aov <- aov(Quantity ~ Plate + Pathogen, data = tdc, contrasts = list(Pathogen = "contr.sum"))
tdc.emml <- emmeans(tdc.aov, consec ~ Pathogen)
tdc$Pathogen=as.factor(tdc$Pathogen)
levels(tdc$Pathogen) #"dc"     "DC3000" "ES4326" "1448A"  "NP29"   
modeltdc <- aov(Quantity~Plate+Pathogen, data = tdc)
tdc_mod <- emmeans(modeltdc, specs = ~ Pathogen)
contrast(tdc_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.8689
contrast(tdc_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.0023**
contrast(tdc_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.0279**
contrast(tdc_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.5324
contrast(tdc_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.1511
######

# # # # # # # contrasts K for local dominance
######
# # # # # # # planned contrasts for focal pathogen vs. all other pathogens on each native host
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

np_t = subset(at, subset = Pathogen =="NP29")
es_t = subset(at, subset = Pathogen =="ES4326")
dc_t = subset(at, subset = Pathogen =="DC3000")
temp = rbind(np_t,es_t,dc_t)
test = aov(Quantity~Plate+Pathogen, data = temp) # no effect of pathogen on K
summary(test)
# > summary(test)
# Df Sum Sq Mean Sq F value Pr(>F)
# Plate        1 0.0198  0.0198   0.057  0.823
# Pathogen     1 0.0006  0.0006   0.002  0.969
# Residuals    4 1.3925  0.3481   

test = t.test(es_t$Quantity,np_t$Quantity) # p = 0.9859 no difference
test = t.test(es_t$Quantity,dc_t$Quantity) # p = 0.8087 no difference



######

# # # # # # # contrasts K for each pathogen vs all others 
######

rad = subset(k_data, subset = Host =="RADISH")
at = subset(k_data, subset = Host =="A.THALIANA")
bean = subset(k_data, subset = Host =="BEAN")
ta9 = subset(k_data, subset = Host =="TOMATO-A9")
tdc = subset(k_data, subset = Host =="TOMATO-DC")

# radish
rad.aov <- aov(Quantity ~ Plate + Pathogen, data = rad, contrasts = list(Pathogen = "contr.sum"))
rad.emml <- emmeans(rad.aov, consec ~ Pathogen)
rad$Pathogen=as.factor(rad$Pathogen)
levels(rad$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelrad <- aov(Quantity~Plate+Pathogen, data = rad)
rad_mod <- emmeans(modelrad, specs = ~ Pathogen)
contrast(rad_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.9568
contrast(rad_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.0674
contrast(rad_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.9343
contrast(rad_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.2709
contrast(rad_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.4058


# At
at.aov <- aov(Quantity ~ Plate + Pathogen, data = at, contrasts = list(Pathogen = "contr.sum"))
at.emml <- emmeans(at.aov, consec ~ Pathogen)
at$Pathogen=as.factor(at$Pathogen)
levels(at$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelat <- aov(Quantity~Plate+Pathogen, data = at)
at_mod <- emmeans(modelat, specs = ~ Pathogen)
contrast(at_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.7160
contrast(at_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.6758
contrast(at_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.4202
contrast(at_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.0322 **
contrast(at_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.3152

# bean
bean.aov <- aov(Quantity ~ Plate + Pathogen, data = bean, contrasts = list(Pathogen = "contr.sum"))
bean.emml <- emmeans(bean.aov, consec ~ Pathogen)
bean$Pathogen=as.factor(bean$Pathogen)
levels(bean$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelbean <- aov(Quantity~Plate+Pathogen, data = bean)
bean_mod <- emmeans(modelbean, specs = ~ Pathogen)
contrast(bean_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.7756
contrast(bean_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.2273
contrast(bean_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.9852
contrast(bean_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.4858
contrast(bean_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.0329 **

# Tom A9
ta9.aov <- aov(Quantity ~ Plate + Pathogen, data = ta9, contrasts = list(Pathogen = "contr.sum"))
ta9.emml <- emmeans(ta9.aov, consec ~ Pathogen)
ta9$Pathogen=as.factor(ta9$Pathogen)
levels(ta9$Pathogen) #"A9"     "DC3000" "ES4326" "1448A"  "NP29"   
modelta9 <- aov(Quantity~Plate+Pathogen, data = ta9)
ta9_mod <- emmeans(modelta9, specs = ~ Pathogen)
contrast(ta9_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.0368**
contrast(ta9_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.2285
contrast(ta9_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.3810
contrast(ta9_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.1089
contrast(ta9_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.2469

# Tom DC
tdc.aov <- aov(Quantity ~ Plate + Pathogen, data = tdc, contrasts = list(Pathogen = "contr.sum"))
tdc.emml <- emmeans(tdc.aov, consec ~ Pathogen)
tdc$Pathogen=as.factor(tdc$Pathogen)
levels(tdc$Pathogen) #"dc"     "DC3000" "ES4326" "1448A"  "NP29"   
modeltdc <- aov(Quantity~Plate+Pathogen, data = tdc)
tdc_mod <- emmeans(modeltdc, specs = ~ Pathogen)
contrast(tdc_mod, method = list(c(-4,1,1,1,1))) # A9 vs all others 0.8689
contrast(tdc_mod, method = list(c(1,-4,1,1,1))) # DC vs all others 0.0023**
contrast(tdc_mod, method = list(c(1,1,-4,1,1))) # ES vs all others 0.0279**
contrast(tdc_mod, method = list(c(1,1,1,-4,1))) # 14 vs all others 0.5324
contrast(tdc_mod, method = list(c(1,1,1,1,-4))) # NP vs all others 0.1511

######


