# Rebecca Satterwhite
# script to check cfu/mL at different ods

library(dplyr) 
library(ggplot2)

data = read.table(".../cfu_in.txt", header=T)
dim(data) # 83  5
head(data)


# cfu/ml = counts/(1/df)*vol plated
data[6] = (data$count/(1/data$df))*data$ml_plated
data[6] = log10(data[6])
head(data)

namevec = c("strain","od","count","df","ml_plated","cfuml")
names(data) =namevec
head(data)


# remove 2 ES obs = 0 (drop as plating errors)
#test = subset(data, subset = strain=="ES4326")
dim(data) # 75  6
data=subset(data, subset = cfuml > 1)
dim(data) # 73  6
data=subset(data, subset = count <2000)
dim(data) # 71  6

head(data)

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


pdf(file = ".../fig_CFU-per-OD.pdf", width = 4, height = 3)
dodge =position_dodge(.3) # how much jitter on the x-axis?
ggplot(data=stats, aes(x=od, y=Mean, color=strain)) + 
  geom_point(size=2, position = dodge) + 
  geom_errorbar(aes(x=od, ymin=Mean-(1/2*CI95), ymax=Mean+(1/2*CI95)),position = dodge) +
  labs(y = "log10(CFU/mL)", color="Pathogen") +
  scale_x_discrete(name ="OD600", limits=c("0.01","0.1","1")) +
  scale_color_manual(values = c("A9" = "sienna1","DC3000" = "skyblue4","ES4326" = "violet","1448A" = "mediumaquamarine", "NP29" = "olivedrab3")) +
  annotate("rect", xmin = -Inf, xmax = "0.1", ymin = 8.1, ymax = Inf, fill = "white", color = NA) +
  annotate("text", x = "0.01", y = max(stats$Mean) + 0.5, label = "1-Way ANOVA", hjust = 0.55, vjust = 1,size = 3) +
  annotate("text", x = "0.01", y = max(stats$Mean) + 0.1,label = "F==0.568*','~italic(p)==0.687", hjust = 0.4, vjust = 0.7,parse = TRUE,size = 3)+
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray80", size = 0.5), panel.grid.minor = element_blank())
dev.off()


# ANOVA
test = (aov(cfuml~strain, data))
summary(test)
#             Df Sum Sq  Mean Sq  F value  Pr(>F)
# strain       4   3.85  0.9623   0.568    0.687
# Residuals   82 138.92  1.6941  



