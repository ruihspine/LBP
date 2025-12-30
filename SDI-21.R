library(tidyverse)
library(patchwork)
rm(list = ls())
load("./data/GBD.Rdata")
df=read.csv("QN 21region 1990-2023.csv",header = T)
colnames(df)



######## ######## ######## ######## ######## 
## Table1
dfx=df %>% filter(measure_name=="Prevalence") %>%    ##!!修改指标 incidence YLDs
  filter(location_id %in% c(1,name_21)) %>% 
  #filter(year %in% c(2021)) %>% 
  filter(sex_name=="Both") %>% 
  filter(#cause_name=="Liver cancer",
    metric_name=="Rate"
  ) 
#%>% 
  #filter(age_id %in% c(27))

dim(dfx)


# Define the specific order for the location names



location_order <- c(
  "Global", "High-income Asia Pacific", "High-income North America", "Western Europe",
  "Australasia", "Andean Latin America", "Tropical Latin America", "Central Latin America",
  "Southern Latin America", "Caribbean", "Central Europe", "Eastern Europe", "Central Asia",
  "North Africa and Middle East", "South Asia", "Southeast Asia", "East Asia", "Oceania",
  "Western Sub-Saharan Africa", "Eastern Sub-Saharan Africa", "Central Sub-Saharan Africa",
  "Southern Sub-Saharan Africa"
)
library(ggsci)
colors <- c(  pal_npg("nrc", alpha = 0.7)(9),
              pal_aaas("default", alpha = 0.7)(9),
              pal_nejm("default", alpha = 0.7)(8),
              pal_jama("default", alpha = 0.7)(7))
######### 1.1  Incidence SDI 2019 in 21
############### ############### ############### ############### ############### 
#bin data
DALY_2017=dfx %>% select(location_id,val,year)
SDI2019=read.csv("SDI 2023.csv",header = T)


df11 = left_join(DALY_2017,SDI2019) 
# get df
df3=df11 %>% filter(sdi>0.2)
#plot
# Filter out rows with NA in 'sdi' or 'val'
df11_filtered <- df11 %>% 
  filter(!is.na(sdi) & !is.na(val)) %>%
  mutate(location_name = factor(location_name, levels = location_order))

ggplot(df11_filtered, aes(x = sdi, y = val, color = location_name, shape = location_name)) +
  geom_point() +  # Plot points
  geom_smooth(method = "loess", se = T, aes(group = 1), color = "#92A8D1") +  # Add a smoothing line
  scale_shape_manual(values = 1:22,breaks = location_order ,labels = location_order) +  # Use manually defined shapes
  scale_color_manual(values = colors,breaks = location_order ,labels = location_order) +
  labs(x = "SDI", 
       shape="",color="",x=paste0("SDI"),
       y=paste0("labx"," Rate per 100,000 population"))+
  theme_bw() +
  theme(legend.key.size = unit(0.03,"line"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_blank(),
        legend.title=element_text(size=12),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=8)) +
  guides(shape = guide_legend(nrow = 22))


# Calculate Spearman's correlation
spearman_cor <- cor.test(df11_filtered$sdi, df11_filtered$val, method = "spearman")
# Extract the correlation coefficient (R) and p-value
r <- spearman_cor$estimate;p <- spearman_cor$p.value
# Print the results
cat("Spearman's Correlation Coefficient (R):", r, "\n")
cat("p-value:", p, "\n")
# Create the combined string using sprintf()
spearx <- sprintf("r=%.4f, p=%.3e", r, p)



