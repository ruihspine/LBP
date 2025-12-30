library(tidyverse)
library(readxl)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) 
library(mapdata)
rm(list = ls())
df=read.csv("QN eapc.csv",header = T)
load("./data/GBD.Rdata")
colnames(df)

labelx="Incidence"


## Filter data

dfx<- df %>%
  filter(measure_name == labelx) %>% 
  filter(sex_name == "Both") %>% 
  #filter(age_id == 27) %>% 
  filter(metric_name == "Rate") %>% 
  filter(cause_name == "Low back pain")

print(dim(dfx))

####  map for ASMR
################################################################################################
dfplot= dfx  %>% filter(location_id %in% namex$location_id ) %>% 
  select(location_id,val) %>% left_join(.,namex)
color1 <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
            "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac")

colorx1=color1


################################################################################################
########################
# 1.1 Prevalence
########################
# Prevalence
index="Incidence"
ASR=dfplot %>%  select(location_id,location,val) %>% 
  mutate(val=val/1)
df_asr=left_join(df_world,ASR)


xmin=min(na.omit(df_asr$val))
# Calculate the quantile breaks
breaks <- quantile(ASR$val, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
breaks=breaks[-1]
print(breaks)
breaks100=1*breaks
# Initialize an empty vector to store the formatted breaks
formatted_breaks = vector("character", length(breaks100))
# Set the first break manually to start from 0
formatted_breaks[1] = paste(paste0(round(xmin,3)," -"), format(breaks100[1], digits = 4, nsmall = 1))
# Iterate over the breaks and format them
for (i in 2:length(breaks100)) {
  formatted_breaks[i] = paste(paste0(format(breaks100[i - 1], digits = 4, nsmall = 1)), "-",paste0( format(breaks100[i], digits = 4, nsmall = 1)))
}
# Print the formatted breaks
#print(formatted_breaks)
# Print the formatted breaks
unb=unique(breaks)
if(unb[1]==xmin){
  unb=unb[-1]
} else{
  unb=unb
}

print(formatted_breaks)
x=cut(df_asr$val, breaks = c(xmin, unb))
print(na.omit(unique(x)))
xlen=length(na.omit(unique(x)))
xstar=length(formatted_breaks)-xlen+1
formatted_breaks=formatted_breaks[xstar:length(formatted_breaks)]
# Cut break for Prevalence
df_asr = df_asr %>%
  mutate(asr = replace_na(val, -99)) %>%
  mutate(asr_cut = cut(asr, breaks = c(xmin, unb), 
                       labels = formatted_breaks)) %>% 
  filter(!is.na(asr_cut))



ggplot(df_asr %>% na.omit()) +
  geom_sf(aes(geometry = geometry, fill = asr_cut),size = 0.1)+
  #labs(title = paste0("Age-standardized ",index," Rate (Per 100,000),both sexes in 2017")) +
  scale_fill_manual(#name="Anaemia prevalence",
    values = rev(colorx1),
    guide = guide_legend(reverse=T))+
  guides(fill = guide_legend(ncol = 2, title = labelx))->p

p1=p+ theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = c(0.13, 0.29),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.title=element_text(size=8),
            legend.text=element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
p1


p2=p1+theme(legend.position = 'top')
p2

ggsave("map.pdf",width = 12,height = 10)