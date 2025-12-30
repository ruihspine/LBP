library(tidyverse)
library(scales)
df=read.csv("QN trend prevalence rate.csv",header = T)
colnames(df)


######## ######## ######## ######## ######## 
## Table1
df1=df %>% filter(measure_name=="Prevalence") %>% 
  filter(location_id %in% c(44634,44635,44636,44637,44639,1)) %>% 
  #filter(cause_name=="Liver cancer") %>% 
  #filter(age_id==27) %>% 
  filter(metric_name=="Rate") %>% 
  filter(sex_name=="Both") 

dfx=df1 %>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = location_name)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = location_name), alpha = 0.2, show.legend = T) + 
  geom_line(aes(color = location_name,linetype=location_name)) +
  geom_point(aes(color = location_name,linetype=location_name)) +
  #facet_wrap(~ sex_name, scales = 'free') +
  theme_classic()



ggplot(dfx, aes(x = year, y = val, fill = location_name)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = location_name), alpha = 0.2, show.legend = T) + 
  geom_line(aes(color = location_name,linetype=location_name)) +
  #facet_wrap(~ sex_name, scales = 'free') +
  theme_classic() +
  scale_y_continuous(labels = label_number(unit = "K")) +
  scale_x_continuous(breaks = seq(min(dfx$year, na.rm = TRUE), max(dfx$year, na.rm = TRUE), by = 3)) + 
  labs(x = "Year", y ="Age standard", color = "", linetype = "") +
  guides(color = guide_legend(title = "", nrow = 1), 
         linetype = guide_legend(title = "", nrow = 1)) +
  theme(strip.background = element_blank(),
        legend.position="top")




### 1. 按照不同地区分面
dfx=df1 #%>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = location_name)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = location_name), alpha = 0.2, show.legend = T) + 
  geom_line(aes(color = location_name,linetype=location_name)) +
  facet_wrap(~ location_name, scales = 'free') +
  theme_classic() +
  scale_y_continuous(labels = label_number(unit = "K")) +
  scale_x_continuous(breaks = seq(min(dfx$year, na.rm = TRUE), max(dfx$year, na.rm = TRUE), by = 5)) + 
  labs(x = "Year", y ="value", color = "", linetype = "") +
  guides(color = guide_legend(title = "", nrow = 1), 
         linetype = guide_legend(title = "", nrow = 1)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "top")




### 1. 按照不同地区分面
dfx=df1 #%>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = location_name)) +
  geom_histogram(stat = "identity", position = "dodge") +
  facet_wrap(~ location_name, scales = 'free') +
  theme_classic() 



df1=df %>% filter(measure_name=="Prevalence") %>% 
  filter(location_id %in% c(44634,44635,44636,44637,44639,1)) %>% 
  #filter(cause_name=="Liver cancer") %>% 
 # filter(age_id==27) %>% 
  filter(metric_name=="Rate") 

dfx=df1 #%>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = sex_name)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ location_name, scales = 'free') +
  theme_classic() 


## change the color
dfx=df1 #%>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = sex_name)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ location_name, scales = 'free') +
  theme_classic() +
  scale_fill_manual(values = c("#e41a1c", "#377eb8","grey")) 


## change ggsci
dfx=df1 #%>% filter(location_id==1)
ggplot(dfx, aes(x = year, y = val, fill = sex_name)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ location_name, scales = 'free') +
  theme_classic() +
  ggsci::scale_fill_nejm()


ggsave("p1.pdf",width = 8,height = 6,dpi = 300)





