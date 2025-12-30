library(tidyverse)
library(patchwork)
rm(list = ls())

df=read.csv("QN 204 region 1990-2023.csv",header = T)
colnames(df)


dfx=df %>% filter(measure_name=="Incidence") %>%   ##可依次换 YLDs (Years Lived with Disability)，Prevalence

  filter(sex_name=="Both") %>% 
  filter(metric_name=="Rate")

print(dim(dfx))



  print(unique(dfx$measure_name))
  
  
  
  
  
  
  ######### 1.1  Incidence SDI 2019 in 204
  ############### ############### ############### ############### ############### 
  GBD1990_2019 = dfx %>%
    select(location_id,val,year)
  
  SDI2023=read.csv("SDI 2023.csv",header = T)
  
  dim(SDI2023)
  
  SDI2019 = SDI2023 %>% 
    filter(year==2023) %>% 
    left_join(.,GBD1990_2019) %>% 
    filter(!is.na(val))
  # link
  dim(SDI2019)
  
  library(ggrepel)
  ## plot
  p1=ggplot(SDI2019, aes(x=sdi, y=val)) + 
    geom_point(aes(col=location_name),size=0.8) + 
    geom_smooth(method="loess",se = T,color="#708090") +
    geom_text_repel(aes(label = location_name,col=location_name),size=3.2, max.overlaps = 60, segment.size = 0.2)+
    labs(#title = paste0(index," ASR and SDI in National level"),
      x=paste0("Socio-demographic Index"),
      y=paste0("Age-standardized prevalence rate per 100,000 population"))+
    theme(legend.position = "none")+
    scale_x_continuous(name="SDI ",limits = c(0.2,1.0),breaks =  seq(0,1,0.1))+
    theme_bw() +
    theme(legend.position = "none",
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  p1
  ggsave("111111.tiff",
         plot = p1,
         width = 8.7,    # 常见期刊单栏宽度
         height = 6.5,
         dpi = 600,      # 高分辨率
         compression = "lzw")  # TIFF压缩
  
  ggsave(paste0("QN  204 SDI.pdf"),width = 10,height = 8)
 
  
   # Calculate Spearman's correlation
  spearman_cor <- cor.test(SDI2019$sdi, SDI2019$val, method = "spearman")
  # Extract the correlation coefficient (R) and p-value
  r <- spearman_cor$estimate
  p <- spearman_cor$p.value
  
  # Print the results
  cat("Spearman's Correlation Coefficient (R):", r, "\n")
  cat("p-value:", p, "\n")
  # Create the combined string using sprintf()
  spearx <- sprintf("r=%.4f, p=%.3e", r, p)
  
  

}

