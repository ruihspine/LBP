library(tidyverse)
library(patchwork)
rm(list = ls())
load("./data/GBD.Rdata")
load("./data/population231.Rdata")
source("source_apc.R")
source("function_year5.R")
df=read.csv("QN APC.csv",header = T)
colnames(df)

######## ######## ######## ######## ######## 
## Table1
df1=df %>% filter(measure_name=="Incidence") %>% 
  filter(location_id %in% c(1)) %>% 
  #filter(cause_name=="Pancreatic cancer") %>% 
  filter(metric_name=="Number") %>% 
  filter(sex_name=="Both") %>% 
  filter( age_id %in% c(1,5:20,30,31,32,235))


dfx=df1 %>% 
  select(age_id,age_name,year,val) %>% 
  spread(year,val) %>% 
  separate(age_name,sep=' ', c("age_name"))


dfx1=dfx %>% select(-1,-2)
rownames(dfx1)=dfx$age_name

#转成5年一组, 从21年开始往前每5年一个组
dfx2 <- function_year5(dfx1, 1990, 2023, 2023)
rownames(dfx2)=dfx$age_name

population=read.csv("population.csv",header = T)

dfpop=population %>% dplyr::filter(location_id %in% 1 )  %>% 
  #filter(cause_name%in% input$tablecause_name10) %>% 
  filter(sex_name %in% "Both") %>% 
  filter( age_id %in% c(1,5:20,30,31,32,235)) %>% 
  filter(metric_name %in% "Number") %>%
  select(age_id,age_name,year,val) %>% 
  spread(year,val) %>% 
  separate(age_name,sep=' ', c("age_name")) %>% as.data.frame()

dfpop1=dfpop %>% select(-1,-2)
rownames(dfpop1)=dfpop$age_name
#转成5年一组, 从21年开始往前每5年一个组
dfpop2 <- function_year5(dfpop1, 1990, 2023, 2023)
rownames(dfpop2)=dfpop$age_name
####合并数据，number+pop
name2 = paste0(names(dfpop2),"p")
dfpop2 <- dfpop2%>% stats::setNames(name2)
print(dfx2 %>% head())

dfx = tibble(cbind(dfx2,dfpop2)) %>% 
  dplyr::select(`1990-1993`,`1990-1993p`,`1994-1998`,`1994-1998p`,`1999-2003`,`1999-2003p`,
                `2004-2008`,`2004-2008p`,`2009-2013`,`2009-2013p`,`2014-2018`,`2014-2018p`,
                `2019-2023`,`2019-2023p`)

#保存表格，用APC网页工具完成APC分析，具体见文章原文。

write.csv(dfx, file = "apc qn")
write.csv(dfx, file = "apc qn.dfx")
