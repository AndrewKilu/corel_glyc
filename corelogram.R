library(readr)
library(readxl)
library(dplyr)
library(Hmisc)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(Hmisc)

results_type <- read.csv("slc34a2_and_glic_AE (1).csv", sep = ";")
  
ae_cancer <- results_type %>%
  filter(tissue_type == "cancer") %>%
  filter(Characteristics.organism.part. == "ovary") %>%
  dplyr::select(-5, -tissue_type) 


ae_healthy <- results_type %>%
  filter(tissue_type == "normal") %>%
  filter(Characteristics.organism.part. == "ovary") %>% 
  dplyr::select(-5, -tissue_type)

convert_number <- function(x){
  x <- as.character(x)
  x <- gsub(pattern = ",", replacement = ".",x = x, fixed = TRUE)
  x <- as.numeric(x)
  return(x)
}

ae_cancer$Expression <- convert_number(ae_cancer$Expression) 
ae_healthy$Expression  <- convert_number(ae_healthy$Expression) 

ae_cancer$Expression  <- as.integer(ae_cancer$Expression)
ae_healthy$Expression  <- as.integer(ae_healthy$Expression)

ae_cancer <- ae_cancer%>% 
  mutate_if(is.numeric, round) 
ae_healthy <- ae_healthy%>% 
  mutate_if(is.numeric, round) 


ae_cancer <- pivot_wider(ae_cancer, names_from = sample, values_from =  Expression )
ae_healthy <- pivot_wider(ae_healthy, names_from = sample, values_from =  Expression )

ae_cancer <- as.data.frame(ae_cancer)

for (i in 2:ncol(ae_cancer)){
  for (m in 1:nrow(ae_cancer)){
    ae_cancer[m,i] <-  convert_number(round((mean(unlist(ae_cancer[m,i])))))
  }
}

ae_healthy <- as.data.frame(ae_healthy)

for (i in 2:ncol(ae_healthy)){
  for (m in 1:nrow(ae_healthy)){
    ae_healthy[m,i] <-  convert_number(round((mean(unlist(ae_healthy[m,i])))))
  }
}

row.names(ae_cancer) <- ae_cancer$HGNC.symbol
ae_cancer[, c(2:534)] <- as.integer(unlist(ae_cancer[, c(2:534)])) 
ae_cancer <- dplyr::select(ae_cancer, -HGNC.symbol)
#на этом этапе все готтово для корреляции, получения коррелограммы 


row.names(ae_healthy) <- ae_healthy$HGNC.symbol
ae_healthy[, c(2:41)] <- as.integer(unlist(ae_healthy[, c(2:41)])) 
ae_healthy <- dplyr::select(ae_healthy, -HGNC.symbol)



library(corrplot)
df <- ae_cancer
df_cor <- cor(df)

corrplot(df_cor)
#слишком громоздко получается, поэтому используем функцию
corr_simple <- function(data=df,sig=0.5){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ", type="upper", tl.cex = 0.2)
}
corr_simple()

file_path= "cor_can.png"
png(height=1800, width=1800, file=file_path, type = "cairo")
corr_simple()

dev.off()



df <- ae_healthy
df_cor <- cor(df_filtered)
# Исключить столбцы с нулевым стандартным отклонением
df_filtered <- df[, sapply(df, function(x) sd(x, na.rm = TRUE) != 0)]

corrplot(df_cor)

corr_simple <- function(data=df,sig=0.5){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ", type="upper", tl.cex = 0.2)
}
corr_simple()


