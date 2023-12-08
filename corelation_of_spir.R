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
  select(-5, -tissue_type)

ae_healthy <- results_type %>%
  filter(tissue_type == "normal") %>%
  filter(Characteristics.organism.part. == "ovary") %>% 
  select(-5, -tissue_type)

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


ae_cancer <- pivot_wider(ae_cancer, names_from = HGNC.symbol, values_from =  Expression )
ae_healthy <- pivot_wider(ae_healthy, names_from = HGNC.symbol, values_from =  Expression )

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



ae_cancer[, c(2:245)] <- as.integer(unlist(ae_cancer[, c(2:245)])) 
ae_cancer <- select(ae_cancer, -sample)
#на этом этапе все готтово для корреляции, получения коррелограммы 


ae_healthy[, c(2:245)] <- as.integer(unlist(ae_healthy[, c(2:245)])) 
ae_healthy <- select(ae_healthy, -sample)


cormat_c <- rcorr(as.matrix(ae_cancer, method = "spearman"))
cormat_h <- rcorr(as.matrix(ae_healthy, method = "spearman"))

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

my_cor_matrix_c <- flat_cor_mat(cormat_c$r, cormat_c$P)
my_cor_matrix_h <- flat_cor_mat(cormat_h$r, cormat_h$P)


#умеренная корреляция между всеми 
my_cor_matrix_filter_c <- my_cor_matrix_c %>%
  filter(p < 0.05) %>% 
  filter(cor < -0.4 | cor > 0.4) 
my_cor_matrix_filter_h <- my_cor_matrix_h %>%
  filter(p < 0.05) %>% 
  filter(cor < -0.4 | cor > 0.4) 

#корреляция именно с SLC34A2
list_for_glicslc34a2_cancer <- my_cor_matrix_c %>%
  filter(row == "SLC34A2") %>% 
  filter(p < 0.05) %>% 
  filter(cor < -0.4 | cor > 0.4) 

list_for_glicslc34a2_health <- my_cor_matrix_h %>% 
  filter(row == "SLC34A2") %>% 
  filter(p < 0.05) %>% 
  filter(cor < -0.4 | cor > 0.4) 

