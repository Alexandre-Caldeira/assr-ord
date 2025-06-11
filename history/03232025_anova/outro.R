
# Limpar o ambiente, gráficos e console
rm(list = ls())       # Remove todas as variáveis do ambiente
graphics.off()        # Fecha todos os gráficos abertos
cat("\014")           # Limpa o console

setwd("C:/PPGEE/Assessing CGST on ASSR/clean_code/assr-ord/history/23032025_anova")

# Carregar pacotes necessários
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(car)
library(multcomp)

dados1   <- read.table(file   = "30db_cgst_fix.csv", 
                       header = TRUE, 
                       sep    = ",")
dados2   <- read.table(file   = "50db_cgst_fix.csv", 
                       header = TRUE, 
                       sep    = ",")

dados <- rbind(dados1,dados2)


# Define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

z_norm <- function(x) {
  (x - mean(x)) / (sd(x))
}

# Apply Min-Max normalization to the first four columns in the iris dataset
#dados[c(1,2,3,4,5,9,10,11,12)] <- as.data.frame(lapply(dados[c(1,2,3,4,5,9,10,11,12)], min_max_norm))
dados[c(1,2,3,4,5)] <- as.data.frame(lapply(dados[c(1,2,3,4,5)], z_norm))




col_names <- names(dados[,c(1,3,4,7,10)])
#col_names <- names(dados[,c(1,3,4,5,6,7)])
#col_names <- names(dados[,c(5,6)])
#col_names <- names(dados[,c(4,6)])
dados[col_names] <- lapply(dados[col_names] , factor)
dados = dados[,c(col_names,"tp")]
summary(dados)

#model <-aov(fp ~ mmin+minwindowsize+ntests+eegchannel+stimlvl , data = dados)
model <-aov(tp ~ ., data = dados)
summary(model)

#mc1 <- glht(model, linfct = mcp(stimlvl = "Tukey"))
#mc1 <- glht(model, linfct = mcp(eegchannel = "Tukey"))
#summary(mc1)

#mc1_CI <- confint(mc1, level = 0.95)



#par(mar=c(5,6,4,1)+2.5)
#plot(mc1_CI)
shapiro.test(model$residuals)
qqPlot(model$residuals)
