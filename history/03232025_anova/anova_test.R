
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


# https://www.r-bloggers.com/2020/09/how-to-convert-continuous-variables-into-categorical-by-creating-bins/

dados<-dados%>%mutate(mminBins = cut(mmin, 
                                     breaks = c(-Inf,20,40,60,80,100,Inf)))

#dados<-dados%>%mutate(minwindowsizeBins = cut(minwindowsize  , 
#                                     breaks = c(-Inf,20,40,60,80,100,Inf)))
number_of_bins = 12
dados<-dados%>%mutate(minwindowsizeBins = cut(minwindowsize, 
        breaks = unique(quantile(minwindowsize,
                                 probs=seq.int(0,1, by=1/number_of_bins))), 
          include.lowest=TRUE))


#col_names <- names(dados[,c(1,3,4,5,6,7)])
col_names <- names(dados[,c(5,6,7,13,14)])
dados[col_names] <- lapply(dados[col_names] , factor)
dados = dados[,c(col_names,"fp")]
summary(dados)



#model <-aov(fp ~ mmin+minwindowsize+ntests+eegchannel+stimlvl , data = dados)
model <-aov(fp ~ ., data = dados)
summary(model)

#mc1 <- glht(model, linfct = mcp(stimlvl = "Tukey"))
mc1 <- glht(model, linfct = mcp(minwindowsizeBins = "Tukey"))
summary(mc1)

# Calcula os intervalos de confiança para as diferenças entre os grupos 
# com nível de confiança de 95%
mc1_CI <- confint(mc1, level = 0.95)

# Plota os intervalos de confiança para cada comparação
# O gráfico ajuda a visualizar quais pares de grupos 
# têm diferenças significativas

png(filename="figure4.png", 
    width=800,height = 800, bg="white")
par(mar=c(5,6,4,1)+5.5)
plot(mc1_CI)
dev.off()
#minwindowsize
#mc1 <- glht(model, linfct = mcp(minwindowsize = "Tukey"))
#summary(mc1)

# Calcula os intervalos de confiança para as diferenças entre os grupos 
# com nível de confiança de 95%
#mc1_CI <- confint(mc1, level = 0.95)

#mc1 <- glht(model, linfct = mcp(mmin = "Tukey"))
#summary(mc1)


#hist(model$residuals)
shapiro.test(model$residuals)
qqPlot(model$residuals)















