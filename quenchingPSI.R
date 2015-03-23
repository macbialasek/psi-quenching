setwd("C:/Users/maciek/Desktop/Studia dok/NP/Chlorophyll fluorescence/Quenching/Fluorcam/13_2_15")
library(plyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(broom)

fl <- list.files(pattern=".TXT")

df = data.frame(
  "param" = NA,
  "plant"= NA,
  "value"= NA,
  "read"=NA,
  "conc"=NA,
  check.names= FALSE
)

for(i in 1:length(fl)){
  temp <- as.data.frame(read.table(fl[i], skip=3, header=F))
  temp <- rbind(temp[42,],temp[47,],temp[52,], temp[64,])
  temp <- melt(temp, id.vars=("V1"))
  temp$read <- str_sub(fl[i], end=4)
  temp$conc <- paste(str_sub(fl[i], start=6, end=7), "mg/l", sep=" ")
  colnames(temp)=c("param", "plant", "value", "read", "conc")
  df = rbind(df, temp)
}
df = df[2:nrow(df),]

#################

con10 <- subset(df, conc=="10 mg/l")
con10$read = factor(con10$read, levels= c("Kont", "Au13", "Au30", "Au40", "Au50"))
con10split <- split(con10, con10$param)

#Test TukeyHSD dla każdego parametru fluorescencji, stężenie 10 ppm:

tk_qmax_c10 <-TukeyHSD(aov(value~read, data = con10split[[4]]))
tk_fvfm_c10 <-TukeyHSD(aov(value~read, data = con10split[[1]]))
tk_qlss_c10 <-TukeyHSD(aov(value~read, data = con10split[[3]]))
tk_npq_c10 <-TukeyHSD(aov(value~read, data = con10split[[2]]))

tk_npq_c1 <-TukeyHSD(aov(value~read, data = con1split[[2]]))

tukeyc10 <- TukeyHSD(aov(value~read*param, data = con10))
tukeyc1 <- TukeyHSD(aov(value~read*param, data = con1))

tk10 <- as.data.frame(tukeyc10[[3]])
tk10$param <- rownames(tk10)
tk10$param <- str_split_fixed(tk10[,5], ":", 3)
tk10$param.1 <- tk10$param[,1]
tk10$param.2 <- tk10$param[,2]
tk10$param.2 <- str_split_fixed(tk10$param.2, "-", 2)
tk10$param.3 <- tk10$param.2[,1]
tk10$param.4 <- tk10$param.2[,2]
tk10$param.5 <- tk10$param[,3]
tk10 <- cbind(tk10[,1:4], tk10$param.1, tk10$param.3, tk10$param.4, tk10$param.5)
tk10 <- subset(tk10, tk10[,6]==tk10[,8] & tk10[,4] < 0.05)


df$param <- factor(df$param, levels=c("Size", "QY_max", "Fv/Fm_Lss","QY_Lss", "NPQ_Lss"))
df$conc = factor(df$conc, levels = c("01 mg/l", "10 mg/l"))
df$read = factor(df$read, levels= c("kont", "ag13", "au13", "au30", "au40", "au50"))

df1 <- ddply(df,~read+conc+param,summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

#######Ploty#########

####Theme####

mytheme <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
  )
#############

myPalette <- c("#95C7FF", "#D070E8", "#FF8048", "#E8DA6A", "#72FF8E", "#7500CC", "#FF3E3E")

facets <- c("QY_max", "Fv/Fm_Lss", "QY_Lss", "NPQ_Lss")

pdf(file="quenching1v2ppm.pdf", width=10, height=4)
p <- ggplot(con1[con1$param %in% facets,], aes(x=read, y=value)) + geom_boxplot(aes(fill = read), size=0.5) + mytheme
p + facet_wrap(~ param, scales="free_y") + scale_fill_manual(limits=c("kontr", "ag_13", "au_13", "au_30", "au_40", "au_50"),values=myPalette)
dev.off()

pdf(file="quenching.pdf", width=10, height=4)
p1 <- ggplot(con10, aes(x=read, y=value)) + geom_boxplot(aes(fill = read), size=0.5) + mytheme
p1 + facet_wrap(~ param, scales="free_y") + scale_fill_manual(limits=c("Kont", "Au13", "Au30", "Au40", "Au50"),values=myPalette) 
dev.off()

myPalette <- c("#95C7FF", "#D070E8", "#FF8048", "#E8DA6A", "#72FF8E", "#7500CC", "#FF3E3E")

#############33

dfp1 <- cbind(con1split$NPQ_Lss[,c(2:5)], con1split$Size[,3])
colnames(dfp1) = c("plant","npq","read","conc" ,"Size")
dfp10 <- cbind(con10split$NPQ_Lss[,c(2:5)], con10split$Size[,3])
colnames(dfp10) = c("plant","npq","read","conc" ,"Size") 

pdf(file="sizevsnpq1ppm.pdf", width=8, height=4)
size <- ggplot(dfp1, aes(x=Size, y=npq)) + geom_point(aes(color=read)) + mytheme
size + facet_wrap(~read, scales="free_y") + stat_smooth(method=lm)
size + scale_colour_manual(values=myPalette) + ggtitle("Concentration 1 ppm")
dev.off()

###Porównanie wg stężenia
concentdf <- split(df, df$param)
sizedf <- as.data.frame(concentdf[1])
pdf(file="sizebyconc.pdf", width=10, height=5)
concent <-  ggplot(sizedf, aes(x=Size.conc, y=Size.value)) + geom_boxplot(aes(fill = Size.read), size=0.5)
concent + facet_wrap(~Size.read, scales="free_y") + scale_fill_manual(values=myPalette)
dev.off()