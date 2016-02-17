
library(plyr)
library(ggplot2)
library(stringr)
library(reshape2)


fl <- list.files(recursive = T, pattern="TXT")

df = data.frame(
  "param" = NA,
  "plant"= NA,
  "value"= NA,
  "genotype"=NA,
  "group"=NA,
  "trt"=NA,
  "al"=NA,
  check.names= FALSE
)

for(i in 1:length(fl)){
  temp <- as.data.frame(read.table(fl[i], skip=3, header=F))
  temp <- rbind(temp[42,],temp[47,],temp[52,], temp[64,])
  temp <- melt(temp, id.vars=("V1"))
  temp$genotype <- (strsplit(fl[i], split = "\\/|_")[[1]])[2]
  temp$group <- (strsplit(fl[i], split ="\\/|_")[[1]])[3]
  temp$trt <- (strsplit(fl[i], split = "\\/|_")[[1]])[4]
  
        if(temp$trt=='1-6.TXT' || temp$trt == '7-12.TXT'){temp$trt="BL"}
        else if(temp$trt=='13-18.TXT' || temp$trt == '19-24.TXT'){temp$trt="RL"}
        else if(temp$trt=='25-30.TXT' || temp$trt == '31-36.TXT'){temp$trt="Control"}
        
  temp$al <- (strsplit(fl[i], split = "\\/|_")[[1]])[1]
  
  colnames(temp)=c("param", "plant", "value", "genotype", "group", "trt", "al")
  df = rbind(df, temp)
}
df = df[2:nrow(df),]



#Test TukeyHSD dla każdego parametru fluorescencji, stężenie 10 ppm:----

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

#####

df$param <- factor(df$param, levels=c("QY_max", "Fv/Fm_Lss","QY_Lss", "NPQ_Lss"))
df$genotype = factor(df$genotype, levels = c("WT", "npq4", "oePsbS"))
df$group = factor(df$group, levels= c("Control", "24h", "48h", "72h"))
df$trt = factor(df$trt, levels = c('Control', 'BL', 'RL'))
df$al = factor(df$al, levels = c("AL200", "AL600"))

df1 <- ddply(df,~genotype+group+param+trt,summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

#######Ploty#########

####Theme####

mytheme <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle=30, vjust =1, hjust = 1)
    #axis.ticks.x = element_blank()
  )
#############

myPalette <- c("#95C7FF", "#D070E8", "#FF8048", "#E8DA6A", "#72FF8E", "#7500CC", "#FF3E3E")

facets <- c("QY_max", "Fv/Fm_Lss", "QY_Lss", "NPQ_Lss")

png(file="quenching.png", width=16, height=10, res=300, units="cm")
p <- ggplot(df, aes(x=group, y=value)) + geom_point(aes(color = genotype), size=1, position = "jitter") + mytheme
p + facet_wrap(~ param, scales="free_y")#+ scale_fill_manual(values=myPalette)
dev.off()

pd <- position_dodge(width = 1)

pdf(file="quenching1.pdf", width=20, height=10)
p <- ggplot(df1, aes(x=group, y=mean, group=genotype)) + geom_line(aes(color = genotype), position = pd, size = 1) + mytheme
p + facet_wrap(trt ~ param,scales= 'free_y') +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=genotype),width=0.3, position = pd, size=1)
dev.off()

pdf(file="quenching_bars.pdf", width=20, height=10)
p <- ggplot(df1, aes(x=group, y=mean)) + geom_bar(aes(fill = genotype),stat="identity",position = position_dodge()) + mytheme
p + facet_wrap(trt ~ param, scales= 'free_y') +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se, color=genotype),position = position_dodge(),width=0.2, size=1)
dev.off()

pdf(file="quenching_box.pdf", width=12, height=8)
p <- ggplot(df, aes(x=group, y=value)) + geom_boxplot(aes(color = genotype), position = position_dodge(1)) + mytheme
p + facet_wrap(trt ~ param, scales= 'free_y', ncol = 4)
dev.off()
###Tukey
df_spl <- split(df,df$param)
tk_npq <- TukeyHSD(aov(value~genotype*group, data = df_spl[[4]]))


###
df[df$group=="control",]$trt <- "Control"

df2 <- subset(df, df$group %in% c('72h'))
df2_spl <- split(df2, df2$param)

tk_npq <- TukeyHSD(aov(value~genotype*trt, data = df2_spl[[4]]))


pdf(file="quenching_box_72_1.pdf", width=12, height=6)
p <- ggplot(df2, aes(x=genotype, y=value)) + geom_boxplot(aes(color = trt)) + mytheme
p + facet_wrap( ~ param, scales = 'free_y')
dev.off()

p <- ggplot(df2, aes(x=group, y=trt)) + geom_tile(aes(fill = value), color = 'white') + mytheme
p + facet_wrap(genotype ~ param) + scale_fill_gradient(low = "white", high = "steelblue")

#rozdzia? wg param
tlabels = c("Control", "Light trt", "24h p. UVC", "48h p. UVC")
######################################
npq <- subset(df, df$param=="NPQ_Lss")
npq <- ddply(npq,~group+genotype+trt+al,
              summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

pd <- position_dodge(width = 0.1)

png(filename = "npq200v600.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(npq, aes(x=group, y=mean, group=al)) + 
        geom_line(aes(color = al), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ genotype) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=al),width=0.3, 
                      position = pd, size=1)
dev.off()#npq al200 v al600
#####################################
qmax <- subset(df, df$param=="QY_max")
qmax <- ddply(qmax,~group+genotype+trt+al,
             summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

pd <- position_dodge(width = 0.1)

png(filename = "qmax200v600.png", width = 15, height = 10, units = "cm", res = 600)
p <- ggplot(qmax, aes(x=group, y=mean, group=al)) + geom_line(aes(color = al), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ genotype) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=al),width=0.3, position = pd, size=1)
dev.off()#qymax al200 v al600
#####################################
qylss <- subset(df, df$param=="QY_Lss")
qylss <- ddply(qylss,~group+genotype+trt+al,
              summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

pd <- position_dodge(width = 0.1)

png(filename = "qylss200v600.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(qylss, aes(x=group, y=mean, group=al)) + geom_line(aes(color = al), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ genotype) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=al),width=0.3, position = pd, size=1)
dev.off()#qylss al200 v al600

####
npq <- subset(df, df$param=="NPQ_Lss")
npq <- ddply(npq,~group+genotype+trt+al,
             summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(2))

pd <- position_dodge(width = 0.1)

png(filename = "npq_genotypes.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(npq, aes(x=group, y=mean, group=genotype)) + 
        geom_line(aes(color = genotype), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=genotype),width=0.3, 
                      position = pd, size=1)
dev.off()

png(filename = "qylss_genotypes.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(qylss, aes(x=group, y=mean, group=genotype)) + 
        geom_line(aes(color = genotype), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=genotype),width=0.3, 
                      position = pd, size=1)
dev.off()

png(filename = "qmax_genotypes.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(qmax, aes(x=group, y=mean, group=genotype)) + 
        geom_line(aes(color = genotype), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(trt ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=genotype),width=0.3, 
                      position = pd, size=1)
dev.off()

####
png(filename = "npq_treatments.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(npq, aes(x=group, y=mean, group=trt)) + 
        geom_line(aes(color = trt), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(genotype ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=trt),width=0.3, 
                      position = pd, size=1)
dev.off()

png(filename = "qylss_treatments.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(qylss, aes(x=group, y=mean, group=trt)) + 
        geom_line(aes(color = trt), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(genotype ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=trt),width=0.3, 
                      position = pd, size=1)
dev.off()

png(filename = "qmax_treatments.png", width = 15, height = 12, units = "cm", res = 600)
p <- ggplot(qmax, aes(x=group, y=mean, group=trt)) + 
        geom_line(aes(color = trt), position = pd, size = 1) + mytheme
p + scale_x_discrete(labels = tlabels) +
        facet_wrap(genotype ~ al, ncol=2) +
        geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color=trt),width=0.3, 
                      position = pd, size=1)
dev.off()
