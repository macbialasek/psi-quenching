library(ggplot2)
library(reshape2)
library(plyr)

mytheme <- theme_bw() +
        theme(
                panel.background = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color="black"),
                axis.title.x = element_blank(),
                #axis.title.y = element_text("Amplitude (mV)"),
                plot.margin = unit(c(0,0,0,0),"mm"),
                #legend.title = element_blank(),
                #legend.position = c(0.2,0.85),
                legend.key.height = unit(0.4, "cm"),
                legend.key.width = unit(0.5, "cm"),
                legend.title = element_blank(),
                strip.text.y = element_blank(),
                strip.text.x = element_blank(),
                legend.position = "none",
                strip.background = element_blank()
        )

##
file_list <- list.files(pattern = "csv")

df = data.frame(
        "yii" = NA,
        "npq" =NA,
        "qp" = NA,
        "time"= NA,
        "read" =NA
)
for(i in 1:length(file_list)){
        
        if(grepl("_1",file_list[i])){
                temp = as.data.frame(read.table(file_list[i], header=F, sep=";",skip=10))[,c(5:7)]
                temp[,4] = seq(from=-60,to=440,by=10)
                temp[,5] = paste("kontr_",file_list[i],sep="")
        }
        if(grepl("_k",file_list[i])) {
                temp = as.data.frame(read.table(file_list[i], header=F, sep=";",skip=10))[,c(5:7)]
                temp[,4] = seq(from=-60,to=440,by=10)
                temp[,5] = paste("lacl3_",file_list[i],sep="")
        }
        
        colnames(temp) = c("yii","npq", "qp", "time", "read")
        
        df = rbind(df, temp)
}

df = df[2:nrow(df),]
df$npq <- df$npq * 4
df$read <- factor(df$read)

df1 <- melt(df, id.vars = c("read", "time"))
baseline <- subset(df1, df$time==0)
bmeans <- ddply(baseline, ~read+variable, summarise, mean=mean(value)) 

byii <-split(bmeans, bmeans$variable)[[1]]
bnpq <- split(bmeans, bmeans$variable)[[2]]
bqp <-split(bmeans, bmeans$variable)[[3]]

peaks <- subset(df1, time > 0 & time < 120)

pyii <- ddply(split(peaks, peaks$variable)[[1]], ~read+variable, summarise,
              peak = min(value))
pnpq <- ddply(split(peaks, peaks$variable)[[2]], ~read+variable, summarise,
              peak = max(value))
pqp <- ddply(split(peaks, peaks$variable)[[3]], ~read+variable, summarise,
              peak = min(value))
              
dfbm <- rbind(byii,bnpq,bqp)
dfpk <- rbind(pyii,pnpq,pqp)

df_final <- merge(dfbm,dfpk)
df_final$amplitude <- df_final$peak-df_final$mean

dfamps <- data.frame(
        read = c(rep("kontr",12),rep("lacl",12)),
        variable = df_final$variable,
        amplitude = df_final$amplitude
)
yiis <- split(lacl_amps,lacl_amps$variable)[[3]]
TukeyHSD(aov(amplitude~read,yiis)) #yii ***
npqs<- split(lacl_amps,lacl_amps$variable)[[1]]
TukeyHSD(aov(amplitude~read,npqs)) # ***
qps <- split(lacl_amps,lacl_amps$variable)[[2]]
TukeyHSD(aov(amplitude~read,qps)) # ***

TukeyHSD(aov(amplitude~read*variable,dfamps)) #***npq, *qp, 

lacl_amps <- read.csv("amplitudy lacl.txt")
lacl_amps$trt <- "lacl"   
kontr_amps <- read.csv("amplitudy mock.txt")
kontr_amps$trt <- "mock"
dfamps <- rbind(lacl_amps,kontr_amps)
dfamps$trt <- factor(dfamps$trt)
#do wykresu
dfamps_gg <- ddply(dfamps, ~ dfamps$variable+dfamps$read+dfamps$trt, summarise, 
                   mean=mean(amplitude), sd=sd(amplitude))
colnames(dfamps_gg) <- c("var","read","trt","mean","sd")
#dfamps_gg[,2] <- factor(dfamps_gg[,1],levels=c("kontr","lacl"))


dfamps_gg <- split(dfamps_gg,dfamps_gg$var)

pd <- position_dodge(width=0.9)

png(filename = "npq1.png", width = 8, height = 5, res = 600, units = "cm")
ggplot(dfamps_gg[[1]], aes(x=dfamps_gg[[1]]$trt, y= mean,fill= dfamps_gg[[1]]$read)) + 
        geom_bar(color="black",  stat = "identity", position=pd) + 
        geom_errorbar(aes(ymin=dfamps_gg[[1]]$mean-dfamps_gg[[1]]$sd, 
                          ymax=dfamps_gg[[1]]$mean+dfamps_gg[[1]]$sd), width=0.2, position=pd) + 
        scale_x_discrete(labels=c(expression(LaCl[3]), "Mock")) +
        scale_fill_manual(labels=c("Control", "Treatment"),values=c("black","white")) +
        mytheme + ylab("Amplitude of NPQ")
dev.off()

png(filename = "yii.png", width = 8, height = 5, res = 600, units = "cm")
ggplot(dfamps_gg[[3]], aes(x=dfamps_gg[[3]]$trt, y= mean,fill= dfamps_gg[[3]]$read)) + 
        geom_bar(color="black",  stat = "identity", position=pd) + 
        geom_errorbar(aes(ymin=dfamps_gg[[3]]$mean-dfamps_gg[[3]]$sd, 
                          ymax=dfamps_gg[[3]]$mean+dfamps_gg[[3]]$sd), width=0.2, position=pd) + 
        scale_x_discrete(labels=c(expression(LaCl[3]), "Mock")) +
        scale_fill_manual(labels=c("Control", "Treatment"),values=c("black","white")) +
        mytheme + ylab("Amplitude of Y(II)")
dev.off()

png(filename = "qp.png", width = 8, height = 5, res = 600, units = "cm")
ggplot(dfamps_gg[[2]], aes(x=dfamps_gg[[2]]$trt, y= mean,fill= dfamps_gg[[2]]$read)) + 
        geom_bar(color="black",  stat = "identity", position=pd) + 
        geom_errorbar(aes(ymin=dfamps_gg[[2]]$mean-dfamps_gg[[2]]$sd, 
                          ymax=dfamps_gg[[2]]$mean+dfamps_gg[[2]]$sd), width=0.2, position=pd) + 
        scale_x_discrete(labels=c(expression(LaCl[3]), "Mock")) +
        scale_fill_manual(labels=c("Control", "Treatment"),values=c("black","white")) +
        mytheme + ylab("Amplitude of qP")
dev.off()
