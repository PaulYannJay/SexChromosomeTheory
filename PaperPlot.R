library(tidyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gghalves)
library(plyr)
library(colorRamps)
library(cmocean)
library(viridis)
library(ggnewscale)
library(tidyverse)
library(directlabels)
options(scipen=999)

#### Fig 2 ###
#### Fig 2_XY ###

### Deterministe ##
### Time simulation , Fig 2a_XY
###

Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)
Data_Y=read.table("~/Analysis/DelSheltering/Output/TimeSimul_h_s_u1e-8_cM_PapData_2MbInv_0.80_Fl_XYsyst_Y_BC.txt", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1) #Focus on two s value
Data_Y=Data_Y[Data_Y$s %in% FocS,]
Data_Y$r=Data_Y$r*100 # Change recombination rate to distance in centimorgan
FocR=c(0,50)
Data_Y=Data_Y[Data_Y$r %in% FocR,] #Focus on fully linked or fully unlinked inversion
Data_Y$Chrom="Y-linked"
Data_Y[Data_Y$r==50,]$Chrom="Autosomal" #Full unlinked inversion are considered to be on autosome

Data=Data_Y
Data$Freq=Data$FYI+Data$FXI #Overall frequency of the inversion

base=ggplot(Data)
PlotA=base+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  geom_line(aes(x=time, y=Freq, linetype=as.factor(h), color=Chrom), size=1, alpha=0.6)+
  facet_grid(paste0("s=", s)~., switch = "y")+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_linetype_manual("h=",values=c("solid", "dashed", "dotted"))+
  scale_color_manual("Chromosome", values=Col, guide=F)+
  theme(panel.border = element_blank(),  
        legend.position = c(0.6,0.98),
        legend.direction = "horizontal",
        legend.key.width = unit(1.2,"cm"),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=20),
        axis.line = element_line(colour = "grey"),
  )+
  scale_y_continuous(breaks=c(0.0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.05,1.1))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))+
  geom_dl(aes(x=time-3000, y=Freq+0.03, label=Chrom, color=Chrom),method=list('last.bumpup', cex = 1.2, hjust = -0.1))
PlotA
save_plot("Fig2a.png", PlotA, nrow=2)

### Fig 2B ###

Simul=read.table(paste("LinkedAndUnlinked_2MbInv_IntroduceInvFromInit_Nmut_Freq_Fit_IndivPlot_MidSDR_XY.txt",sep=""), stringsAsFactors = F)
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv","Chrom", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Simul=Simul[!(Simul$DebInv>10000000 & Simul$Chrom=="Y"),] # 20000 inversion are run in autosome. Removing here all the inversion that were introduce in Y-bearing slim genome but totally unlinked to the Y allele, to get 10000 inversion
Simul[Simul$DebInv>10000000,]$Chrom="Autosome" #Inversion on chromosome 2 are unlinked to the Sex-determining loci, therefore in an autosome
Simul=Simul[!Simul$Chrom=="X",] #For the main plot, do not consider inversion linked to the X allele
Simul[Simul$Chrom=="Y",]$Freq=Simul[Simul$Chrom=="Y",]$Freq*4 #Multiply Y-inversion frequency by four to get the overall frequency
Simul$InvSize=Simul$FinInv - Simul$DebInv  #Inversion size
SimulSub=Simul
SimulSub=SimulSub[!(SimulSub$Freq==0.0),] #Remove generation were inversion fall at 0 frequency
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv,SimulSub$Chrom, SimulSub$Rep, sep="_") #Define a code for every simulation
FocS=c(-0.01,-0.1) #Focus on two s value
SimulSub=SimulSub[SimulSub$s %in% FocS,] 
SimulSub20G=unique(SimulSub[SimulSub$Gen>15020,]$Code) #Consider only inversion not lost after 20 generation.
SimulSub=SimulSub[SimulSub$Code %in% SimulSub20G,] 
summarySub=ddply(SimulSub, .(Code), summarize, LastFreq=last(Freq), maxGen=max(Gen)) #Summarize what happened to inversion maximum frequency or the generation when they were lost or fixed
#Since we ended simulation when the inversions were fixed of lost (this allows to save much time), for plotting purpose we need to recreate the end of these simulation --> in both case, the frequency of these fixed or lost simulationmust have been the same for the rest of the simulations, being either 1 or 0
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Code of Inversion that have been lost or fixed (since we ended simulation when this occured)
GoodSimulSub=SimulSub #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15010),]#For lost or fixed inversion, grep their initial state (for instance a generation 10)
FalseEndGoodSimulSub$Gen=25000 #Create a false simulation end frequency
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
if (length(summarySub[(summarySub$LastFreq>0.95 & summarySub$maxGen != 25000) ,]$Code)>0) #If their last recorded frequency is above 0.95, they have fixed
{
  FixedSimul=summarySub[summarySub$LastFreq>0.95,]$Code #Grep the code of the inversion that have fixed
  FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
}
GoodSimulSub=rbind(GoodSimulSub,FalseEndGoodSimulSub) #Concatenate the dataset (by including the false end)
DataAll=GoodSimulSub
DataAll$s=-DataAll$s #reverse the s parameter as it is defined differently in SLiM than in the deterministic analyses

DataAllEnd=subset(DataAll, DataAll$Gen==25000) #Summarise what happens at the end
DataAllEnd$state="Segregate"
DataAllEnd[DataAllEnd$Freq==0.0,]$state="Lost"
DataAllEnd[DataAllEnd$Freq==1.0,]$state="Fixed"
DataAll$Gen=DataAll$Gen-15000 #Since the simulation started at generation 15000 (after initialisation), correct this for plotting

SumEnd=DataAllEnd %>% count(h,s,Chrom,state, sort = TRUE) #Count number of simulation that ended up in each state
SumEnd$Pos=0.0 #Plotting position for these count ("n=""")
SumEnd[SumEnd$state=="Segregate",]$Pos=0.5
SumEnd[SumEnd$state=="Fixed",]$Pos=1.0

Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)

base=ggplot(DataAll[DataAll$h==0.01,]) #Focus on case with h=0.01
PlotB=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chrom)), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("InversionFrequency")+
  geom_text(data=SumEnd[SumEnd$h==0.01,], aes(x=8000, y=Pos, label=paste0("n=",n), color=as.factor(Chrom)), vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=20),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2, alpha=1.0)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(s ~ Chrom)
PlotB

save_plot("Fig2b.png", PlotB, ncol=2, nrow=2)

###

### Fig 2c_XY Proba Spread ###
Simul=read.table(paste("AllInv_LinkedAndUnlinked_IntroduceInvFromInit_Nmut_Freq_Fit_IndivPlot_MidSDR_XY.txt",sep=""), stringsAsFactors = F)
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv","Chrom", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")

#We do the same thing as for Figure 2b
Simul=Simul[!(Simul$DebInv>10000000 & Simul$Chrom=="Y"),] # 20000 inversion are run in autosome. Removing here all the inversion that were introduce in Y-bearing slim genome but totally unlinked to the Y allele, to get 10000 inversion
Simul[Simul$DebInv>10000000,]$Chrom="Autosome" #Inversion not linked to the Y allele are in autosome (Chromosome 2)
Simul=Simul[!Simul$Chrom=="X",] #Not consider X-linked inversion
Simul$InvSize=Simul$FinInv - Simul$DebInv #inversion size
SimulSub=Simul
SimulSub$cM=SimulSub$r * (SimulSub$DebInv-1) * 100 # Change the unity of the recombination (cM)
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv,SimulSub$Chrom, SimulSub$Rep, sep="_") #Define a Code for each simulation
FocS=c(-0.01,-0.1) #Focus on two selection coefficiennt

SimulSub=SimulSub[SimulSub$s %in% FocS,] # Focus only on linked or unlinked inversion
summarySub=ddply(SimulSub, .(Code), summarize, maxFreq=max(Freq), maxGen=max(Gen)) # For each simulation, grep its max frequency and it end generation (which can be differet from 25000 in case the inversion is lost or fixed)
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Inversion that have been lost or fixed
NonLostSimulSub=summarySub[summarySub$maxGen==25000,]$Code #Inversion still segregating at the end
GoodSimulSub=SimulSub[SimulSub$Code %in% NonLostSimulSub,] #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15010),]#For lost inversion, grep their initial state
FixedSimul=summarySub[summarySub$maxFreq>0.90,]$Code #Grep the code of the inversion that have fixed
FalseEndGoodSimulSub$Gen=25000 #Modify their initial state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
GoodSimulSubEnd=rbind(GoodSimulSub[GoodSimulSub$Gen==25000,],FalseEndGoodSimulSub) #Concatenate the dataset
DataAllEnd=GoodSimulSubEnd
DataAllEnd2=DataAllEnd # Save, just un case...
DataAllEnd2$State=0 ##Lost inversion
DataAllEnd2[DataAllEnd2$Freq>0,]$State=1 #Fixed or segregating inversion have their state=1
DataSummary=ddply(DataAllEnd2, .(u,r,h,s,DebInv,FinInv,Chrom), summarize, ProbSpread=mean(State)) #ProbSpread is the fraction of inversion fixed or segregating
DataSummary$InvSize=DataSummary$FinInv - DataSummary$DebInv #Inversion size
FocSize=c(1000000,2000000,5000000) 
DataSummary=DataSummary[DataSummary$InvSize %in% FocSize,] #Focus on three inversion size

Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)
base=ggplot(DataSummary)
PlotC=base+geom_point(aes(x=h,y=ProbSpread, color=Chrom,shape=as.factor(InvSize)),
                      size=5, alpha=0.8)+
  geom_line(aes(x=h, y=ProbSpread, color=Chrom,linetype=as.factor(InvSize)), size=1.5, alpha=0.6)+
  geom_dl(aes(x=h, y=ProbSpread, label=Chrom, color=Chrom),method=list('last.bumpup', cex = 1.5,vjust=0.1, hjust = -0.2))+
  scale_x_log10(breaks=c(0.5, 0.1,0.01, 0.001), expand = expansion(add = c(0.1, 0.8)))+
  ylab("Fraction of inversions spreading or fixed")+
  xlab("Dominance (h)")+
  scale_color_manual("", values=Col, guide=F)+
  scale_linetype("Inversion size:", labels=c("1Mb", "2Mb", "5Mb"), guide=F)+
  scale_shape_manual("Inversion size:", values=c(15,16,17,18), labels=c("1Mb", "2Mb", "5Mb"))+
  facet_grid(fct_relevel(as.factor(s), "-0.01", "-0.1")~., switch="y")+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    legend.position = c(0.50,0.98),
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(size=20),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.text = element_blank())+
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
PlotC

save_plot("Fig2c.png", PlotC, nrow=2, base_aspect_ratio = 2)

### Figure 3 
###Heatmap Fig3b
Data=read.table("N=1000_u=1e-08_r=1e-05_MaxSizeInv=20000000_Rep_9_Nrecomb_IndivSimulation_XY.Parsed.txt", stringsAsFactors = F, header = F)
Pos=seq(1:(length(colnames(Data))-1))
colnames(Data)=c("Generation", Pos)
DataLong=gather(Data, Position, Value, "1":"200")
DataLong$Position=as.numeric(as.character(DataLong$Position))
DataLong$Chrom="Chromosome 1 (X/Y)"
DataLong[DataLong$Position>100,]$Chrom="Chromosome 2"
DataLong$NormVal=DataLong$Value/mean(DataLong[DataLong$Generation<15000,]$Value) #Normalize recombination rate
DataLong[DataLong$Chrom=="Chromosome 2",]$Position=DataLong[DataLong$Chrom=="Chromosome 2",]$Position - 100 #The position of Chromosome 2  need to be redifined
base=ggplot(DataLong[DataLong$Generation < 101000,])
PlotRecomb=base+geom_tile(aes(x=Position, y=Generation, fill=NormVal),colour="white",size=0.02)+
  scale_y_reverse(expand=c(0.005,0),breaks=pretty(DataLong[DataLong$Generation < 101000,]$Generation, 10))+
  geom_vline(data=DataSub, aes(xintercept = 50, color="Sex-determining locus"))+
  scale_color_manual("", values=c("red"))+
  scale_x_continuous(expand=c(0.01,0.01))+
  scale_fill_viridis("Recombination rate", option = "plasma", direction = -1, limits=c(0,max(DataLong$NormVal)))+
  facet_wrap(.~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x")+
  xlab("Position (Mb)")+
  theme(panel.border = element_blank(),  
        panel.grid.major = element_line(color="grey", size=0.1),
        panel.grid.minor = element_line(color="grey", size=0.08),
        panel.background = element_blank(),
        strip.text = element_blank())
PlotRecomb
### InvEvol ## Fig3a
DataInv=read.table("N=1000_u=1e-08_r=1e-05_MaxSizeInv=20000000_Rep_9_InvFreq_IndivSimulation_XY.Parsed.txt", stringsAsFactors = F, header = F)
colnames(DataInv)=c("Generation", "InvStart", "InvEnd", "p__", "pI_", "pII", "Yfreq", "Xfreq")
DataInv$InvStart=DataInv$InvStart/1000000 #Position in Mb
DataInv$InvEnd=DataInv$InvEnd/1000000
DataInv$Yfreq=DataInv$Yfreq/500 #Frequency of Y inversions (their is 500 proto Y chromosome)
DataInv$Xfreq=DataInv$Xfreq/1500 #Frequency of X inversions (their is 1500 proto X chromosome)
DataInv$Chrom="Chromosome 1 (X/Y)"
DataInv[DataInv$InvStart>100,]$Chrom="Chromosome 2"
DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart=DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart - 100
DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd=DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd - 100
Col=scales::viridis_pal(begin=0.20, end=0.80, option="cividis")(2)
DataInvSub=subset(DataInv, (DataInv$Chrom=="Chromosome 1 (X/Y)" & DataInv$Generation %% 10000 == 0 & DataInv$Generation <101000))
base=ggplot(DataInv[(DataInv$Generation %% 10000 == 0 & DataInv$Generation <101000), ])
PlotFreq=base+geom_rect(aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Yfreq, fill="Y"),color="black", alpha=0.6, size=0.2)+
  geom_rect(aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Xfreq, fill="X"),color="black", alpha=0.6, size=0.2)+
  scale_fill_manual("Inversion on chromosome:", values = Col)+
  ylab("Inversion Frequency")+
  scale_color_manual("", values=c("red"))+
  geom_vline(data=DataInvSub, aes(xintercept = 50, color="Sex-determining locus"))+
  scale_x_continuous(expand=c(0.01,0.01))+
  scale_y_continuous(expand=c(0,0),breaks = c(0.0,0.5,1.0), position = "right")+
  facet_grid(Generation~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x", drop=F, switch="y")+
  theme(#strip.background.y = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    panel.border = element_blank(),  
    panel.grid.major = element_line(color="grey", size=0.1),
    panel.grid.minor = element_line(color="grey", size=0.08),
    panel.background = element_blank(),
    axis.text.y.right = element_text(size=7))
PlotFreq

plots <- align_plots(PlotFreq, PlotRecomb, align = 'v', axis = 'lr')

MergedPlot=plot_grid(plots[[1]],plots[[2]],  ncol=1, labels = c('a', 'b'))


save_plot("Fig3.png", MergedPlot, ncol = 2, nrow=2)



