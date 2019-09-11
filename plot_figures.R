#-------------------------------------------------------------
# Code to generate figures for clonal nephrogenesis paper
# Tim Coorens - September 2019
#-------------------------------------------------------------
# This script includes most of the plotting code for the paper (or at least a way to get to the data)
# Any other code or data can be requested by contacting us

#-------------------------------------------------------------
# Libraries and data
#-------------------------------------------------------------

options(stringsAsFactors = F)
library(openxlsx)
library(beeswarm)
library(viridis)

sample_info=read.xlsx("Supplementary_tables_S1-S8_R2.xlsx",sheet=1,startRow=4,colNames=T) #Table S1
normal_muts=read.xlsx("Supplementary_tables_S1-S8_R2.xlsx",sheet=2,startRow=4,colNames=T) #Table S2
tumour_mutations=read.xlsx("Supplementary_tables_S1-S8_R2.xlsx",sheet=4,startRow=4,colNames=T) #Table S4
methylation_data=read.xlsx("Supplementary_tables_S1-S8_R2.xlsx",sheet=6,startRow=4,colNames=T) #Table S6
kidney_muts=read.xlsx("Supplementary_tables_S1-S8_R2.xlsx",sheet=8,startRow=4,colNames=T) #Table S8

#-------------------------------------------------------------
# Figure 1C
#-------------------------------------------------------------
#Get the VAFs of the first five mutations in Fig 1C
embryonic=normal_muts[normal_muts$Patient=="PD37272"&normal_muts$Tumour&normal_muts$PD37272g!="NS"&normal_muts$PD37272d!="NS",grepl("PD37272",colnames(normal_muts))]
nephrogenic=normal_muts[normal_muts$Patient=="PD37272"&normal_muts$Tumour&normal_muts$Type=="Nephrogenic",grepl("PD37272",colnames(normal_muts))]
nephrogenic[nephrogenic=="NS"]=0
muts=rbind(embryonic,nephrogenic)
muts=apply(muts,2,as.numeric)
muts=muts[order(rowMeans(muts),decreasing = T),]

#-------------------------------------------------------------
# Figure 2 
#-------------------------------------------------------------

rownames(sample_info)=sample_info$Sample
normal_samples_wilms=sample_info$Sample[grepl("kidney",sample_info$Sample_type_2)&sample_info$Cohort%in%c("Discovery","Extension")]
clonal_nephrogenesis_samples=normal_samples_wilms[ #Get samples with clonal nephrogenesis
  unlist(lapply(normal_samples_wilms,function (x) !all(normal_muts[normal_muts$Type=="Nephrogenic"& #Select the nephrogenic mutations
                                                                     normal_muts$Tumour,x]%in%c("NS","-"))))] #Select the samples that have a significant VAF for them


#Get clone sizes for samples with clonal nephrogenesis
clone_sizes_df=data.frame(Sample=clonal_nephrogenesis_samples,
                          ID=NA,
                          Est=NA,
                          Min=NA,
                          Max=NA)
for (k in 1:nrow(clone_sizes_df)){
  sample=clone_sizes_df$Sample[k]
  muts=kidney_muts[kidney_muts$Sample==sample&kidney_muts$Tumour,c("ID","VAFs","NV","NR")] 
  mut=muts[which.max(muts$VAFs),]
  test = binom.test(x=mut$NV,n=mut$NR) #Use binom.test to get confidence intervals
  clone_sizes_df$ID[k]=mut$ID
  clone_sizes_df[k,c("Min","Max")] = test$conf.int[1:2]
  clone_sizes_df$Est[k] = test$estimate
}
clone_sizes_df$PD=substr(clone_sizes_df$Sample,1,7)
# Add methylation values to clone sizes
clone_sizes_df$H19_methylation=sample_info[clone_sizes_df$Sample,"H19.Methylation.(beta.score)"]
clone_sizes_df$KvDMR1_methylation=sample_info[clone_sizes_df$Sample,"KvDMR1.Methylation.(beta.score)"]

#get background mean methylation levels from normal kidney for plotting purposes
#Note: statistical testing is done using whole distribution, rather than just the mean
background_means_H19=sample_info$'H19.Methylation.(beta.score)'[grepl("kidney",sample_info$Sample_type_2)&
                                                        !is.na(sample_info$'H19.Methylation.(beta.score)')&
                                                        !sample_info$Sample%in%clone_sizes_df$Sample]

# Which normal kidney samples (with clonal nephrogenesis) are hypermethylated compared to normal background?
normal_background_samples=sample_info$Sample[grepl("kidney",sample_info$Sample_type_2)&
                                               !is.na(sample_info$'H19.Methylation.(beta.score)')&
                                               !sample_info$Sample%in%clone_sizes_df$Sample]

#Sample with clonal nephrogenesis and methylation data
CN_methylation_samples=colnames(methylation_data)[colnames(methylation_data)%in%clone_sizes_df$Sample] 
#All probes in H19 belonging to normal, polyclonal renal tissues
H19_background_full=methylation_data[methylation_data$Region=="H19",colnames(methylation_data)%in%normal_background_samples] 
#Same, but KvDMR1
KvDMR1_background_full=methylation_data[methylation_data$Region=="KvDMR1",colnames(methylation_data)%in%normal_background_samples] #Same, but KvDMR1

H19_qval=KvDMR1_qval=rep(1,length(CN_methylation_samples))
for (n in 1:length(CN_methylation_samples)){
  s=CN_methylation_samples[n]
  #Use Wilcoxon rank sum test on entire distribution
  H19_qval[n]=p.adjust(wilcox.test(methylation_data[methylation_data$Region=="H19",s],
                                   unlist(H19_background_full))$p.value,n=length(CN_methylation_samples))
  KvDMR1_qval[n]=p.adjust(wilcox.test(methylation_data[methylation_data$Region=="KvDMR1",s],
                                unlist(KvDMR1_background_full))$p.value,n=length(CN_methylation_samples))
  print(paste0("H19 in sample ",s,": ",H19_qval[n]))
  print(paste0("KvDMR1 in sample ",s,": ",KvDMR1_qval[n]))
}

signif_H19=CN_methylation_samples[H19_qval<0.05] 
#Should be 7: "PD40729c" "PD40735h" "PD40735b" "PD40738g" "PD40727c" "PD36159g" "PD36159h"
signif_KvDMR1=CN_methylation_samples[KvDMR1_qval<0.05] 
#Should be none

signif=clone_sizes_df$Sample%in%signif_H19

#----------------------
#Figure 2A
#----------------------
dens_H19=density(background_means_H19,n=250)
cols = unlist(lapply(dens_H19$y/max(dens_H19$y)*0.25+0.02,function(x) adjustcolor("dodgerblue3",alpha.f=x)))

plot(1,xlim=c(-60,215),ylim=c(-2,nrow(clone_sizes_df)+2),
     type = "n",xlab="", ylab="",yaxt="n",xaxt='n',bty="n")
freq=table(clone_sizes_df$PD)
start=nrow(clone_sizes_df)+0.5
for(n in 1:(length(freq)-1)){
  start=start-freq[n]
  segments(x0=-5,x1=200,y0=start,col='grey60')
}
rect(ybottom=0.5,xleft=-5,xright=105,ytop=nrow(clone_sizes_df)+0.5)
points(x=200*clone_sizes_df$Est,y=nrow(clone_sizes_df):1,pch=19)
segments(x0=200*clone_sizes_df$Min,
         x1=200*pmin(0.5,clone_sizes_df$Max),
         y0=nrow(clone_sizes_df):1)
segments(y0=seq(nrow(clone_sizes_df)+0.2,1.2,by=-1),y1=seq(nrow(clone_sizes_df)-0.2,0.8,by=-1),x0=200*clone_sizes_df$Min)
segments(y0=seq(nrow(clone_sizes_df)+0.2,1.2,by=-1),y1=seq(nrow(clone_sizes_df)-0.2,0.8,by=-1),x0=200*pmin(0.5,clone_sizes_df$Max))
segments(y0=0.5,y1=0.1,x0=c(0,50,100))
segments(y0=0.5,y1=0.25,x0=c(25,75))
text(y=-0.25,x=seq(0,100,25)+2.5,labels=paste0(seq(0,100,25),"%"),pos=1,srt=-90)
segments(y0=0.5,y1=0.25,x0=c(25,75))
text(x=-5.5,y=nrow(clone_sizes_df):1,labels=clone_sizes_df$Sample,pos=2)

rect(ybottom=0.5,xleft=105,xright=200,ytop=nrow(clone_sizes_df)+0.5)
segments(x0=105+100*((dens_H19$x-0.50)/0.35),y0=nrow(clone_sizes_df)+0.5,y1=0.5,col=cols)
colvec_signif=rep('grey50',nrow(clone_sizes_df))
colvec_signif[signif]="olivedrab1"
points(x=105+100*(clone_sizes_df$H19_methylation-0.5)/0.35,y=nrow(clone_sizes_df):1,pch=21, bg=colvec_signif)
segments(y0=0.5,y1=0.1,x0=105+100*(c(0.6,0.75)-0.5)/0.35)
text(y=-0.25,x=107.5+100*(c(0.6,0.75)-0.5)/0.35,labels=c(0.60,0.75),pos=1,srt=-90)

text(y=-0.25,x=212.5+100*(c(0.75,0.85)-0.7)/0.2,labels=c(0.75,0.85),pos=1,srt=-90)
text(x=50,y=nrow(clone_sizes_df)+1,label="Clone sizes")
text(x=150,y=nrow(clone_sizes_df)+1,label="H19 methylation")

#----------------------
#Figure 2B
#----------------------

kidney_muts$Group=factor(kidney_muts$Group,levels = c("Clonal Nephrogenesis - Tumour",
                                                              "Clonal Nephrogenesis - Non-Tumour",
                                                              "Wilms - Non-Clonal Nephrogenesis",
                                                              "CMN/MRT",
                                                              "RCC",
                                                              "Transplant Donor"))
colvec=c("firebrick","dodgerblue3",rep("grey40",4))
light_colvec=adjustcolor(colvec,alpha.f = 0.4)
beeswarm(2*VAFs ~ Group, data = kidney_muts, method='hex',xlim = c(0,17), 
         ylim = c(0,1),ylab="Contribution to bulk tissue (VAF x 2)",xlab="",
         corral='gutter', pch = 21, col = colvec,bg=light_colvec,
         corralWidth = 2,cex=1.3,at =seq(from=1,by=3,length.out = 6),las=2)

#test for difference in VAF between clonal nephrogenic mutations and all other muts from non-clonal nephrogenic samples
wilcox.test(kidney_muts$VAFs[kidney_muts$Group=="Clonal Nephrogenesis - Tumour"],
            kidney_muts$VAFs[!grepl("Clonal",kidney_muts$Group)])
#Should be p-value = 2.548e-05

#----------------------
#Figure 2D
#----------------------
library(viridis)
y_start=2
plot(1,xlim=c(-0.1,1.1),ylim=c(-0.5,3.2),
     type = "n",xlab="", ylab="",yaxt="n",xaxt='n',bty="n")

segments(x0=c(0.25,0.5,0.75),y0=0,y1=3,col='grey80')
segments(x0=0,x1=1,y0=1:2,lwd=1.5)
rect(ybottom=0,xleft=0,xright=1,ytop=y_start+1,lwd=1.5)

cols=viridis(51)
names(cols)=seq(0,0.5,0.01)
for(p in c(0.5,0.25,0.1)){
  covg=rpois(n=1000,lambda = 30) #simulate coverage
  NV=unlist(lapply(covg,function(x) rbinom(n=1,size=x,prob=as.numeric(p)))) #simulate number of reads supporting variant
  dens=density(NV[NV>3]/covg[NV>3]) #Truncate binomial to reflect minimum number of reads needed by CaVEMan
  lines(x=dens$x,y=y_start+0.9*dens$y/max(dens$y),xlim=c(0,1),xlab="VAF",lwd=2)
  polygon(x=dens$x,y=y_start+0.9*dens$y/max(dens$y), col=adjustcolor(cols[as.character(p)],alpha.f = 0.65), border="black")
  y_start=y_start-1 
}
segments(x0=c(0,0.5,1),y1=0,y0=-0.2)
segments(x0=c(0,0.25,0.75),y1=0,y0=-0.15)
text(y=-0.3,x=seq(0,1,0.25),labels=seq(0,1,0.25),pos=1,cex=0.8)

#----------------------
#Figure 2E
#----------------------
#Example
tubules=c("PD36720b_CA_9__F1","PD28690hk_CA_9_A3","PD28690hk_CA_9_B3","PD36723c_CA_9_B4","PD36723b_CA_9_C4")
tubules=c(tubules,read.table("Normal_Kidney_LCM/samples.txt")[,1])

library(viridis)
cols=viridis(51)
names(cols)=seq(0,0.5,0.01)

y_start=length(tubules)-1
plot(1,xlim=c(-0.75,1.5),ylim=c(-1,length(tubules)+1),
     type = "n",xlab="", ylab="",yaxt="n",xaxt='n',bty="n")

rect(ybottom=0,xleft=0,xright=1,ytop=y_start+1,lwd=1.5)
segments(x0=seq(0,1,0.25),y0=0,y1=length(tubules),col='grey80')

segments(x0=0,x1=1,y0=1:(length(tubules)-1),lwd=1.5)

for(sample in tubules){
  muts=read.table("LCM_example_muts.txt",header=T)
  res=readRDS("LCM_example_binom_mix.Rdata") #result from truncated binomial mixture model
  dens=density(muts$VAF[muts$Chr%in%c(1:22)])
  lines(x=dens$x,y=y_start+0.8*dens$y/max(dens$y),xlim=c(0,1),xlab="VAF",lwd=2)
  colname=as.character(round(max(res$p),digits=2))
  if(max(res$p)>0.5)colname=as.character(round(min(res$p),digits=2))
  polygon(x=dens$x,y=y_start+0.8*dens$y/max(dens$y), col=adjustcolor(cols[colname],alpha.f = 0.65), border="black")
  
  y_start=y_start-1 
}
segments(x0=c(0,0.5,1),y1=0,y0=-0.2)
segments(x0=c(0,0.25,0.75),y1=0,y0=-0.15)
text(y=-0.3,x=seq(0,1,0.25),labels=seq(0,1,0.25),pos=1,cex=0.8)
text(y=(length(tubules)-0.5):0.5,x=-0.02,labels=tubules,pos=2,cex=0.8)
# 
horiz_start=seq(8,length(tubules)-8,length.out = 52)
for (n in 1:51){
  rect(xleft=1.175,xright=1.325,ybottom=horiz_start[n],ytop=horiz_start[n+1],col=adjustcolor(cols[n],alpha.f = 0.65),border=F)
}
rect(xleft=1.175,xright=1.325,ytop=length(tubules)-8,ybottom=8)

text(labels="Fitted peak of\nautosomal variant\nallele frequency",x=1.25,y=length(tubules)-7.25,cex=0.8)
text(labels="0",x=1.4,y=8,cex=0.8)
text(labels="0.25",x=1.4,y=length(tubules)/2,cex=0.8)
text(labels="0.5",x=1.4,y=length(tubules)-8,cex=0.8)
text(labels="Variant Allele Frequency",x=0.5,y=-0.75,cex=0.8,pos=1)

#----------------------
#Figure 2F
#----------------------
H19_methylation=data.frame(Sample=sample_info$Sample[!is.na(sample_info$H19.Methylation..beta.score.)],
                           Meth=sample_info$H19.Methylation..beta.score.[!is.na(sample_info$H19.Methylation..beta.score.)],
                           Sample.type=sample_info$Sample_type_2[!is.na(sample_info$H19.Methylation..beta.score.)],
                           Group="Tumours")

H19_methylation$Group[H19_methylation$Sample.type=="Blood"]="Blood"
H19_methylation$Group[grepl("kidney",H19_methylation$Sample.type)&H19_methylation$Sample%in%clone_sizes_df$Sample]="Clonal Nephrogenesis"
H19_methylation$Group[grepl("kidney",H19_methylation$Sample.type)&!H19_methylation$Sample%in%clone_sizes_df$Sample]="Normal Kidney"
H19_methylation$Group=factor(H19_methylation$Group,levels=c("Blood","Normal Kidney","Clonal Nephrogenesis","Tumours"),ordered=T)

boxplot(data=H19_methylation,as.numeric(Meth)~Group,ylab="Methylation (beta)",main="Methylation of H19",
        col=c(rgb(60,60,60,maxColorValue = 255),
              rgb(197,197,197,maxColorValue = 255),
              rgb(89,158,224,maxColorValue = 255),
              rgb(173,58,58,maxColorValue = 255)),cex=0.9,las=2)

wilcox.test(H19_methylation$Meth[H19_methylation$Group=="Clonal Nephrogenesis"],
            H19_methylation$Meth[H19_methylation$Group=="Normal Kidney"])

#Should be p-value = 0.02834

#Figure 2h

pdf("Methylation_KvDMR1_boxplot.pdf",useDingbats = F)
boxplot(data=df_sub,as.numeric(Meth)~Type,ylab="Methylation (beta)",main="Methylation of KvDMR1",
        col=c(rgb(60,60,60,maxColorValue = 255),
              rgb(197,197,197,maxColorValue = 255),
              rgb(89,158,224,maxColorValue = 255),
              rgb(173,58,58,maxColorValue = 255)),cex=0.9)
dev.off()
#----------------------
#Figure 2G
#----------------------
newx=seq(-0.2,1.2,0.01)
clone_sizes_df$Clone_size=pmin(1,clone_sizes_df$Est*2)
clone_sizes_meth=clone_sizes_df[!is.na(clone_sizes_df$H19_methylation),]
clone_size_meth_somatic=clone_sizes_meth[clone_sizes_meth$Sample!="PD40738g",]

corr=lm(data=clone_size_meth_somatic,H19_methylation~Clone_size)
conf_interval <- predict(corr, newdata=data.frame(Clone_size=newx), interval="confidence",
                         level = 0.95)
plot(x=clone_sizes_meth$Clone_size,
     y=clone_sizes_meth$H19_methylation,xlim=c(0,1),col='white',bg="white",
     xlab="Proportion of\nNephrogenic Clone",
     ylab="Methylation score (Beta)",
     main="Methylation of H19 in samples with clonal nephrogenesis")
polygon(x=c(newx,rev(newx)),y=c(conf_interval[,3],rev(conf_interval[,2])),col=adjustcolor("lightblue",alpha.f = 0.35),border=F)
abline(corr,col="firebrick",lwd=2,lty="dashed")

colour=rep("dodgerblue",nrow(clone_sizes_df))
colour[clone_sizes_meth$Sample=="PD40738g"]="black"
points(x=clone_sizes_meth$Clone_size,
       y=clone_sizes_meth$H19_methylation,
       col=colour,pch=19,cex=2.5)
text(labels=paste0("r=",round(cor(clone_size_meth_somatic$H19_methylation,clone_size_meth_somatic$Clone_size),digits=2)),
     y=0.57,x=0.9) #Should be r=0.86

#Figure 2j
LOH11p=cbind(c(1,5),c(17,4))
colnames(LOH11p)=c("LOH","WT")
rownames(LOH11p)=c("With clonal nephrogenesis","Without clonal nephrogenesis")

barplot(t(LOH11p),horiz = T,las=2,col=c("black","white"))
fisher.test(rbind(c(17,1),c(4,5)))

#------------------------------------
#Figure 3A (Trees - branch length)
#------------------------------------
# This bit of code will generate the numbers used for the trees

tumour_mutations$ID=paste(tumour_mutations$Chr,tumour_mutations$Pos,tumour_mutations$Reference_base,tumour_mutations$Mutation,sep='_')
#make tree for patient PD40735 (Bilateral)

patient="PD40735"

tumour_sub=tumour_mutations[tumour_mutations$Case=="PD40735"&tumour_mutations$Mutation_type=="Substitution",]
tumour_samples_left=sample_info$Sample[sample_info$Case==patient&sample_info$Sample_type_1=="Neoplasm"&grepl("left",sample_info$Sample_type_2)]
tumour_somatic_left=sum(table(tumour_sub$ID[tumour_sub$Sample%in%tumour_samples_left])==length(tumour_samples_left)) 
tumour_samples_right=sample_info$Sample[sample_info$Case==patient&sample_info$Sample_type_1=="Neoplasm"&grepl("right",sample_info$Sample_type_2)]
tumour_somatic_right=sum(table(tumour_sub$ID[tumour_sub$Sample%in%tumour_samples_right])==length(tumour_samples_right)) 

#Now for the shared/nephrogenic branches
normal_muts_sub=normal_muts[normal_muts$Patient==patient,]
nephrogenic_muts_left_all=sum(rowSums(normal_muts_sub[,tumour_samples_left]!="NS")==length(tumour_samples_left)&
                                normal_muts_sub$Type=="Nephrogenic")
nephrogenic_muts_right_all=sum(rowSums(normal_muts_sub[,tumour_samples_right]!="NS")==length(tumour_samples_right)&
                                 normal_muts_sub$Type=="Nephrogenic")

nephrogenic_muts_shared=sum(rowSums(normal_muts_sub[,c(tumour_samples_right,tumour_samples_left)]!="NS")==
                              length(c(tumour_samples_right,tumour_samples_left))&
                                    normal_muts_sub$Type=="Nephrogenic") 
nephrogenic_muts_right=nephrogenic_muts_right_all-nephrogenic_muts_shared 
nephrogenic_muts_left=nephrogenic_muts_left_all-nephrogenic_muts_shared

embryonic_muts=sum(rowSums(normal_muts_sub[,c(tumour_samples_right,tumour_samples_left)]!="NS")==
                              length(c(tumour_samples_right,tumour_samples_left))&
                              normal_muts_sub$Type=="Embryonic")
