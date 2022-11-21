remove(list = ls())
source('Code/Defined_Functions.R')

############################ Calculate TRIM's ssGSEA score ############################
setwd('expression/ssgsea')

TRIM=list(Trimdata$ENSG)
names(TRIM)='TRIM'

myssgsea_all<-data.frame(Immune_cell=character(),Sample=character(),
                         Score=numeric(),Sample_Type=character(),Cancer=character())
for(i in 1:length(cancer_names)){
   cancer0=cancer_names[i]
   cancerpath=paste0(workdir3,'expression/expression_clear/')
   
   cancerdata=paste0(cancer0,'_cancer.txt')
   expdat_cancer=read.csv(paste0(cancerpath,cancerdata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
   expdat_cancer=expdat_cancer %>% dplyr::select(-Hugo_Symbol)
   expdat=expdat_cancer
   
   normaldata=paste0(cancer0,'_normal.txt')
   if(normaldata %in% list.files(cancerpath)){
   expdat_normal=read.csv(paste0(cancerpath,normaldata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
   expdat_normal=expdat_normal %>% dplyr::select(-Hugo_Symbol)
   expdat<-cbind(expdat_cancer,expdat_normal)
   }
   
   expdat2=as.matrix(expdat)
   myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Poisson',abs.ranking=TRUE)
   write.csv(myssgsea,paste0(cancer0,'_ssgsea_score.csv'),quote = F)

   myssgsea2=reshape2::melt(myssgsea)
   colnames(myssgsea2)=c('Immune_cell','Sample','Score')
   myssgsea2$Sample_Type='Tumor'
   if(normaldata %in% list.files(cancerpath)){
   myssgsea2[which(myssgsea2$Sample %in% colnames(expdat_normal)),'Sample_Type']<-'Normal'
   }
   myssgsea2$Cancer=cancer0
   myssgsea_all<-rbind(myssgsea_all,myssgsea2)
}

write.csv(myssgsea_all,paste0('all_ssgsea_score.csv'),quote = F,row.names = F)

#----------- boxplot-TRIMs score

myssgsea_all=readtxt_1(workdir3,paste0('expression/ssgsea/all_ssgsea_score.csv'),sep0 = ',')
myssgsea_all$Cancer=factor(myssgsea_all$Cancer,levels = cancer_names)
myssgsea_all=myssgsea_all[order(myssgsea_all$Cancer),]

myssgsea_all_compare0= unique(myssgsea_all[which(myssgsea_all$Sample_Type == 'Normal'),'Cancer'])
myssgsea_all_compare=myssgsea_all %>% filter(Cancer %in% myssgsea_all_compare0)

colors=c(brewer.pal(3,"BrBG")[1],
         brewer.pal(9,"Blues")[c(6,9)],
         brewer.pal(9,"Blues")[4],
         brewer.pal(11,"PiYG")[c(5:4)],brewer.pal(9,"Set3")[9],
         brewer.pal(11,"PiYG")[7],
         brewer.pal(11,"BrBG")[4:2],
         brewer.pal(9,"PuRd")[2:7],brewer.pal(9,"Purples")[8:9],
         brewer.pal(9,"GnBu")[5:9],
         brewer.pal(9,"Greens")[4:7],
         brewer.pal(9,"Oranges")[4:6],
         brewer.pal(9,"Greys")[7],
         brewer.pal(9,"Set1")[3]
)
names(colors)<-cancer_names

myssgsea_all2=myssgsea_all %>% dplyr::filter(Sample_Type=='Tumor')
Data_summary <- summarySE(myssgsea_all2, measurevar="Score", groupvars=c('Cancer','Sample_Type'))
colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')
p=ggplot(data=myssgsea_all2,aes(x=Cancer,y=Score,color=Cancer))+
  stat_boxplot(geom="errorbar",width=0.15,aes(x=Cancer,color=Cancer),position = position_dodge(0.9))+
  geom_point(data = Data_summary,aes(x=Cancer,y=median,color='black'),alpha = 0.8,pch = 19,position = position_dodge(0.9),size = 1.5)+
  geom_jitter(aes(color=Cancer),size=1)+
  ylab('TRIMs score')+xlab('')+
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
        axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = 'black',
                                  size=10,
                                  hjust=0.5)
        )+
  guides(colour = 'none')

pdf('TRIMs_score_allcancer.pdf',width = 5,height = 8)
p+coord_flip()
dev.off()

#----------- boxplot-TRIMs score TvsN

Data_summary <- summarySE(myssgsea_all_compare, measurevar="Score", groupvars=c("Cancer","Sample_Type"))
colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')

P3 <- ggplot(data=Data_summary, aes(x=Cancer,color=Sample_Type)) + 
  geom_errorbar(data = Data_summary,aes(ymin = quantile_25, ymax=quantile_75,group=Sample_Type),
                width=0.2, 
                position=position_dodge(0.9), 
                alpha = 1,
                size=1) + 
  theme_bw()+
  
  geom_point(data = Data_summary,aes(x=Cancer, y=median),pch=19,position=position_dodge(0.9),size=6) +
  scale_color_manual(values = c("#1F78B4","#CB181D"))+
  theme(
    axis.text.x=element_text(size=10), 
    axis.text.y=element_text(size=10), 
    axis.title.y=element_text(size=15), 
    axis.title.x=element_text(size=15), 
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank())+  
  xlab("Cancer")+ylab("TRIMs score")+ 
  geom_vline(xintercept=c(seq(1.5,17.5,by=1)), linetype="dotted") 

pdf('TRIMs_score_TvsN.pdf',width = 12,height = 5)
print(P3)
dev.off()

########################################################################

############################ TRIMs score group ############################
library("survival")
library("survminer")

covariates2='Score'

clinicdata_all_score_cut0=data.frame()


result_score=data.frame(p.val=numeric(),HR=numeric(),up95=numeric(),low95=numeric())

cancer=sort(unique(clinicdata_all_score$Cancer))

for(i in 1:length(cancer)){
cancer0=as.character(cancer[i])
data0=clinicdata_all_score %>% filter(Cancer==cancer0)

res.cut <- surv_cutpoint(data0, 
                         time = "days_to_last_follow_up", 
                         event = "vital_status",
                         variables = covariates2)

res.cat <- surv_categorize(res.cut)
res.cat[,covariates2]
res.cat[,covariates2]<-factor(res.cat[,covariates2],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$Tumor_Sample_Barcode
res.cat$Cancer=cancer0

clinicdata_all_score_cut0=rbind(clinicdata_all_score_cut0,res.cat)

univ_formulas=as.formula(paste('Surv(days_to_last_follow_up, vital_status)~', 'Score_cut'))
fit <- survfit(univ_formulas, data = res.cat)
data.survdiff <- survdiff(univ_formulas,data = res.cat)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

result=c(p.val,HR,up95,low95,cancer0)
result_score=rbind(result_score,result)

}


colnames(result_score)=c('Pvalue','HR','Up95','Low95','Cancer')
result_score[,1:4]=apply(result_score[,1:4],2,as.numeric)
write.table(result_score,'TRIM_score_Cox.txt',sep='\t',quote = F,row.names = F)


clinicdata_all_score_cut0=clinicdata_all_score_cut0[,c('Score_cut','Sample')]
clinicdata_all_score_cut=merge(clinicdata_all_score,clinicdata_all_score_cut0,
                               by.x='Tumor_Sample_Barcode',by.y='Sample')
clinicdata_all_score_cut=clinicdata_all_score_cut[order(clinicdata_all_score_cut$Cancer),]
write.table(clinicdata_all_score_cut,'clinicdata_all_score_cut.txt',sep='\t',quote = F,row.names = F)

save.image('ssGSEA_cox_file.RData')



#------------boxplot-COX

cox_result_score=readtxt_1('./','TRIM_score_Cox.txt')
cox_result_score$Type='Sig'
cox_result_score[which(cox_result_score$Pvalue > 0.05),'Type']='Not-Sig'
cox_result_score2=cox_result_score
cox_result_score2[,2:4]=log2(cox_result_score2[,2:4])
cox_result_score2$Cancer=factor(cox_result_score2$Cancer,levels = cancer_names)
cox_result_score2=cox_result_score2[order(cox_result_score2$Cancer),]


plotcox=function(){
  P3 <- ggplot(data=cox_result_score2, aes(x=Cancer,color=Type)) + 
    geom_errorbar(data = cox_result_score2,aes(ymin = round(Low95,4), ymax=round(Up95,4)),
                  width=0.3, 
                  position=position_dodge(0.9), 
                  alpha = 1,
                  size=1) + 
    theme_bw()+
    
    geom_point(data = cox_result_score2,aes(x=Cancer, y=round(HR,4)),pch=19,position=position_dodge(0.9),size=6) +
    scale_color_manual(values = c('#999999','#D6604D') )+
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank()
    )+  
    xlab("Cancer")+ylab("log2 TRIM score HR(95% CI)")+ 
    geom_hline(yintercept=0, linetype="dotted")
  return(P3)
}
P3=plotcox()
pdf('TRIM_score_COX.pdf',width = 15,height = 5)
P3
dev.off()

#----------boxplot-TRIMs score ratio

library(tidyverse)
library(dplyr)
library(viridis)
tumorOrder=read.xlsx(paste0(work_path,'tumorOrder.xlsx'),sheet=2,colNames =F)

tumorClass=read.csv(paste0(work_path,'Cancer_Class.csv'),header=T)
tumorClass=tumorClass %>% dplyr::select(Cancer_System,Cohort)

clinicdata_all_score_cut=read.csv('clinicdata_all_score_cut.txt',sep='\t',header=T)
clinicdata_all_score_cut2= clinicdata_all_score_cut %>% 
  dplyr::select(Cancer,Score_cut,Tumor_Sample_Barcode) %>% 
  group_by(Cancer,Score_cut) %>% dplyr::summarize(value=n())

clinicdata_all_score_cut3=merge(clinicdata_all_score_cut2,tumorClass,by.x='Cancer',by.y='Cohort')
clinicdata_all_score_cut3$Cancer=factor(clinicdata_all_score_cut3$Cancer,levels = tumorOrder$X1)
clinicdata_all_score_cut3=clinicdata_all_score_cut3[order(clinicdata_all_score_cut3$Cancer),]

write.csv(clinicdata_all_score_cut3,'TRIMs_score_ratio.csv',row.names = F)


data=clinicdata_all_score_cut3
data$value=log2(data$value)
data$Cancer=factor(data$Cancer,levels = tumorOrder$X1)
data=data[order(data$Cancer),]
ord=unique(data$Cancer_System)
data$Cancer_System=factor(data$Cancer_System,levels = ord)
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
nObsType <- nlevels(as.factor(data$Score_cut))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Cancer_System)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Cancer_System <- rep(levels(data$Cancer_System), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Cancer_System, Cancer)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- data %>% group_by(id, Cancer) %>% dplyr::summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust2 <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

label_data2=label_data %>% dplyr::select(Cancer,hjust2,angle)
label_data2=as.data.frame(na.omit(label_data2))
label_data2=merge(label_data2,tumorClass,by.x='Cancer',by.y='Cohort')
label_data2=label_data2 %>% dplyr::select(Cancer_System,hjust2,angle)
label_data3=label_data2 %>% 
  group_by(Cancer_System) %>% 
  dplyr::summarize(angle2=mean(angle),hjust2=mean(hjust2))
label_data3$hjust2=ceiling(label_data3$hjust2)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Cancer_System) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

base_data2=merge(base_data,label_data3,by='Cancer_System')

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=Score_cut), stat="identity", alpha=0.7) +
  scale_fill_brewer(palette='Set1')+
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-10,max(label_data$tot, na.rm=T)+20) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot+2, label=Cancer, hjust=hjust2), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data2, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data2, aes(x = title, y = -3, label=Cancer_System,hjust=hjust2),angle= base_data2$angle2-90 , colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)

pdf('TRIMs_score_ratio_incancer.pdf',width = 8,height = 8)
print(p)
dev.off()


#-------------Differential expression analysis between TRIMs score high and low group
clinicdata_all_score_cut=readtxt_1(workdir3,filename = 'expression/ssgsea/clinicdata_all_score_cut.txt')
score_cut=clinicdata_all_score_cut%>% dplyr::select(Tumor_Sample_Barcode,disease,Score_cut)
gene_symbol=readtxt_1(workdir3,'ENSGall_Symbol.txt')

diff_result_all<-data.frame()
for(i in 1:length(cancer_names)){
   cancer0=cancer_names[i]
   cancerpath=paste0(workdir3,'expression/expression_clear/')
   expdat_cacner=readtxt_1(cancerpath,paste0(cancer0,'_cancer.txt'))
   rownames(expdat_cacner)=expdat_cacner$ENSG
   expdat_cacner2=expdat_cacner[,-c(1,2)]
   colnames(expdat_cacner2)=lapply(colnames(expdat_cacner2),splitsample) %>% unlist()
   expdat=expdat_cacner2
   expdat2=expdat
   cleargene<-rownames(expdat2)[which(apply(expdat2,1,function(v){return((sum(v==0)/length(v))>=0.7)}))]
   if(length(cleargene)>0){
   expdat2=expdat2[-which(rownames(expdat2) %in% cleargene),]
   }
   
   score_cut0=score_cut %>% filter(disease == cancer0)
   sample_high=score_cut0[which(score_cut0$Score_cut=='high'),'Tumor_Sample_Barcode']
   sample_low=score_cut0[which(score_cut0$Score_cut=='low'),'Tumor_Sample_Barcode']
   
   high<-rowMeans(expdat2[,sample_high])
   low<-rowMeans(expdat2[,sample_low])
   foldchange<- high/low
   foldchange2<-as.data.frame(foldchange)
   foldchangelog2<-log2(foldchange)
   
   lg2highexp<-log2(expdat2[,sample_high]+1) 
   lg2lowexp<-log2(expdat2[,sample_low]+1) 
   
   list0<-cbind(c(1:nrow(lg2highexp)),c(1:nrow(lg2lowexp)))
   wilcox_test<-function(list0){
      high0<-as.numeric(lg2highexp[list0[1],])
      low0<-as.numeric(lg2lowexp[list0[2],])
      wtestvalue<-wilcox.test(high0, low0)
      pvalue0<-wtestvalue$p.value
      return(pvalue0)
   }
   pvalue<-apply(list0,1,wilcox_test)
   names(pvalue)<-rownames(lg2highexp)
   pvalue_adjust<-p.adjust(pvalue,"fdr")
   diff_result<-cbind(foldchange2,foldchangelog2,pvalue,pvalue_adjust)
   diff_result$ENSG<-rownames(lg2highexp)
   diff_result$Cancer=cancer0
   diff_result_all=rbind(diff_result_all,diff_result)
   cat(cancer0,sep = "\n")
}
diff_result_all2=merge(diff_result_all,gene_symbol,by.x='ENSG')
write.csv(diff_result_all2,paste0('TRIMscoreLH_DE_allcancer.csv'),quote = F,row.names = F)

########################################################################


############################ FunctionEnrichment ############################
#-------------------GO
library(simplifyEnrichment)
library(GOSemSim)
library(rrvgo)
library(magick)
diff_result_all2=readtxt_1(paste0('TRIMscoreLH_DE_allcancer.csv'),sep0 = ',')
cancerall=unique(diff_result_all2$Cancer)

lapply(cancerall,function(cancer0){
diff_result_all2_0=diff_result_all2[which(diff_result_all2$Cancer == cancer0),]
diff_result_all2_0_up=diff_result_all2_0 %>% filter(pvalue_adjust < 0.05 & foldchangelog2 > 1)
diff_result_all2_0_down=diff_result_all2_0 %>% filter(pvalue_adjust < 0.05 & foldchangelog2 < -1)

eg_up <- bitr(diff_result_all2_0_up$Hugo_Symbol, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")
head(eg_up)
go_up <- enrichGO(unique(eg_up$ENTREZID), 
               OrgDb = org.Hs.eg.db, 
               ont='BP',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID')
saveRDS(go_up,paste0('result/',cancer0,'_up.Rds'))

eg_down <- bitr(diff_result_all2_0_down$Hugo_Symbol, 
              fromType="SYMBOL", 
              toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
              OrgDb="org.Hs.eg.db")
head(eg_down)
go_down <- enrichGO(unique(eg_down$ENTREZID), 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,
                  keyType = 'ENTREZID')
saveRDS(go_down,paste0('result/',cancer0,'_down.Rds'))

})


#----------------GOplot
fileall=list.files('result','Rds',recursive =T,full.names=T)
fileall_up=grep('up',fileall,value =T)
fileall_down=grep('down',fileall,value =T)

newdatat_up <- c()
for (i in 1:length(fileall_up)) {
  fileall0=fileall_up[i]
  fileall0_1=gsub('result/','',fileall0)
  fileall0_1=unlist(strsplit(fileall0_1,'\\.'))[1]
  tryCatch({
  temp <- readRDS(fileall0)
  dat2=as.data.frame(temp) 
  dat2$Cancer=fileall0_1
  newdatat_up=rbind(newdatat_up,dat2)
  },error = function(e){
    print(fileall0_1)
  }
  )
}
pdf(paste0('cancer_up.pdf'),width=12,height = 9)

dat2=newdatat_up %>% 
  filter(!is.na(qvalue) & qvalue <= 0.01 & Count>5 & Count <500) %>% arrange(desc(Count))
length(unique(dat2$ID))

mat = GO_similarity(unique(dat2$ID),ont='BP')
df = simplifyGO(mat,column_title='cancer_up',fontsize=c(6,20))

dev.off()

newdatat_down <- c()
for (i in 1:length(fileall_up)) {
  fileall0=fileall_down[i]
  fileall0_1=gsub('result/','',fileall0)
  fileall0_1=unlist(strsplit(fileall0_1,'\\.'))[1]
  tryCatch({
    temp <- readRDS(fileall0)
    dat2=temp %>% arrange(desc(Count))
    getrow=min(nrow(dat2),50)
    dat2=dat2[1:getrow,]
    dat2$Cancer=fileall0_1
    newdatat_down=rbind(newdatat_down,dat2)
  },error = function(e){
    print(fileall0_1)
  }
  )
}

pdf(paste0('cancer_down.pdf'),width=12,height = 9)
dat2=newdatat_down %>% 
  filter(!is.na(qvalue) & qvalue <= 0.0001 & Count>5 & Count <500) %>% arrange(desc(Count))
length(unique(dat2$ID))

mat = GO_similarity(unique(dat2$ID),ont='BP')
df = simplifyGO(mat,column_title='cancer_down',fontsize=c(10,20))

dev.off()


#---------------------GSEA

library(clusterProfiler)
kegmt<-read.gmt("h.all.v7.4.symbols.gmt")
ENSGall_Symbol=readtxt_1(workdir3,'ENSGall_Symbol.txt',sep0 = '\t',header0 = T)

cancerall_spearman=data.frame(Gene=as.character(),Cancer=as.character(),
                              Rho=as.numeric(),Pvalue=as.numeric(),Padjust=as.numeric())
for(i in 1:length(cancer_names)){
   cancer0=cancer_names[i]
   cancerpath=paste0(workdir3,'expression/expression_clear/')
   
   cancerdata=paste0(cancer0,'_cancer.txt')
   expdat_cancer=read.csv(paste0(cancerpath,cancerdata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
   expdat_cancer2=expdat_cancer[,-1]
   expdat=expdat_cancer2
   
   expdat2=apply(expdat, 2, as.numeric)
   rownames(expdat2)=rownames(expdat)
   expdat3=log2(expdat2+1)
   
   myssgsea=read.csv(paste0(workdir3,'expression/ssgsea/',cancer0,'_ssgsea_score.csv'),sep=',',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
   myssgsea2=myssgsea[1,-1] %>% as.numeric()
   names(myssgsea2)=colnames(myssgsea)[-1]
   
   sample_intersect=intersect(colnames(expdat3),names(myssgsea2))
   myssgsea2=myssgsea2[sample_intersect]
   expdat3=expdat3[,sample_intersect]
   
   cancer0_spearman=data.frame()
   for(j in 1:nrow(expdat2)){
   cor_result=cor.test(expdat3[j,],myssgsea2,method='spearman')
   pvalue=cor_result$p.value
   rho=cor_result$estimate
   cancer0_spearman0=c(rownames(expdat3)[j],cancer0,rho,pvalue)
   cancer0_spearman=rbind(cancer0_spearman,cancer0_spearman0)
   }
   colnames(cancer0_spearman)=c('Gene','Cancer','R','Pvalue')
   cancer0_spearman$Padjust=p.adjust(cancer0_spearman$Pvalue,method ='fdr')

   dat=merge(cancer0_spearman,ENSGall_Symbol,by.x='Gene',by.y='ENSG')
   dat_1=dat %>% dplyr::select(Hugo_Symbol,R)
   dat_1=dat_1[order(dat_1$R,decreasing = T),]
   if(nrow(dat_1)<1){
      cat(paste0(genetrp[j],'_Rnot\n'),file='notattend_all.txt',append=T)
   }else{
      dupid = unique(dat_1$Hugo_Symbol[which(duplicated(dat_1$Hugo_Symbol))])
      
      if(length(dupid)>0){
         dat_2= dat_1[-which(dat_1$Hugo_Symbol %in% dupid),]
      }else{
         dat_2=dat_1
      }
      dat_2$R=as.numeric(dat_2$R)
      dat_2=dat_2[order(dat_2$R,decreasing = T),]
      dat_2=dat_2[which(dat_2$Hugo_Symbol!=''),]
      dat_2=na.omit(dat_2)
      
      print(max(dat_2$R))
      kgene=dat_2
      geneList<-kgene[,2] 
      names(geneList)=kgene[,1]
      geneList=sort(geneList,decreasing = T) 
      KEGG<-tryCatch({
         clusterProfiler::GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 1,eps=0,
                               minGSSize=5,maxGSSize=10000)
      },error = function(e){
         cat(paste0(cancer0,'_KEGGnot1\n'),file='notattend_all.txt',append=T)
      })
      if(!is.null(KEGG)){
         if(nrow(KEGG@result)>0){
            saveRDS(KEGG,paste0(cancer0,'_GSEA_all.Rds'))
            cat(paste0(cancer0,'\n'),file='attend_all.txt',sep='\t',append=T)
         }else{ cat(paste0(cancer0,'_KEGGnot2\n'),file='notattend_all.txt',sep='\t',append=T)}
      }
   }

   cancerall_spearman=rbind(cancerall_spearman,cancer0_spearman)
}

colnames(cancerall_spearman)=c('Gene','Cancer','R','Pvalue','Padjust')
write.table(cancerall_spearman,'TRIMscore_other_spearman.txt',sep='\t',quote = F,row.names = F)

cancerGSEA_all=data.frame()
for (i in 1:length(cancer_names)) {
   cancer0=cancer_names[i]
   cancerGSEA0=readRDS(paste0(cancer0,'_GSEA_all.Rds'))
   cancerGSEA=cancerGSEA0@result %>% dplyr::select('ID','enrichmentScore','NES','pvalue','p.adjust','qvalues')
   cancerGSEA$Cancer=cancer0
   cancerGSEA_all=rbind(cancerGSEA,cancerGSEA_all)
}

write.table(cancerGSEA_all,'TRIMscore_GSEAresult.txt',sep='\t',quote = F,row.names = F)

#----------------------GSEAplot
library(stringr)
cancerGSEA_all=readtxt_1(filename='TRIMscore_GSEAresult.txt')
cancerGSEA_all$ID=gsub('HALLMARK_','',cancerGSEA_all$ID)
cancerGSEA_all$ID=gsub('_',' ',cancerGSEA_all$ID)
cancerGSEA_all$ID=str_to_title(cancerGSEA_all$ID)
cancerGSEA_all$ID=factor(cancerGSEA_all$ID,levels = sort(unique(cancerGSEA_all$ID)))
cancerGSEA_all$Cancer=factor(cancerGSEA_all$Cancer,levels = cancer_names)
cancerGSEA_all=cancerGSEA_all[order(cancerGSEA_all$ID),]
cancerGSEA_all_sig=cancerGSEA_all %>% filter(qvalues < 0.05)

result_GSEA_point <- ggplot()+
  geom_tile(data=cancerGSEA_all,aes(Cancer,ID),fill="white",color=NA)+
  geom_point(data=cancerGSEA_all,aes(Cancer,ID,size=-log10(qvalues),color=NES),
             shape=20)+#scale_size_continuous(range = c(4,6))+
  scale_color_gradient2(low = 'blue', mid = "white", high = 'red',midpoint = 0,limits = c(-4,4.5))+
  geom_point(data=cancerGSEA_all_sig,aes(Cancer,ID,size=-log10(qvalues)),#alpha=2,
             shape=21)+#scale_size_continuous(range = c(4,6))+
  xlab('Cancer')+ylab('Hallmark')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust = 1),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_line(linetype=1,color='grey'),
        axis.ticks = element_line(linetype=2,color='grey'),
        axis.text=element_text(face="bold"),
        panel.grid=element_line(linetype =2),
        legend.position = 'top'
  )

pdf('TRIM_score_GSEA_merge.pdf',
                width = 8, height = 12)
print(result_GSEA_point)
dev.off()


dat=cancerGSEA_all
dat_2=dat %>% filter(qvalues < 0.05)
dat_H=dat_2 %>% filter(NES>0)
hist_H=as.data.frame(table(dat_H$ID))
colnames(hist_H)=c('Pathway','Positive')
dat_L=dat_2 %>% filter(NES<0)
hist_L=as.data.frame(table(dat_L$ID))
colnames(hist_L)=c('Pathway','Negative')
hist=merge(hist_H,hist_L,by='Pathway',all=T)
hist[is.na(hist)]=0
hist2=melt(hist)

result_GSEA_bar=ggplot(hist2, aes(
  x = factor(Pathway,levels = unique(Pathway)),             
  y = ifelse(variable == "Positive", value, -value), 
  fill = variable)) +
  geom_bar(stat = 'identity')+                                
  ylab('Number of Cancer')+
  xlab('Pathway')+
  scale_fill_manual(values = c('#D6604D','#4393C3') )+  
  scale_y_continuous(                                         
    labels = abs,                                             
    expand = expansion(mult = c(0.1, 0.1)))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),

        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        legend.title = element_blank()
  )

result_GSEA_bar=result_GSEA_bar+coord_flip()+
  geom_text(                                                  
    aes(label=value,                                          
        # vjust = ifelse(variable == "Up", -0.5, 1),          
        hjust = ifelse(variable == "Positive", -0.4, 1.1)
    ),
    size=2
  )
pdf('GSEA_TRIMscore.pdf',width = 10, height = 10)
print(result_GSEA_bar)
dev.off()

pdf('TRIM_score_GSEA_merge_bar_point.pdf',
                width = 12, height = 12)
result_GSEA_bar2=result_GSEA_bar+labs(x='')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
result_GSEA_point+result_GSEA_bar2+plot_layout(widths = c(3, 1))

dev.off()

dat=cancerGSEA_all_sig

#  Positive
condition='Positive'
dat_H=dat %>% filter(NES >0)
dat_H2=dat_H %>% group_by(ID) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
top=dat_H2$ID
for(k in 1:length(top)){
  top1=top[k]
  plotH=MyGSEAplot(dat_H,top1,condition)
  top2=gsub(' ','_',top1)
  pdf(paste0('GSEAPLOT/',top2,'/',top2,'_',condition,'.pdf'),height = 4,width=6)
  print(plotH)
  dev.off()
  try(dev.off())
}


#  Negative
condition='Negative'
dat_L=dat %>% filter(NES < 0)
dat_L2=dat_L %>% group_by(ID) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
top=dat_L2$ID
for(k in 1:length(top)){
  top1=top[k]
  plotL=MyGSEAplot(dat_L,top1,condition)
  top2=gsub(' ','_',top1)
  pdf(paste0('GSEAPLOT/',top2,'/',top2,'_',condition,'.pdf'),height = 4,width=6)
  print(plotL)
  dev.off()
  try(dev.off())
}

########################################################################


############################ TRIMs score - Immune_cell ############################

immune_data<-readtxt_1(filename = paste0('infiltration_estimation_for_tcga.csv'),sep0=',',header0=T)
rownames(immune_data)=immune_data$cell_type
Sample0=immune_data$cell_type
colnames(immune_data)=gsub('\\.','_',colnames(immune_data))
CIBERSORTcol=grep('CIBERSORT$',colnames(immune_data))
immune_data_CIBERSORT=immune_data[,CIBERSORTcol]
colnames(immune_data_CIBERSORT)=gsub('_CIBERSORT','',colnames(immune_data_CIBERSORT))

cancertype=readtxt_1(workdir3,'all_cancer_sample.txt',sep0='\t',header0=T)
cancertype$Sample=lapply(cancertype$Sample,function(a){
a1=unlist(strsplit(a,'-'))[1:4]
a1[4]=substr(a1[4],1,nchar(a1[4])-1)
b=paste(a1,collapse='-')
return(b)
}) %>% unlist()

cancertype_1=cancertype %>% filter(Sample %in% rownames(immune_data_CIBERSORT))
cancertype_1=cancertype_1[order(cancertype_1$Cancer),] %>% distinct()
immune_data_CIBERSORT_1=t(immune_data_CIBERSORT)
immune_data_CIBERSORT_1=immune_data_CIBERSORT_1[,cancertype_1$Sample]


immune_data_CIBERSORT2=reshape2::melt(as.matrix(immune_data_CIBERSORT))
colnames(immune_data_CIBERSORT2)=c('Sample','Immune_cell','Value')


immune_data_CIBERSORT3=merge(immune_data_CIBERSORT2,cancertype,by='Sample')

myssgsea=read.csv(paste0(workdir3,'expression/ssgsea/all_ssgsea_score.csv'),sep=',',header=T)

cancerall_spearman=data.frame(Gene=as.character(),Cancer=as.character(),
                              Rho=as.numeric(),Pvalue=as.numeric(),Padjust=as.numeric())

for(i in 1:length(cancer_names)){
cancer0=cancer_names[i]

immune_data_CIBERSORT3_cancer0=immune_data_CIBERSORT3 %>% filter(Cancer==cancer0) %>% distinct()
myssgsea_cancer0=myssgsea %>% filter(Cancer==cancer0) %>% distinct() %>% arrange(Sample)
Immune_cell_all=as.character(unique(immune_data_CIBERSORT3_cancer0$Immune_cell)) %>% sort()
# j=1
cancer0_spearman=data.frame()
for(j in 1:length(Immune_cell_all)){
Immune_cell0=Immune_cell_all[j] %>% as.character()

immune_data_CIBERSORT3_cancer_1=immune_data_CIBERSORT3_cancer0 %>% 
  filter(Immune_cell==Immune_cell0 & Sample %in% myssgsea_cancer0$Sample) %>% 
  distinct() %>% 
  arrange(Sample)
myssgsea_cancer_1=myssgsea_cancer0 %>% filter(Sample %in% immune_data_CIBERSORT3_cancer_1$Sample) %>% distinct()

dupsample=c(immune_data_CIBERSORT3_cancer_1$Sample[which(duplicated(immune_data_CIBERSORT3_cancer_1$Sample))],
myssgsea_cancer_1$Sample[which(duplicated(myssgsea_cancer_1$Sample))])
if(length(dupsample)>0){
immune_data_CIBERSORT3_cancer_1=immune_data_CIBERSORT3_cancer_1[-which(immune_data_CIBERSORT3_cancer_1$Sample %in% dupsample),]
myssgsea_cancer_1=myssgsea_cancer_1[-which(myssgsea_cancer_1$Sample %in% dupsample),]
}

rownames(immune_data_CIBERSORT3_cancer_1)=immune_data_CIBERSORT3_cancer_1$Sample
rownames(myssgsea_cancer_1)=myssgsea_cancer_1$Sample

sampleorder=as.character(myssgsea_cancer_1$Sample)

cor_result=cor.test(immune_data_CIBERSORT3_cancer_1[sampleorder,'Value'],myssgsea_cancer_1[sampleorder,'Score'],method='spearman')
pvalue=cor_result$p.value
rho=cor_result$estimate
cancer0_spearman0=c(Immune_cell0,cancer0,rho,pvalue)
cancer0_spearman=rbind(cancer0_spearman,cancer0_spearman0)
}
colnames(cancer0_spearman)=c('Immune_cell','Cancer','R','Pvalue')
cancer0_spearman$Padjust=p.adjust(cancer0_spearman$Pvalue,method ='fdr')
cancerall_spearman=rbind(cancerall_spearman,cancer0_spearman)

}
colnames(cancerall_spearman)=c('Immune_cell','Cancer','R','Pvalue','Padjust')
cancerall_spearman[,3:5]=apply(cancerall_spearman[,3:5],2,as.numeric)
write.table(cancerall_spearman,'TRIMscore_CIBERSORT_spearman.txt',sep='\t',quote = F,row.names = F)


#--------------dotplot
cancerall_spearman=readtxt_1(filename = 'TRIMscore_CIBERSORT_spearman.txt')
cancerall_spearman$Cancer=factor(cancerall_spearman$Cancer,levels = cancer_names)
cancerall_spearman_sig=cancerall_spearman %>% filter(Padjust <= 0.05)

result_spearman_immune_point <- ggplot()+
   geom_tile(data=cancerall_spearman,aes(Cancer,Immune_cell),fill="white",color=NA)+
   scale_color_gradient2(low = 'blue', mid = "white", high = 'red',midpoint = 0,limits = c(-0.5,0.5))+#
   geom_point(data=cancerall_spearman_sig,aes(Cancer,Immune_cell,size=-log10(Padjust)),#alpha=2,
              shape=16)+scale_size_continuous(range = c(4,8))+
  geom_point(data=cancerall_spearman,aes(Cancer,Immune_cell,size=-log10(Padjust),color=R),
             shape=20)+scale_size_continuous(range = c(4,8))+
   ggtitle('')+
   theme_bw()+
   theme(axis.text.x = element_text(angle=90,hjust = 1,size=15),
   axis.text.y = element_text(size=15),
   axis.title.y=element_text(size = 18),
   axis.title.x=element_text(size = 18),
         plot.title = element_text(hjust = 0.5,size = 18),
		 legend.text=element_text(size = 15),
		 legend.title=element_text(size = 18),
         axis.line=element_line(linetype=1,color='grey'),
         axis.ticks = element_line(linetype=2,color='grey'),
         panel.grid=element_line(linetype =2),
         legend.position = 'right'
   )


pdf('TRIMscore_CIBERSORT_cor_dotplot.pdf',
                width = 12, height = 13)
print(result_spearman_immune_point)
dev.off()


#----------------boxplot
immune_data_CIBERSORT4=reshape2::melt(as.matrix(immune_data_CIBERSORT))
colnames(immune_data_CIBERSORT4)=c('Sample','Immune_cell','CIBERSORT')
immune_data_CIBERSORT4$Sample=lapply(immune_data_CIBERSORT4$Sample,function(a){
a=as.character(a)
a1=unlist(strsplit(a,'-'))[1:3]
b=paste(a1,collapse='-')
return(b)
}) %>% unlist()


TCGAall_clinic_cut=readtxt_1('../','clinicdata_all_score_cut.txt',sep0='\t',header0=T)

mergedata=merge(TCGAall_clinic_cut,immune_data_CIBERSORT4,by.x='Tumor_Sample_Barcode',by.y='Sample')

colmy=c("#CB181D","#1F78B4")
names(colmy)<-c('high','low')
p <- ggplot(mergedata,aes(x=Immune_cell.y,y=CIBERSORT,color=Score_cut))+
stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
geom_boxplot(position = position_dodge(0.9),outlier.shape=NA)+
ylab('CIBERSORT')+xlab('Immune_cell')+
scale_color_manual(values = colmy)+
theme_bw()+
   theme(
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())


pdf(paste0('TRIMscoreLH_CIBERSORT_boxplot.pdf'),width=15,height=6)
p + ylim(0,0.75) +stat_compare_means(aes(group = Score_cut), label = "p.signif",label.y=0.7)
dev.off()


#----------------scatterplot
cancerall_spearman=readtxt_1(filename = 'TRIMscore_CIBERSORT_spearman.txt')
cancerall_spearman$Cancer=factor(cancerall_spearman$Cancer,levels = cancer_names)
cancerall_spearman=cancerall_spearman[order(cancerall_spearman$Cancer),]
immune_cell_all=unique(cancerall_spearman$Immune_cell)

for(i in 1:length(immune_cell_all)){

immune_cell_0=as.character(immune_cell_all[i])

cancerall_spearman_sig=cancerall_spearman %>% 
  filter(Padjust <= 0.05 & Immune_cell %in% c(immune_cell_0))%>% arrange(desc(R))
cancer_sig=sort(unique(cancerall_spearman_sig$Cancer))

immune_data_CIBERSORT3_M1=immune_data_CIBERSORT3 %>% 
  filter(Immune_cell == immune_cell_0 & Cancer %in% cancer_sig)
CIBERSORT3_myssgsea_M1=merge(immune_data_CIBERSORT3_M1,myssgsea,by=c('Sample','Cancer'))
CIBERSORT3_myssgsea_M1$Cancer=factor(CIBERSORT3_myssgsea_M1$Cancer,levels = cancer_names)
CIBERSORT3_myssgsea_M1=CIBERSORT3_myssgsea_M1[order(CIBERSORT3_myssgsea_M1$Cancer),]


p=ggplot(CIBERSORT3_myssgsea_M1, aes(Score, Value,fill=Cancer,color=Cancer))+
  geom_point(alpha=0.6)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=4)+
  theme_bw() + 
  scale_color_manual(values = colorsmy[cancer_sig])+
  scale_fill_manual(values = colorsmy[cancer_sig])+
  ylab(immune_cell_0)+xlab('TRIM score')+ggtitle(immune_cell_0)+
  theme(plot.title = element_text(hjust = 0.5))

pdf(paste0('scatterplot/Immune_',immune_cell_0,'.pdf'),width=9,height = 8)
print(p)
dev.off()
}
########################################################################




############################ TRIMs score - Immunomodulators ############################

Immunomodulators=read.xlsx('Immunomodulators.xlsx',sheet=1)
Immunomodulators=Immunomodulators %>% dplyr::select(Gene,HGNC.Symbol,Friendly.Name,Entrez.ID,Gene.Family,
                                                    Super.Category,Immune.Checkpoint,Function)

library(org.Hs.eg.db)

columns(org.Hs.eg.db)
ensembls <- mapIds(org.Hs.eg.db, keys = as.character(unique(Immunomodulators$Entrez.ID)), 
                   keytype = "ENTREZID", column="ENSEMBL")
ensembls2=as.data.frame(unlist(ensembls))
ensembls2$Entrez.ID=rownames(ensembls2)
colnames(ensembls2)[1]='ENSEMBL'

Immunomodulators=merge(Immunomodulators,ensembls2,by='Entrez.ID')
Immunomodulators=Immunomodulators %>% arrange(Super.Category)
Immunomodulators$HGNC.Symbol=factor(Immunomodulators$HGNC.Symbol,levels = Immunomodulators$HGNC.Symbol)

cancerall_spearman=readtxt_1(workdir3,'expression/ssgsea/GSEA/TRIMscore_other_spearman.txt')

cancerall_spearman$Cancer=factor(cancerall_spearman$Cancer,levels = cancer_names)
cancerall_spearman_immune=cancerall_spearman %>% filter(Gene %in% Immunomodulators$ENSEMBL)
cancerall_spearman_immune=merge(cancerall_spearman_immune,Immunomodulators,by.x='Gene',by.y='ENSEMBL')
cancerall_spearman_immune=cancerall_spearman_immune %>% arrange(HGNC.Symbol)
cancerall_spearman_immune_sig=cancerall_spearman_immune %>% filter(Padjust < 0.05)

#-----------heatmap
cancerall_spearman_immune2=cancerall_spearman_immune %>% dplyr::select(Cancer,HGNC.Symbol,R)
cancerall_spearman_immune2<- cancerall_spearman_immune2 %>% tidyr::spread(Cancer,R)
rownames(cancerall_spearman_immune2)=cancerall_spearman_immune2$HGNC.Symbol
cancerall_spearman_immune2=cancerall_spearman_immune2 %>% dplyr::select(-HGNC.Symbol)


annotation_row3=data.frame(Class=factor(Immunomodulators$Super.Category))
library(RColorBrewer)
mycols<-brewer.pal(8,"Set2")
annotation_colors=mycols[1:7]
names(annotation_colors)=sort(unique(annotation_row3$Class))
annotation_colors=list(Class=annotation_colors)
rownames(annotation_row3)<-Immunomodulators$HGNC.Symbol

cancerall_spearman_immune3=cancerall_spearman_immune %>% dplyr::select(Cancer,HGNC.Symbol,Padjust)
cancerall_spearman_immune3<- cancerall_spearman_immune3 %>% spread(Cancer,Padjust)
rownames(cancerall_spearman_immune3)=cancerall_spearman_immune3$HGNC.Symbol
cancerall_spearman_immune3=cancerall_spearman_immune3 %>% dplyr::select(-HGNC.Symbol)


data_norm<-cancerall_spearman_immune3
data_norm2<-cancerall_spearman_immune2
for(i in 1:dim(data_norm)[1]){
   for(j in 1:dim(data_norm)[2]){
      ddp<-as.numeric(data_norm[i,j])
      ddf<-as.numeric(data_norm2[i,j])
      if(is.na(ddp)){
         data_norm[i,j]=''
      }else if(ddp < 0.05){
         data_norm[i,j]='*'
      }else {
         data_norm[i,j]=''
      }
   }
}
library(pheatmap)
library(RColorBrewer)
bk <- c(seq(-0.3,-0.01,by=0.002),seq(0,0.3,by=0.002))

result_spearman_immune_heatmap =pheatmap::pheatmap(cancerall_spearman_immune2,scale = "none",fontsize = 10,
                                                   breaks=bk,
                              annotation_row=annotation_row3,annotation_colors=annotation_colors,
                              main='TRIMscore_Immunomodulators',cluster_rows=FALSE,
                              cluster_cols = FALSE,
                              angle_col = "45",
                              fontsize_number = 10,
                              color=colorRampPalette(rev(c("#CB181D","white","#1F78B4")))(length(bk)),
                              number_color=1,cellheight=15,cellwidth=15)

pdf('TRIMscore_Immunomodulators_cor_pheatmap.pdf',
                width = 12, height = 18)
print(result_spearman_immune_heatmap)
dev.off()


#----------------circos
library('circlize')
dat=cancerall_spearman_immune %>% filter(Padjust < 0.05 & R > 0)
dat$Node_1='TRIMs score'
dat$Type_1='TRIMs score'

df=dat %>% dplyr::select(Node_1,Type_1,HGNC.Symbol,Super.Category,Cancer)
df=df %>% group_by(Node_1,Type_1,HGNC.Symbol,Super.Category) %>% summarise(N=n())
df$HGNC.Symbol=as.character(df$HGNC.Symbol)

pdf(paste0('Spearman_circos.pdf'),width = 10,height = 10)

brand = c(structure(df$Super.Category,names= df$HGNC.Symbol), 
          structure(df$Type_1, names=df$Node_1))

brand = brand[!duplicated(names(brand))]
brand = brand[order(brand, names(brand))]
brand=factor(brand,levels = c('Antigen presentation','Cell adhesion',
                              'Co-inhibitor','Co-stimulator','Ligand',
                              'Receptor','Other','TRIMs score'))
brand=brand[order(brand, names(brand))]

mycols<-c(brewer.pal(7,"Set2"),brewer.pal(3,"Set1")[2])
names(mycols)=c(unique(cancerall_spearman_immune$Super.Category),'TRIMs score')
mycols=mycols[levels(brand)]
prismatic::color(mycols)

brand_color = structure(mycols, names = levels(brand))
model_color = structure(rep(mycols,as.numeric(table(brand))), names = names(brand))
prismatic::color(brand_color)


circos.clear()

kk=table(brand)
kk=kk[which(kk!=0)]
library(circlize)
gap.after = do.call("c", lapply(kk, function(i) c(rep(2, i-1), 8)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

chordDiagram(df[, c('HGNC.Symbol','Node_1','N')], order = names(brand), grid.col = model_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02)
             ),
             link.lty=0.5
)

circos.track(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.index, col = "black", cex = 1, 
              facing = "reverse.clockwise", niceFacing = T)
}, bg.border = NA)

for(b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(sector.index = model, track.index = 1, col = brand_color[b],
                   text = b, text.vjust = -1, niceFacing = TRUE
  )
}
circos.clear()
dev.off()

#----------------barplot
dat=cancerall_spearman_immune
dat_2=dat %>% filter(Padjust < 0.05)
dat_H=dat_2 %>% filter(R>0)
hist_H=dat_H %>% group_by(Cancer) %>% summarise(n=n()) %>% arrange(desc(n))
colnames(hist_H)=c('Cancer','Positive')
dat_L=dat_2 %>% filter(R<0)
hist_L=dat_L %>% group_by(Cancer) %>% summarise(n=n()) %>% arrange(desc(n))
colnames(hist_L)=c('Cancer','Negative')
hist=merge(hist_H,hist_L,by='Cancer',all=T)
hist[is.na(hist)]=0
hist2=melt(hist)
result_spearman_immune_bar=ggplot(hist2, aes(
   x = factor(Cancer,levels = unique(Cancer)),
   y = ifelse(variable == "Positive", value, -value),  
   fill = variable)) +
   geom_bar(stat = 'identity')+
   # coord_flip()+
   ylab('Number of Immunomodulators')+
   xlab('Cancer')+
   scale_fill_manual(values = c('#D6604D','#4393C3') )+
   scale_y_continuous(
      labels = abs,
      expand = expansion(mult = c(0.1, 0.1)))  +
   theme_bw()+
   theme(plot.title = element_text(hjust = 0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'top',
         legend.title = element_blank()
   )
result_spearman_immune_bar=result_spearman_immune_bar+theme(axis.text.x = element_text(angle=45,hjust = 1,colour = 'black'))
pdf('TRIMscore_Immunomodulators_barplot_cancer.pdf',width = 6, height = 3)
print(result_spearman_immune_bar)
dev.off()

dat=cancerall_spearman_immun
dat_2=dat %>% filter(Padjust < 0.05)
dat_H=dat_2 %>% filter(R>0)
hist_H=dat_H %>% group_by(HGNC.Symbol) %>% summarise(n=n()) %>% arrange(desc(n))
colnames(hist_H)=c('HGNC.Symbol','Positive')
dat_L=dat_2 %>% filter(R<0)
hist_L=dat_L %>% group_by(HGNC.Symbol) %>% summarise(n=n()) %>% arrange(desc(n))
colnames(hist_L)=c('HGNC.Symbol','Negative')
hist=merge(hist_H,hist_L,by='HGNC.Symbol',all=T)
hist[is.na(hist)]=0
hist2=melt(hist)
result_spearman_immune_bar=ggplot(hist2, aes(
   x = factor(HGNC.Symbol,levels = rev(unique(HGNC.Symbol))),
   y = ifelse(variable == "Positive", value, -value),
   fill = variable)) +
   geom_bar(stat = 'identity')+
   ylab('Number of Cancer')+
   xlab(NULL)+
   scale_fill_manual(values = c('#D6604D','#4393C3') )+
   scale_y_continuous(
      labels = abs,
      expand = expansion(mult = c(0.1, 0.1)))  +
   theme_bw()+
   theme(plot.title = element_text(hjust = 0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'top',
         legend.title = element_blank()
   )
result_spearman_immune_bar=result_spearman_immune_bar+coord_flip()

pdf('TRIMscore_Immunomodulators_barplot_gene.pdf',width = 3, height =10)
print(result_spearman_immune_bar)
dev.off()

########################################################################



############################ TRIMs score ROC(Normal vs Tumor) ############################
myssgsea=read.csv(paste0(workdir3,'expression/ssgsea/all_ssgsea_score.csv'),sep=',',header=T)

samplenum=read.csv(paste0(workdir3,'lenall_cancer.csv'),sep=',',header=T)
cancer_selected=samplenum[which(samplenum$normal_number>15),'cancername']


pdf('TRIMscore_TvsN_ROC.pdf',width=15,height=15)

ROC.scoreTN=list()
ROC.scoreTN_legend=c()

i=1
cancer0=cancer_selected[i]
df2=myssgsea %>% filter(Cancer==cancer0)
df2$Sample_Type2=df2$Sample_Type
df2$Sample_Type2=gsub('Tumor',1,df2$Sample_Type2)
df2$Sample_Type2=gsub('Normal',0,df2$Sample_Type2)
df2$Sample_Type2=as.numeric(df2$Sample_Type2)

if(length(unique(df2$Sample_Type2))>1){

library(pROC)  
proc<-roc(df2$Sample_Type2,df2$Score,ci=T,smooth=T);proc$auc

ROC.scoreTN[[cancer0]]=proc
ROC.scoreTN_legend=c(ROC.scoreTN_legend,paste0(cancer0,": ",round(proc$auc,2),'(',paste(round(as.numeric(proc$ci)[c(1,3)],2),collapse='-'),')'))

#first plot
plot(proc, col=colorsmy[cancer0], lwd=2, title = "ROC for TRIMscore")
}

for(i in 2:length(cancer_selected)){
cancer0=cancer_selected[i]
df2=myssgsea %>% filter(Cancer==cancer0)
df2$Sample_Type2=df2$Sample_Type
df2$Sample_Type2=gsub('Tumor',1,df2$Sample_Type2)
df2$Sample_Type2=gsub('Normal',0,df2$Sample_Type2)
df2$Sample_Type2=as.numeric(df2$Sample_Type2)

if(length(unique(df2$Sample_Type2))>1){
  library(pROC)  
  proc<-roc(df2$Sample_Type2,df2$Score,ci=T,smooth=T);proc$auc
  
  ROC.scoreTN[[cancer0]]=proc
  ROC.scoreTN_legend=c(ROC.scoreTN_legend,paste0(cancer0,": ",round(proc$auc,2),'(',paste(round(as.numeric(proc$ci)[c(1,3)],2),collapse='-'),')'))

  plot(proc, col=colorsmy[cancer0], lwd=2, add = T)
}else{
  print(cancer0)
}

}

legend("bottomright",
	 ROC.scoreTN_legend,
	 col=colorsmy[cancer_selected],
	 lty=1, lwd=2,bty = "n",cex = 2)

dev.off()

ROC.scoreTN$legend=ROC.scoreTN_legend
saveRDS(ROC.scoreTN,'ROC.scoreTN.Rds')


########################################################################


############################ TRIMs score - drug response ############################
#download drug response data
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library("clusterProfiler")
TCGAall=getGDCprojects()$project_id
TCGAall=grep('TCGA',TCGAall,value=T)

for(i in 1:length(TCGAall)){
  TCGAall[i]
  TCGAbiolinks:::getProjectSummary(TCGAall[i])
  query <- GDCquery(project = TCGAall[i],data.category = 'Clinical',file.type = 'xml')
  GDCdownload(query)
  clinical.drug <- GDCprepare_clinic(query,"drug")
  write.table(clinical.drug,file=paste0(TCGAall[i],"_drug.txt"),sep='\t',quote=F,row.names = F)
}

cancer_names_drug=cancer_names[-which(cancer_names %in% c('LAML'))]
clinicdata_all=data.frame()
# i=1
for (i in 1:length(cancer_names_drug)) {
  cancer0=cancer_names_drug[i]
  clinicdata=readtxt_2(workdir1,paste0('TCGA_clinical','/TCGA-',cancer0,'_drug.txt'))
  clinicdata2=clinicdata %>% dplyr::select(bcr_patient_barcode,project,
                                           therapy_types ,drug_name ,
                                           measure_of_response 
                                           )
  clinicdata_all=rbind(clinicdata_all,clinicdata2)
}
clinicdata_all=as.data.frame(distinct(clinicdata_all))
clinicdata_all$project=gsub('TCGA-','',clinicdata_all$project)
clinicdata_all=as.data.frame(distinct(clinicdata_all))


vroom_write(clinicdata_all,'all_cancer_drugdata.txt')


scorecutdata=read.csv(paste0(workdir3,'expression/ssgsea/clinicdata_all_score_cut.txt'),
                        sep='\t',header=T)
clinic_drug=read.csv(paste0(workdir3,'all_cancer_drugdata.txt'),
                     sep='\t',header=T)
clinic_drug_score=merge(clinic_drug,scorecutdata,by.x=c('bcr_patient_barcode','project'),
                        by.y=c('Tumor_Sample_Barcode','disease'))
clinic_drug_score=clinic_drug_score[!is.na(clinic_drug_score$measure_of_response),]
clinic_drug_score=clinic_drug_score %>% distinct()
as.character(unique(clinic_drug_score$measure_of_response))

clinic_drug_score$Response_class=str_replace_all(clinic_drug_score$measure_of_response,
                                      c(
                                        'Clinical Progressive Disease'='PD',
                                        'Stable Disease'='SD',
                                        'Partial Response'='PR',
                                        'Complete Response'='CR'
                                      )
                                                      )
clinic_drug_score$Response_status=str_replace_all(clinic_drug_score$Response_class,
                                 c(
                                   'PD'='NR',
                                   'SD'='NR',
                                   'PR'='R',
                                   'CR'='R'
                                 )
)

clinicalall2=clinic_drug_score

Rcancer=as.character(unique(clinicalall2$Cancer))
Rcancerfilter=c()
for(i in 1:length(Rcancer)){
Rcancer0=Rcancer[i]
clinicalall2_cancer0=clinicalall2 %>% filter(Cancer==Rcancer0)
print(Rcancer0)
print(dim(clinicalall2_cancer0))
print(table(clinicalall2_cancer0$Response_status))
kk=min(table(clinicalall2_cancer0$Response_status))
if(length(unique(clinicalall2_cancer0$Response_status)) > 1 & kk > 4){
  Rcancerfilter=c(Rcancerfilter,Rcancer0)
}
}

Rcancerfilter=Rcancerfilter[which(Rcancerfilter!='CHOL')]
clinicalall3=clinicalall2 %>% filter(Cancer %in% Rcancerfilter)


#----------barplot
clinicalall3$Response_class=factor(clinicalall3$Response_class,levels = c('PD','SD','PR','CR'))

clinicalall3_ratio=clinicalall3 %>% 
  dplyr::select(Response_class,Score_cut,Cancer) 
 
pvalue_all$Cancer=rownames(pvalue_all)
Rcancerfilter_x=pvalue_all[which(pvalue_all$pvalue_all < 0.05),'Cancer']

clinicalall3_ratio[which(clinicalall3_ratio$Cancer %in% Rcancerfilter_x),'Cancer']=
  paste0(clinicalall3_ratio[which(clinicalall3_ratio$Cancer %in% Rcancerfilter_x),'Cancer'],'*')


pdf(paste0('NRvsR/TRIMscore_category_ratio.pdf'),width = 5,height = 6)
ggplot(clinicalall3_ratio,aes(x=Score_cut))+  
  geom_bar(aes(fill=Response_class),position='fill',width = 0.8)+
  facet_wrap(~Cancer,scales = 'free')+
  theme_bw()+scale_fill_brewer(palette='Set1',direction =-1)+
  theme(legend.position = 'top',strip.background = element_blank(),strip.placement = "outside",strip.switch.pad.grid = unit(1, "inch"))
dev.off()

pvalue_all=c()
Rcancerfilter=sort(Rcancerfilter)
for (i in 1:length(Rcancerfilter)) {
cancer=Rcancerfilter[i]
# print(cancer)
clinicalall3$Response_class=factor(clinicalall3$Response_class,levels = c('PD','SD','PR','CR'))
clinicalall3_ratio=clinicalall3 %>% 
  filter(Cancer == cancer) %>%
  dplyr::select(Response_class,Score_cut) 

tryCatch({
ratio2=as.data.frame(table(clinicalall3_ratio))
a1=ratio2 %>% filter(Score_cut == 'low') %>% dplyr::select(Freq) %>% sum()
a2=ratio2 %>% filter(Score_cut == 'high') %>% dplyr::select(Freq) %>% sum()
ratio2$ratio=1
ratio2[which(ratio2$Score_cut=='low'),'ratio']=ratio2[which(ratio2$Score_cut=='low'),'Freq']/a1
ratio2[which(ratio2$Score_cut=='high'),'ratio']=ratio2[which(ratio2$Score_cut=='high'),'Freq']/a2
ratio2[is.nan(ratio2$ratio),'ratio']=0
p=fisher.test(cbind(ratio2[which(ratio2$Score_cut=='low'),'Freq'],ratio2[which(ratio2$Score_cut=='high'),'Freq']))

pvalue=p$p.value
pvalue_all=c(pvalue_all,pvalue)
},error = function(e){
  print(cancer)
})
}

Rcancerfilter2=Rcancerfilter[which(Rcancerfilter != "CHOL")]
names(pvalue_all)=Rcancerfilter2
pvalue_all=as.data.frame(pvalue_all)

write.csv(pvalue_all,'NRvsR/TRIMscore_category_ratio_pvalue_fisher.test.csv',quote = F)

########################################################################
