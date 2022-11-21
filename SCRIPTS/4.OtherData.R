source('Code/Defined_Functions.R')
work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'



#################Differential expression analysis-additional cohorts#################

cancer_all <- dir('./')

cancer_all_1 <- c('bladder','pancreas','prostate','renal',
                  'ovarian','melanoma','liver','gastric')
cancer_all_2 <- c('breast','colorectal','lung')


## Probe name
probe_name <- read.csv(paste0('bladder/A-AFFY-44.adf.txt'),sep='\t',header = T,check.names = F)
intersect_TRIM <- intersect(Trimdata$ENSG,probe_name[,4])
probe_name_TRIM <- probe_name[which(probe_name[,4] %in% intersect_TRIM),c(1,4)]
colnames(probe_name_TRIM) <- c('Probe','ENSG')
probe_name_TRIM <- merge(probe_name_TRIM,Trimdata[,c(2,5,6)],by = 'ENSG')
#whether one probe corresponds to multiple genes, and delete if so
kk <- table(probe_name_TRIM$Probe)
kk <- kk[kk>1]
if(length(kk) != 0){
   probe_name_TRIM <- probe_name_TRIM[-which(probe_name_TRIM$Probe %in% names(kk)),]
}


k=8# 1:8
cancer0 <- cancer_all_1[k]
cancer0
cancer0_path <- dir(cancer0)
cancer0_path
#Sample_type
sample_path <- grep('.sdrf.txt',cancer0_path,value = T)
clinic <- read.csv(paste0(cancer0,'/',sample_path),sep='\t',header = T,check.names = F)
clinic <- clinic[,c(1,3)]
colnames(clinic) <- c('Sample','Type')
clinic <- clinic[order(clinic$Type),]
unique(clinic$Type)
#"renal" "ovarian":clinic <- clinic[which(clinic$Type %in% c('tumor','normal')),]

## Expression matrix
expr <- read.csv(paste0(cancer0,'/esets_',cancer0,'_exprs.txt'),
                 sep = '\t',header = T,check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- expr[,clinic$Sample]

## Expression TRIM
expr2 <- expr[probe_name_TRIM$Probe,]
expr2$Symbol <- probe_name_TRIM$Symbol
#Multiple probes correspond to one gene, choose their average expression
duplication<-table(expr2$Symbol)
dupl_name <- names(duplication[duplication > 1])
length(dupl_name)
if(length(dupl_name) != 0){
final_all=data.frame()
for(i in 1:length(dupl_name)){
   duplicate=expr2[which(expr2$Symbol == dupl_name[i]),]
   duplicate=apply(duplicate[,-ncol(duplicate)],2,as.numeric)
   final <- t(apply(duplicate,2,mean))
   final<- merge(dupl_name[i],final)
   final_all=rbind(final_all,final)
}
rownames(final_all)=final_all[,1]
final_all2=final_all[,-1]
single_name <- expr2[-which(expr2$Symbol %in% dupl_name),]
rownames(single_name) <- single_name[,'Symbol']
single_name <- single_name[,-ncol(single_name)]
Expression4 <- rbind(single_name,final_all2)
}else{
   Expression4 <- expr2
   rownames(Expression4) <- Expression4$Symbol
   Expression4 <- Expression4[,-ncol(Expression4)]
}
write.table(Expression4,paste0(cancer0,'/',cancer0,'_TRIM_expr.txt'),sep = '\t',quote = F,row.names = T)

## Differential expression analysis
library(limma)
# library(edgeR)
Group <- factor(clinic$Type,levels = c('tumor','normal'))
design <- model.matrix(~0+Group)
colnames(design)=levels(Group)
rownames(design)=clinic$Sample


constrasts = paste(levels(Group),collapse = "-") 
constrasts
# [1] "tumor-normal"
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit <- lmFit(Expression4,design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
DEG$Gene <- rownames(DEG)
DEG$Cancer <- cancer0
DEG <- DEG[order(DEG$logFC,decreasing = T),]
write.table(DEG,paste0(cancer0,'/',cancer0,'_DE_result.txt'),sep = '\t',quote = F,row.names = F)
write.table(DEG,paste0('../Validation_DE_result.txt'),sep = '\t',quote = F,
            row.names = F,col.names = F,append = T)


cancer_all_2 <- c('breast','colorectal','lung')

k=3# 1:3
cancer0 <- cancer_all_2[k]
cancer0
cancer0_path <- dir(cancer0)
cancer0_path
#Sample_type
sample_path <- grep('.sdrf.txt',cancer0_path,value = T)
clinic <- read.csv(paste0(cancer0,'/',sample_path),sep='\t',header = T,check.names = F)
clinic <- clinic[,c(1,3)]
colnames(clinic) <- c('Sample','Type')
clinic <- clinic[order(clinic$Type),]
unique(clinic$Type)
#"colorectal":clinic <- clinic[which(clinic$Type %in% c('tumor','normal')),]

## Expression matrix
expr <- read.csv(paste0(cancer0,'/esets_',cancer0,'_exprs_genes.txt'),
                 sep = '\t',header = T,check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
#"lung":colnames(expr) <- gsub('\\.CEL','',colnames(expr))
# colnames(expr) <- lapply(colnames(expr),function(v){
# a <- unlist(strsplit(v,'_'))[1]
# return(a)
# })
# intersect_Sample <- intersect(clinic$Sample,colnames(expr))
expr <- expr[,clinic$Sample]

intersect_TRIM <- intersect(rownames(expr),Trimdata$Symbol)
Expression4 <- expr[intersect_TRIM,]
write.table(Expression4,paste0(cancer0,'/',cancer0,'_TRIM_expr.txt'),sep = '\t',quote = F,row.names = T)

## Differential expression analysis
library(limma)
Group <- factor(clinic$Type,levels = c('tumor','normal'))
design <- model.matrix(~0+Group)
colnames(design)=levels(Group)
rownames(design)=clinic$Sample


constrasts = paste(levels(Group),collapse = "-") 
constrasts
# [1] "tumor-normal"
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit <- lmFit(Expression4,design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
DEG$Gene <- rownames(DEG)
DEG$Cancer <- cancer0
DEG <- DEG[order(DEG$logFC,decreasing = T),]
write.table(DEG,paste0(cancer0,'/',cancer0,'_DE_result.txt'),sep = '\t',quote = F,row.names = F)
write.table(DEG,paste0('../Validation_DE_result.txt'),sep = '\t',quote = F,
            row.names = F,col.names = F,append = T)


#--------heatmap
library(dplyr)
library(pheatmap)
diff_result_all=read.csv(paste0('../Validation_DE_result.txt'),sep='\t',header=F)
colnames(diff_result_all) <- c('logFC','AveExpr','t',
                               'P.Value','adj.P.Val','B',
                               'Gene','Cancer')
diff_result_all_TRIM=as.data.frame(diff_result_all)
# diff_result_all_TRIM$Cancer=factor(diff_result_all_TRIM$Cancer,levels = cancer_names_analysis)
# diff_result_all_TRIM=diff_result_all_TRIM[order(diff_result_all_TRIM$Cancer),]

diff_result_all_TRIM2_lgfc=diff_result_all_TRIM %>% dplyr::select(Gene,Cancer,logFC)
diff_result_all_TRIM2_lgfc = diff_result_all_TRIM2_lgfc %>% tidyr::spread(Cancer,logFC)
rownames(diff_result_all_TRIM2_lgfc)=diff_result_all_TRIM2_lgfc$Gene
diff_result_all_TRIM2_lgfc=diff_result_all_TRIM2_lgfc[,-1]
diff_result_all_TRIM2_lgfc=diff_result_all_TRIM2_lgfc[Trimdata$Symbol,]
diff_result_all_TRIM2_lgfc=diff_result_all_TRIM2_lgfc[-which(apply(diff_result_all_TRIM2_lgfc,1,function(x){all(is.na(x))})),]

annotation_row3=data.frame(Class=factor(Trimdata$Family))
rownames(annotation_row3)=Trimdata$Symbol

library(RColorBrewer)
mycols<-brewer.pal(12,"Paired")
annotation_colors=list(Class=c('C-I'=mycols[1],'C-II'=mycols[2],'C-III'=mycols[3],
                               'C-IV'=mycols[4],'C-V'=mycols[5],'C-VI'=mycols[6],
                               'C-VII'=mycols[7],'C-VIII'=mycols[8],'C-IX'=mycols[9],
                               'C-X'=mycols[10],'C-XI'=mycols[11],'UC'=mycols[12]))         

diff_result_all_TRIM2_p=diff_result_all_TRIM %>% dplyr::select(Gene,Cancer,P.Value)
diff_result_all_TRIM2_p = diff_result_all_TRIM2_p %>% tidyr::spread(Cancer,P.Value)
rownames(diff_result_all_TRIM2_p)=diff_result_all_TRIM2_p$Gene
diff_result_all_TRIM2_p=diff_result_all_TRIM2_p[,-1]
diff_result_all_TRIM2_p=diff_result_all_TRIM2_p[Trimdata$Symbol,]
diff_result_all_TRIM2_p=diff_result_all_TRIM2_p[-which(apply(diff_result_all_TRIM2_p,1,function(x){all(is.na(x))})),]


data_norm<-diff_result_all_TRIM2_p
data_norm2<-diff_result_all_TRIM2_lgfc
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
bk <- c(seq(-0.5,-0.1,by=0.06),seq(0,0.5,by=0.06))
pdf(paste0('../Validation_heatmap_Log2FC.pdf'),
    width=4,height=10)
p_heatmap =pheatmap(diff_result_all_TRIM2_lgfc,scale = "none",fontsize = 10,breaks=bk,
                    annotation_row=annotation_row3,annotation_colors=annotation_colors,
                    main='log2Fold Change',cluster_rows=FALSE,cluster_cols = FALSE,
                    display_numbers =data_norm,angle_col = "45",
                    fontsize_number = 10,na_col='grey',
                    color=colorRampPalette(rev(c("#CB181D","white","#1F78B4")))(length(bk)),
                    number_color=1)
dev.off()

#---------barplot
library(ggplot2)
dat=diff_result_all_TRIM
dat$Gene=factor(dat$Gene,levels = unique(dat$Gene))
dat_2=dat %>% filter(P.Value < 0.05)

dat_H=dat_2 %>% filter(logFC>0)
hist_H=as.data.frame(table(dat_H$Gene))
colnames(hist_H)=c('Symbol','Positive')

dat_L=dat_2 %>% filter(logFC<0)
hist_L=as.data.frame(table(dat_L$Gene))
colnames(hist_L)=c('Symbol','Negative')

hist=merge(hist_H,hist_L,by='Symbol',all=T)

hist[is.na(hist)]=0

hist2=reshape2::melt(hist)

hist2$Symbol=factor(hist2$Symbol,levels = rev(Trimdata$Symbol ))
hist2=hist2[order(hist2$Symbol),]

result_GSEA_bar=ggplot(hist2, aes(
   x = factor(Symbol,levels = unique(Symbol)),             
   y = ifelse(variable == "Positive", value, -value),  
   fill = variable)) +
   geom_bar(stat = 'identity')+                                
   ylab('Number of Cancer')+
   xlab('TRIM')+
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
   )+coord_flip()

pdf(paste0('../Validation_barplot.pdf'),
    width=2,height=10)

print(result_GSEA_bar)

dev.off()

#------boxplot_expression
library(ggpubr)
diff_result_all_TRIM3=diff_result_all_TRIM %>% dplyr::filter(P.Value < 0.05)
Symbol=unique(diff_result_all_TRIM3$Gene)


cancer_all <- dir('./')
expr_trim_all <- data.frame()
for(m in seq_along(cancer_all)){
   cancer0 <- cancer_all[m]
   #exp
   expr_trim <- read.csv(paste0(cancer0,'/',cancer0,'_TRIM_expr.txt'),
                         sep='\t',header = T,check.names = F)
   expr_trim_2 <- reshape2::melt(as.matrix(expr_trim))
   colnames(expr_trim_2) <- c('Gene','Sample','Value')
   #clinic
   cancer0_path <- dir(cancer0)
   cancer0_path
   #Sample_type
   sample_path <- grep('.sdrf.txt',cancer0_path,value = T)
   clinic <- read.csv(paste0(cancer0,'/',sample_path),sep='\t',header = T,check.names = F)
   clinic <- clinic[,c(1,3)]
   colnames(clinic) <- c('Sample','Type')
   clinic <- clinic[order(clinic$Type),]
   expr_trim_2 <- merge(expr_trim_2,clinic,
                        by='Sample')
   expr_trim_2$Cancer <- cancer0
   expr_trim_all <- rbind(expr_trim_all,expr_trim_2)
}
write.csv(expr_trim_all,
          paste0('../Validation_Expression_allcancer_TRIM_melt.csv'),
          quote = F,row.names = F)

expr_trim_all <- read.csv(
          paste0('../Validation_Expression_allcancer_TRIM_melt.csv'),
          seq = T,header = T)

Symbol <- c('TRIM28','TRIM58','TRIM59','TRIM23')#as.character(unique(expr_trim_all$Gene))

for(k in Symbol){
   # k=Symbol[1]
   up_num=as.numeric(hist2[which(hist2$Symbol == k & hist2$variable == 'Positive'),'value'])
   down_num=as.numeric(hist2[which(hist2$Symbol == k & hist2$variable == 'Negative'),'value'])
   
   
   lg2exp_TRIM_all_1<-as.data.frame(expr_trim_all[which(expr_trim_all$Gene==k),])
   Data_summary <- summarySE(lg2exp_TRIM_all_1, measurevar="Value", groupvars=c('Cancer','Type'))
   colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')
   
   P3 <- ggplot(data=Data_summary, aes(x=Cancer,color=Type)) + 
      geom_errorbar(data = Data_summary,aes(ymin = quantile_25, ymax=quantile_75,group=Type), 
                    width=0.2, 
                    position=position_dodge(0.9), 
                    alpha = 1,
                    size=1) + 
      theme_bw()+
      geom_point(data = Data_summary,aes(x=Cancer, y=median),pch=19,position=position_dodge(0.9),size=6) +
      scale_color_manual(values = c("#1F78B4","#CB181D"))+
      labs(title = paste0(k,' UP(',up_num,') Down(',down_num,')'))+
      ylab(paste('Expression'))+
      theme(
         axis.ticks.x = element_blank(),
         axis.text.y=element_text(size=10), 
         axis.title.y=element_text(size = 12), 
         axis.title.x=element_text(size = 12), 
         plot.title = element_text(hjust = 0.5,size=12,),
         legend.text=element_text(colour="black",  
                                  size=10),
         legend.title=element_text( colour="black", 
                                    size=12),
         panel.grid.major = element_blank(),  
         panel.grid.minor = element_blank()
         
      )+ 
      xlab("Cancer")+
      geom_vline(xintercept=c(seq(1.5,10.5,by=1)), linetype="dotted") 
   
   pdf(paste0('../Validation_boxplot/',k,'_boxplot.pdf'),width=8,height=3)
   print(P3)
   dev.off()
   
}

#####################################################################################






################################# Survival analysis #################################

#--------------GSE176307
GEO='GSE176307'
setwd(paste0(work_path,GEO))
##Expression
Expression=read.csv(paste0(GEO,'_baci_rsem_RS_BACI_headers_tab.txt'),sep='\t',header=T,check.names = F)
head(Expression)
Expression=Expression[-1,]
rownames(Expression)=Expression[,1]
Expression=Expression[,-1]
Expression4 <- Expression
write.table(Expression4,paste0(GEO,'_exp_clear.txt'),quote=F,sep='\t')

##Clinic
library(GEOquery)
gse <- getGEO(GEO,destdir='./')
metdata<-pData(gse[[1]])
head(metdata)
metdata=metdata %>% dplyr::select(
  title,
  'age:ch1',
  'gender:ch1',
  'alive:ch1',
  'overall survival:ch1',
  'd.stage.at.diagnosis:ch1',
  'io.response:ch1'
)
colnames(metdata)=c('Sample','Age','Gender','OS_status','OS','Stage','Response_class')
metdata$Response_status=str_replace_all(metdata$Response_class,
                                        c(
                                          "PD"="NR",
                                          "SD"="NR",
                                          "PR"="R",
                                          "CR"="R"
                                        ))
metdata$Sample=gsub('Patient sample ','',metdata$Sample)
metdata$Sample=gsub('_[0-9]','',metdata$Sample)

metdata$OS_status=str_replace_all(metdata$OS_status,
                                        c(
                                          "No"='0',
                                          "Yes"='1'
                                        ))

IDmap=read.csv('GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv',header = T)

metdata=merge(metdata,IDmap,by.x='Sample',by.y='Sample.ID')
colnames(metdata)[c(1,ncol(metdata))]=c('Sample2','Sample')
metdata=metdata %>% separate_rows(Sample, sep = ",") %>% distinct()
metdata=as.data.frame(metdata)
metdata[,c(2,4,5)]=apply(metdata[,c(2,4,5)], 2, as.numeric)
write.table(metdata,paste0(GEO,'_clinical.txt'),sep='\t',quote = F,row.names = F)

##TRIMs score
metdata=read.table(paste0(GEO,'_clinical.txt'),sep='\t',header=T)
TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'

expdat2=as.matrix(Expression4)
myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Poisson',abs.ranking=TRUE)
myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(metdata,myssgsea2,by='Sample')
write.csv(GEO_clinical_ssgsea,paste0('ssgsea_',GEO,'_clinical.csv'),quote = F,row.names=F)

library("survival")
library("survminer")
data0=GEO_clinical_ssgsea
data0$Score=as.numeric(data0$Score)

res.cut <- surv_cutpoint(data0, 
                         time = "OS", 
                         event = "OS_status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$Sample
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','OS','OS_status'))

##Survival curve
df2=GEO_clinical_ssgsea
survresult=survdiff(Surv(OS, OS_status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value

fit1 <- survfit(Surv(OS, OS_status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
                       data=df2,
                       pval = TRUE, 
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       linetype = "strata",
                       ggtheme = theme_bw(),
                       palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(GEO,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()

#--------------PMID:32472114 Braun et.al
GEO='PMID32472114'
setwd(paste0(work_path,GEO))
##Expression
Expression=read.xlsx(paste0('NIHMS1611472-supplement-Table_S4__RNA_expression__normalized_expression_matrix_.xlsx'),sheet = 1,startRow =2)
head(Expression)
duplication<-table(Expression$gene_name)
duplication <-as.data.frame(duplication)
dupl_name <-duplication[which(duplication[,2]>1),1]
dupl_name <-as.character(dupl_name)
length(dupl_name)

final_all=data.frame()
for(i in 1:length(dupl_name)){
  duplicate=Expression[which(Expression$gene_name == dupl_name[i]),]
  duplicate=apply(duplicate[,-c(1)],2,as.numeric)
  final <- t(apply(duplicate,2,mean))
  final<- merge(dupl_name[i],final)
  final_all=rbind(final_all,final)
}
rownames(final_all)=final_all[,1]
final_all2=final_all[,-1]

single_name <- Expression[-which(Expression$gene_name %in% dupl_name),]
rownames(single_name) <- single_name[,'gene_name']
single_name <- single_name[,-c(1)]

Expression4 <- rbind(single_name,final_all2)

write.table(Expression4,paste0(GEO,'_exp_clear.txt'),quote=F,sep='\t')

##Clinic
metdata<-read.xlsx('NIHMS1611472-supplement-Table_S1__Clinical_and_immune_phenotype_data_for_the_CheckMate_cohorts.xlsx',sheet = 1,startRow = 2)

metdata=metdata[,1:42]
head(metdata)

metdata_1=metdata %>% dplyr::select(
  SUBJID,
  Sex,
  Age,
  OS,
  OS_CNSR,
  PFS,
  PFS_CNSR,
  ORR,
  ImmunoPhenotype
)
colnames(metdata_1)=c('Sample','Sex','Age','OS','OS_status','PFS','PFS_status','Response_class','ImmunoPhenotype')
metdata_2=metdata %>% dplyr::select(
  MAF_Tumor_ID,
  Sex,
  Age,
  OS,
  OS_CNSR,
  PFS,
  PFS_CNSR,
  ORR,
  ImmunoPhenotype
)
colnames(metdata_2)=c('Sample','Sex','Age','OS','OS_status','PFS','PFS_status','Response_class','ImmunoPhenotype')
metdata_2=metdata_2[!is.na(metdata_2$Sample),]
metdata_3=metdata %>% dplyr::select(
  MAF_Normal_ID,
  Sex,
  Age,
  OS,
  OS_CNSR,
  PFS,
  PFS_CNSR,
  ORR,
  ImmunoPhenotype
)
colnames(metdata_3)=c('Sample','Sex','Age','OS','OS_status','PFS','PFS_status','Response_class','ImmunoPhenotype')
metdata_3=metdata_3[!is.na(metdata_3$Sample),]
metdata_4=metdata %>% dplyr::select(
  CNV_ID,
  Sex,
  Age,
  OS,
  OS_CNSR,
  PFS,
  PFS_CNSR,
  ORR,
  ImmunoPhenotype
)
colnames(metdata_4)=c('Sample','Sex','Age','OS','OS_status','PFS','PFS_status','Response_class','ImmunoPhenotype')
metdata_4=metdata_4[!is.na(metdata_4$Sample),]
metdata_5=metdata %>% dplyr::select(
  RNA_ID,
  Sex,
  Age,
  OS,
  OS_CNSR,
  PFS,
  PFS_CNSR,
  ORR,
  ImmunoPhenotype
)
colnames(metdata_5)=c('Sample','Sex','Age','OS','OS_status','PFS','PFS_status','Response_class','ImmunoPhenotype')
metdata_5=metdata_5[!is.na(metdata_5$Sample),]

metdata2=rbind(metdata_1,metdata_2,metdata_3,metdata_4,metdata_5)
metdata2$Response_status=str_replace_all(metdata2$Response_class,
                                        c(
                                          "PD"="NR",
                                          "SD"="NR",
                                          "PR"="R",
                                          "CR"="R"
                                        ))

metdata2=metdata2[which(metdata2$Response_class %in% c('PD','SD','PR','CR')),]
metdata2=as.data.frame(metdata2)

##TRIMs score
TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'

expdat2=as.matrix(Expression4)
myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Gaussian',abs.ranking=TRUE)
myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(metdata2,myssgsea2,by='Sample')


library("survival")
library("survminer")

data0=GEO_clinical_ssgsea
res.cut <- surv_cutpoint(data0, 
                         time = "OS", 
                         event = "OS_status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$Sample
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','OS','OS_status'))

write.table(GEO_clinical_ssgsea,paste0('ssgsea_',GEO,'_clinical.txt'),sep='\t',quote = F,row.names=F)

##Survival curve
df2=GEO_clinical_ssgsea
survresult=survdiff(Surv(OS, OS_status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value

fit1 <- survfit(Surv(OS, OS_status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
                       data=df2,
                       pval = TRUE, 
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       linetype = "strata",
                       ggtheme = theme_bw(),
                       palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(GEO,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()


#-----------------GSE76019
GEO='GSE76019'
setwd(paste0(work_path,GEO))
##Expression
library(affy)
dirmy='data/'
cel.files<-list.files(path=dirmy,pattern='.+\\.CEL$',ignore.case=TRUE,full.names=TRUE,recursive=TRUE)
basename(cel.files)
data.raw=ReadAffy(filenames=cel.files)
eset.rma <- rma(data.raw)
exprSet <- exprs(eset.rma)
colnames(exprSet)=gsub('_(.*).CEL','',colnames(exprSet))
write.table(exprSet, paste0(GEO,"_matrix.txt"), quote = F, sep = "\t")

Expression=readtxt(work_path,paste0(GEO,'/',GEO,'_matrix.txt'),sep0='\t',header0=T)
head(Expression)
library(GEOquery)
GPL=getGEO("GPL13158",destdir = ".")
g1=Table(GPL)
head(g1$ID)

g1=g1 %>% dplyr::select(ID,'Gene Symbol')
colnames(g1)[2]='Gene_Symbol'
g2=g1[grep('///',g1$Gene_Symbol,invert=T),]
g2=g2[-which(g2$Gene_Symbol==''),]

Expression2=merge(Expression,g2,by.x='row.names',by.y='ID')
rownames(Expression2)=Expression2[,1]
Expression3=Expression2[,-1]

duplication<-table(Expression3[,ncol(Expression3)])
duplication <-as.data.frame(duplication)
dupl_name <-duplication[which(duplication[,2]>1),1]
dupl_name <-as.character(dupl_name)
length(dupl_name)

final_all=data.frame()
for(i in 1:length(dupl_name)){
  duplicate=Expression3[which(Expression3$Gene_Symbol == dupl_name[i]),]
  duplicate=apply(duplicate[,-ncol(duplicate)],2,as.numeric)
  final <- t(apply(duplicate,2,mean))
  final<- merge(dupl_name[i],final)
  final_all=rbind(final_all,final)
}
rownames(final_all)=final_all[,1]
final_all2=final_all[,-1]

single_name <- Expression3[-which(Expression3$Gene_Symbol %in% dupl_name),]
rownames(single_name) <- single_name[,ncol(single_name)]
single_name <- single_name[,-ncol(single_name)]

Expression4 <- rbind(single_name,final_all2)
write.table(Expression4,paste0(GEO,'_exp_clear.txt'),quote=F,sep='\t')

##Clinic
library(GEOquery)
gse <- getGEO(GEO,destdir='./')
metdata<-pData(gse[[1]])
head(metdata)
metdata=metdata %>% dplyr::select(
  geo_accession,
  'efs.time:ch1',
  'efs.event:ch1',
  'Stage:ch1'
)
colnames(metdata)=c('Sample','OS','os_status','Stage')
metdata=as.data.frame(metdata)
write.table(metdata,paste0(GEO,'_clinical.txt'),sep='\t',quote = F,row.names = F)

##TRIMs score
TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'

expdat2=as.matrix(Expression4)

myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Gaussian',abs.ranking=TRUE)

myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(metdata,myssgsea2,by='Sample')

library("survival")
library("survminer")
data0=GEO_clinical_ssgsea
data0[,c('OS','os_status')]=apply(data0[,c('OS','os_status')],2,as.numeric)
res.cut <- surv_cutpoint(data0, 
                         time = "OS", 
                         event = "os_status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$Sample
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','OS','os_status'))
write.table(GEO_clinical_ssgsea,paste0(work_path,GEO,'/ssgsea_',GEO,'_clinical.txt'),sep='\t',quote = F,row.names=F)

##Survival curve
df2=GEO_clinical_ssgsea
df2[,c('OS','os_status')]=apply(df2[,c('OS','os_status')],2,as.numeric)

survresult=survdiff(Surv(OS, os_status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value

fit1 <- survfit(Surv(OS, os_status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
                       data=df2,
                       pval = TRUE,
                       risk.table = T,
                       risk.table.col = "strata",
                       linetype = "strata",
                       ggtheme = theme_bw(),
                       title=GEO,
                       palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(GEO,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()

#------------IMvigor210
work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
Geo='IMvigor210'
setwd(paste0(work_path,Geo))
library(pacman)

##Expression
library(IMvigor210CoreBiologies)
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))

eset <- IOBR::anno_eset(eset = eset,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")
tumor_type <- "blca"
if(max(eset)>100) eset <- log2(eset+1)

##Clinic
library(IOBR)
pdata <- pData(cds)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
pdata <- rownames_to_column(pdata[,c("binaryResponse",
									 'Best_Confirmed_Overall_Response',
                                     "FMOne_mutation_burden_per_MB",
                                     "Neoantigen_burden_per_MB",
                                     "censOS","os")],var = "ID")
colnames(pdata)<-c("ID","BOR_binary","Overall_Response","TMB","TNB","status","time")
pdata<-pdata[!is.na(pdata$BOR_binary),]
pdata$BOR_binary<-ifelse(pdata$BOR_binary=="CR/PR","R","NR")
save(eset,pdata,file = "expcli_IMvigor210.Rdata")

##TRIMs score
TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'
expdat2=as.matrix(eset)

myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Gaussian',abs.ranking=TRUE)

myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(pdata,myssgsea2,by.y='Sample',by.x='ID')
library("survival")
library("survminer")

data0=GEO_clinical_ssgsea
res.cut <- surv_cutpoint(data0, 
                         time = "time", 
                         event = "status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=GEO_clinical_ssgsea$ID
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','time','status'))

write.table(GEO_clinical_ssgsea,paste0(work_path,Geo,'/ssgsea_',Geo,'_clinical.txt'),sep='\t',quote = F,row.names=F)

##Survival curve
df2=GEO_clinical_ssgsea
survresult=survdiff(Surv(time, status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value

fit1 <- survfit(Surv(time, status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
           data=df2,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(Geo,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()



#Add
cox_result_score=read.csv('../Surv/Survival_HR.txt',sep = '\t',
                          header = F)
colnames(cox_result_score) <- c('Pvalue','HR','Up95','Low95','Cancer')
cox_result_score$Type='Sig'
cox_result_score[which(cox_result_score$Pvalue > 0.05),'Type']='Not-Sig'
cox_result_score2=cox_result_score
cox_result_score2[,2:4]=log2(cox_result_score2[,2:4])
P3 <- ggplot(data=cox_result_score2, aes(x=Cancer,color=Type)) + 
   geom_errorbar(data = cox_result_score2,aes(ymin = round(Low95,4), ymax=round(Up95,4)),
                 width=0.2, 
                 position=position_dodge(0.9), 
                 alpha = 1,
                 size=1) + 
   theme_bw()+
   
   geom_point(data = cox_result_score2,aes(x=Cancer, y=round(HR,4)),pch=19,position=position_dodge(0.9),size=5) +
   scale_color_manual(values = c('#999999','#D6604D') )+
   theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank()
   )+  
   xlab("Cancer")+ylab("log2 TRIM score HR(95% CI)")+ 
   geom_hline(yintercept=0, linetype="dotted")
P3+coord_flip()+annotate("text",x=2,y=-2,label=paste(signif(cox_result_score$Pvalue,digits = 3),collapse = '\n'))

pdf('../Surv/Validation_TRIM_score_COX.pdf',width = 5,height = 6)
P3+coord_flip()+annotate("text",x=2,y=-2,label=paste(signif(cox_result_score$Pvalue,digits = 3),collapse = '\n'))
dev.off()


######################################################################################################


################################# Drug response #################################
#------------GSE78220
```Linux

##download fastq data
wget -c -i fastq_ftp.csv
md5sum *.fastq.gz > file_sum.md5
md5sum -c file_sum.md5

work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
cd ${work_path}
Geo='GSE78220'
cd ${Geo}

##loading software path
export PATH=/boot3/share/software/fastp:$PATH
export PATH=/boot3/share/software:$PATH
ls data/*fastq.gz |cut -d"/" -f 2 |cut -d"." -f 1|grep '_1$'|cut -c 1-10 | sort -u|xargs -I {} fastp -i data/{}_1.fastq.gz -I data/{}_2.fastq.gz -o outdata/out_{}_1.fastq.gz -O outdata/out_{}_2.fastq.gz -w 8 --html outdata/out_{}_paired.html --json outdata/out_{}_paired.json
genomeDir0='/boot2/bio_gaoyueying/study/three/RNA-seq'

##Build Index
STAR --runThreadN 10 --runMode genomeGenerate \
--genomeDir ${work_path}/hg19 \
--genomeFastaFiles ${genomeDir0}/GRCh37.p13.genome.fa \
--sjdbGTFfile ${genomeDir0}/gencode.v19.annotation.gtf \
--sjdbOverhang 149 

##STAR
ls outdata/*fastq.gz | cut -d"/" -f 2 | cut -d"." -f 1 | grep '_1$' | cut -c 1-14 | sort -u |
xargs -I {} echo 'STAR --runThreadN 10 --genomeDir '${work_path}hg19' --readFilesCommand zcat --readFilesIn 'outdata/{}_1.fastq.gz outdata/{}_2.fastq.gz' --outFileNamePrefix 'STARresult/{}_paired/{}_paired_ '--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --quantMode TranscriptomeSAM GeneCounts' > ${Geo}_STARalign_paired.sh

nohup sh ${Geo}_STARalign_paired.sh > ${Geo}_STARalign_paired.log 2>&1 &

##HTseq
ls outdata/*fastq.gz | cut -d"/" -f 2 | cut -d"." -f 1 | grep '_1$' | cut -c 1-14 | sort -u |
xargs -I {} echo 'htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m union' \
STARresult/{}_paired/{}_paired_Aligned.sortedByCoord.out.bam \
${genomeDir0}'/gencode.v19.annotation.gtf > 'STARresult/{}_paired/{}_paired_'counts.txt' > ${work_path}${Geo}/${Geo}_HTSEQ_paired.sh

nohup sh ${Geo}_HTSEQ_paired.sh > ${Geo}_HTSEQ_paired.log 2>&1 &

```

##Merge sample expression matrix
work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
Geo='GSE78220'

Dir1=paste0(work_path,Geo,'/STARresult')
File0="counts.txt"
myfiles <- list.files(path=Dir1,pattern =File0 ,all.files=T,full.names=T,recursive=T)

library(dplyr)
myfiles_SRR=lapply(myfiles,function(a){
a1=unlist(strsplit(a,'/'))[10]
b=unlist(strsplit(a1,'_'))[2]
return(b)
} 
)%>% unlist()

mergedata<-read.csv(myfiles[1],sep='\t',header=F)
colnames(mergedata)=c('ENSG',myfiles_SRR[1])
for(i in 2:length(myfiles)){
   mergedata0=read.table(myfiles[i],sep='\t',header=F)
   colnames(mergedata0)=c('ENSG',myfiles_SRR[i])
   mergedata = merge(mergedata,mergedata0,by='ENSG')
}
mergedata$ENSG=gsub('\\.[0-9]*','',mergedata$ENSG)

write.table(mergedata,paste0(Geo,'_count.txt'),row.names=F,quote=F,sep='\t')

##Clinic
samplemy=colnames(mergedata)[-1]
clinical_1=read.csv('GSE78220_clinical.txt',sep=',',header=T) %>% filter(Run %in% samplemy)

library(GEOquery)
gse <- getGEO("GSE78220",destdir='./') 
metdata<-pData(gse[[1]])

metdata=metdata %>% dplyr::select(title,geo_accession,
 characteristics_ch1.4,characteristics_ch1.6)
colnames(metdata)=c('Patient','Sample.Name','Age','OS')

clinical_2=merge(clinical_1,metdata,by='Sample.Name')

clinical_2$Age=lapply(clinical_2$Age,function(a){
a1=unlist(strsplit(a,': '))[2]
b=as.numeric(a1)
return(b)
}) %>% unlist()

clinical_2$OS=lapply(clinical_2$OS,function(a){
a1=unlist(strsplit(a,': '))[2]
b=as.numeric(a1)
return(b)
}) %>% unlist()

library(stringr)
clinical_2$Response=str_replace_all(clinical_2$anti.pd.1_response, 
           c("Progressive Disease" = "PD",
		   "Partial Response" = "PR",
		   "Complete Response" = "CR"
		   ))

clinical_2$Response_status=str_replace_all(clinical_2$Response, 
           c("PD" = "nonresponse",
		   "PR" = "response",
		   "CR" = "response"
		   ))
		   
clinical_2$os_status=str_replace_all(clinical_2$vital_status, 
           c("Dead" = '1',
		   "Alive" = '0'
		   ))
clinical_2$os_status=as.numeric(clinical_2$os_status)

write.table(clinical_2,'GSE78220_clinical_all.txt',sep='\t',quote=F,row.names=F)

##Expression

work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
Geo='GSE78220'
setwd(paste0(work_path,Geo))

Expression=readtxt(work_path,paste0(Geo,'/',Geo,'_count.txt'),sep0='\t',header0=T)
head(Expression)

rownames(Expression)=Expression$ENSG
Expression=Expression %>% dplyr::select(-ENSG)

ENSG=unique(rownames(Expression))
library(org.Hs.eg.db)

ensembls <- mapIds(org.Hs.eg.db, keys = ENSG, keytype = "ENSEMBL", column="SYMBOL")
ensembls=as.matrix(ensembls)
ensembls=cbind(ensembls,rownames(ensembls))
colnames(ensembls)=c('Hugo_Symbol','ENSG')
ensembls=apply(ensembls, 2, as.character)
ensembls=as.data.frame(ensembls)
ensembls2=ensembls[!is.na(ensembls$Hugo_Symbol),]

Expression2=merge(ensembls2,Expression,by.x='ENSG',by.y='row.names')

duplication<-table(Expression2$Hugo_Symbol)
duplication <-as.data.frame(duplication)
dupl_name <-duplication[which(duplication[,2]>1),1]
dupl_name <-as.character(dupl_name)

final_all=data.frame()
for(i in 1:length(dupl_name)){
  duplicate=Expression2[which(Expression2$Hugo_Symbol == dupl_name[i]),]
  duplicate=apply(duplicate[,-c(1,2)],2,as.numeric)
  final <- t(apply(duplicate,2,mean))
  final<- merge(dupl_name[i],final)
  final_all=rbind(final_all,final)

}
rownames(final_all)=final_all[,1]
final_all=final_all[,-1]

single_name <- Expression2[-which(Expression2$Hugo_Symbol %in% dupl_name),]
rownames(single_name) <- single_name$Hugo_Symbol
single_name <- single_name[,-c(1,2)]
Expression3 <- rbind(single_name,final_all)

write.table(Expression3,paste0(work_path,Geo,'/',Geo,'_count_clear.txt'),quote=F,sep='\t',row.names=T)

##TRIMs score
GEO_clinical=readtxt_1('./','GSE78220_clinical_all.txt',sep0='\t',header0=T)
TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'
Expression3=read.csv(paste0(work_path,Geo,'/',Geo,'_count_clear.txt'),sep='\t',header=T)

expdat2=as.matrix(Expression3)

myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Poisson',abs.ranking=TRUE)

myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(GEO_clinical,myssgsea2,by.y='Sample',by.x='Run')
head(GEO_clinical_ssgsea)

library("survival")
library("survminer")
data0=GEO_clinical_ssgsea
res.cut <- surv_cutpoint(data0, 
                         time = "OS", 
                         event = "os_status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$Run
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','OS','os_status'))
write.table(GEO_clinical_ssgsea,paste0(work_path,Geo,'/ssgsea_',Geo,'_clinical.txt'),sep='\t',quote = F,row.names=F)

##boxplot
GEO_clinical_ssgsea$Response=factor(GEO_clinical_ssgsea$Response,levels=c('PD','SD','PR','CR'))
GEO_clinical_ssgsea=GEO_clinical_ssgsea[order(GEO_clinical_ssgsea$Response),]

compre=unique(GEO_clinical_ssgsea$Response)

compre=list(
c('PD','PR'),
c('PD','CR'),
c('PR','CR'))

Data_summary <- summarySE(GEO_clinical_ssgsea, measurevar="Score", groupvars=c('Response'))
colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')

p=ggplot(data=GEO_clinical_ssgsea,aes(x=Response,y=Score,color=Response))+
   stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
   geom_boxplot(position = position_dodge(0.9),width=0.6)+
   geom_point(data = Data_summary,aes(x=Response,y=median,color=Response),alpha = 1,pch = 19,position = position_dodge(0.1),size = 3)+
   geom_line(data = Data_summary,aes(x=Response,y=median,group=1) ,size=1,colour=brewer.pal(4, 'Set1')[3:1])+
   geom_jitter(aes(color=Response),size=1)+
   ylab('TRIM score')+xlab('response_category')+labs(title = Geo)+
   scale_color_brewer(palette='Set1',direction =-1)+
   theme_bw()+
   theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.y = element_text(size = 10),
         legend.position = 'none',
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(color = 'black',
                                   size=10,
                                   hjust=0.5))

pdf(paste0(Geo,'_TRIMscore_category.pdf'),width=4,height=6)
p+stat_compare_means(comparisons = compre)
dev.off()

##barplot
GEO_clinical_ssgsea_ratio=GEO_clinical_ssgsea %>% dplyr::select(Response,Score_cut)

ratio2=as.data.frame(table(GEO_clinical_ssgsea_ratio))
a1=ratio2 %>% filter(Score_cut == 'low') %>% dplyr::select(Freq) %>% sum()
a2=ratio2 %>% filter(Score_cut == 'high') %>% dplyr::select(Freq) %>% sum()

ratio2$ratio=1
ratio2[which(ratio2$Score_cut=='low'),'ratio']=ratio2[which(ratio2$Score_cut=='low'),'Freq']/a1
ratio2[which(ratio2$Score_cut=='high'),'ratio']=ratio2[which(ratio2$Score_cut=='high'),'Freq']/a2

write.table(ratio2,paste0(Geo,'_TRIMscore_category_ratio.txt'),sep='\t',quote=F,row.names=F)

fishresult=fisher.test(cbind(ratio2[which(ratio2$Score_cut=='low'),'Freq'],ratio2[which(ratio2$Score_cut=='high'),'Freq']))
fishresult_p=signif(fishresult$p.value,2)

pdf(paste0(Geo,'_TRIMscore_category_ratio.pdf'),width=4,height=6)
ggplot(GEO_clinical_ssgsea_ratio,aes(x=Score_cut))+  
  geom_bar(aes(fill=Response),position='fill',width = 0.8)+
  theme_bw()+labs(title = Geo)+scale_fill_brewer(palette='Set1',direction =-1)+
  geom_text(aes(x=1.5,y=1.02),label=paste0('p=',fishresult_p))
dev.off()

##survival curve
df2=GEO_clinical_ssgsea
survresult=survdiff(Surv(OS, os_status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value
fit1 <- survfit(Surv(OS, os_status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
           data=df2,
           pval = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = "strata", 
           ggtheme = theme_bw(), 
           title=Geo,
           palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(Geo,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()

#------------PMID:27956380
##Clinic
clinical_1=read.xlsx('cohort.xlsx')
library(stringr)
clinical_1$Response=str_replace_all(clinical_1$Benefit, 
                                    c("TRUE" = "R",
                                      "FALSE" = "NR"
                                    ))

clinical_1$os_status=str_replace_all(clinical_1$Alive, 
                                     c("1" = 'alive',
                                       "0" = 'dead'
                                     ))
clinical_1$os_status=str_replace_all(clinical_1$os_status, 
                                     c("alive" = '0',
                                       "dead" = '1'
                                     ))
clinical_1$os_status=as.numeric(clinical_1$os_status)
write.table(clinical_1,'PMID27956380_clinical_all.txt',sep='\t',quote=F,row.names=F)

##Expression
work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
Geo='PMID27956380'
setwd(paste0(work_path,Geo))

Expression=readtxt_1(filename=paste0('cufflinks.csv'),sep0=',',header0=T)
head(Expression)
Expression2=Expression %>% dplyr::select(gene_short_name,sample,FPKM)
Expression2=as.data.frame(Expression2)
Expression2=Expression2 %>%
  group_by(gene_short_name,sample) %>%
  summarise(FPKM_mean = mean(FPKM))

Expression3=Expression2 %>% spread(sample,FPKM_mean)
Expression3=as.data.frame(Expression3)
rownames(Expression3)=Expression3$gene_short_name
Expression3=Expression3 %>% dplyr::select(-gene_short_name)
write.table(Expression3,paste0(Geo,'_count_clear.txt'),quote=F,sep='\t',row.names=T)

##TRIMs score
Expression3=read.csv(paste0(Geo,'_count_clear.txt'),sep='\t',header = T,check.names = F)
GEO_clinical=readtxt_1('./','PMID27956380_clinical_all.txt',sep0='\t',header0=T)

TRIM=list(Trimdata$Symbol)
names(TRIM)='TRIM'

expdat2=as.matrix(log2(Expression3+0.05))
myssgsea<-gsva(expdat2,TRIM,method='ssgsea',min.sz =1.01,kcdf='Gaussian',abs.ranking=TRUE)

myssgsea2=reshape2::melt(myssgsea)
colnames(myssgsea2)=c('Immune_cell','Sample','Score')
GEO_clinical_ssgsea=merge(GEO_clinical,myssgsea2,by.y='Sample',by.x='sample')

library("survival")
library("survminer")

data0=GEO_clinical_ssgsea
res.cut <- surv_cutpoint(data0, 
                         time = "OS", 
                         event = "os_status",
                         variables = 'Score')

res.cat <- surv_categorize(res.cut)

res.cat[,'Score']<-factor(res.cat[,'Score'],levels=c("low","high"))
colnames(res.cat)[3]='Score_cut'
res.cat$Sample=data0$sample
colnames(GEO_clinical_ssgsea)[1]='Sample'
GEO_clinical_ssgsea=merge(GEO_clinical_ssgsea,res.cat,by=c('Sample','OS','os_status'))
write.table(GEO_clinical_ssgsea,paste0('ssgsea_',Geo,'_clinical.txt'),sep='\t',quote = F,row.names=F)

##boxplot
colmy=c("#CB181D","#1F78B4")
names(colmy)<-c('NR','R')

p=ggplot(data=GEO_clinical_ssgsea,aes(x=Response,y=Score,color=Response))+
  stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
  geom_boxplot(position = position_dodge(0.9),width=0.6)+
  geom_jitter(aes(color=Response),size=1)+
  ylab('TRIM score')+xlab('response_status')+labs(title = Geo)+
  scale_color_manual(values = colmy)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = 'black',
                                  size=10,
                                  hjust=0.5))

pdf(paste0(Geo,'_TRIMscore_NRvsR.pdf'),width=3,height=6)
p+stat_compare_means()
dev.off()

##barplot
GEO_clinical_ssgsea_ratio=GEO_clinical_ssgsea %>% dplyr::select(Response,Score_cut)

ratio2=as.data.frame(table(GEO_clinical_ssgsea_ratio))
a1=ratio2 %>% filter(Score_cut == 'low') %>% dplyr::select(Freq) %>% sum()
a2=ratio2 %>% filter(Score_cut == 'high') %>% dplyr::select(Freq) %>% sum()

ratio2$ratio=1
ratio2[which(ratio2$Score_cut=='low'),'ratio']=ratio2[which(ratio2$Score_cut=='low'),'Freq']/a1
ratio2[which(ratio2$Score_cut=='high'),'ratio']=ratio2[which(ratio2$Score_cut=='high'),'Freq']/a2

write.table(ratio2,paste0(Geo,'_TRIMscore_category_ratio.txt'),sep='\t',quote=F,row.names=F)


fishresult=fisher.test(cbind(ratio2[which(ratio2$Score_cut=='low'),'Freq'],ratio2[which(ratio2$Score_cut=='high'),'Freq']))
fishresult_p=signif(fishresult$p.value,2)

pdf(paste0(Geo,'_TRIMscore_category_ratio.pdf'),width=4,height=6)
ggplot(GEO_clinical_ssgsea_ratio,aes(x=Score_cut))+  
  geom_bar(aes(fill=Response),position='fill',width = 0.8)+
  labs(title = Geo)+scale_fill_brewer(palette='Set1',direction =-1)+
  geom_text(aes(x=1.5,y=1.02),label=paste0('p=',fishresult_p))+
  theme_bw()
dev.off()

##survival curve
df2=GEO_clinical_ssgsea
survresult=survdiff(Surv(OS, os_status) ~ Score_cut, data = df2)
p.value <- 1 - pchisq(survresult$chisq, length(survresult$n) -1)
p.value

fit1 <- survfit(Surv(OS, os_status) ~ Score_cut, data = df2)
mysurvplot<-ggsurvplot(fit1,
                       data=df2,
                       pval = TRUE, 
                       risk.table = TRUE, 
                       risk.table.col = "strata",
                       linetype = "strata", 
                       ggtheme = theme_bw(), 
                       title = Geo,
                       palette = c('Score_cut=high'="#CB181D",'Score_cut=low'="#1F78B4")
)

pdf(paste0(Geo,'_TRIMscore_LH_survival.pdf'),width=8,height=8,onefile = FALSE)
print(mysurvplot)
dev.off()

######################################################################################################

################################# ROC #################################
work_path='/boot2/bio_gaoyueying/study/one2020/Data4/other_data/'
setwd(paste0(work_path))

Geo1='GSE78220'
Geo2='PMID27956380'

df1=read.csv(paste0(Geo1,'/ssgsea_',Geo1,'_clinical.txt'),sep='\t',header=T)
df2=read.csv(paste0(Geo2,'/ssgsea_',Geo2,'_clinical.txt'),sep='\t',header=T)

pdf('merge_ROC.pdf')
ROC.scoreTN_legend=c()

library(pROC)  
proc<-roc(df1$Response_status,df1$Score,ci=T,smooth =TRUE);proc$auc
ROC.scoreTN_legend=c(ROC.scoreTN_legend,paste0(Geo1,": ",round(as.numeric(proc$auc),2)))
plot(proc, col=brewer.pal(9,"Set1")[1], lwd=2, title = "ROC for TRIMscore")

proc<-roc(df2$Response,df2$Score,ci=T,smooth =TRUE);proc$auc
ROC.scoreTN_legend=c(ROC.scoreTN_legend,paste0("Tavi et al.,2017: ",round(as.numeric(proc$auc),2)))
plot(proc, col=brewer.pal(9,"Set1")[3], lwd=2,add = T )

legend("bottomright",
       ROC.scoreTN_legend,
       col=brewer.pal(9,"Set1")[c(1,3)],
       lty=1, lwd=2,bty = "n")

dev.off()

######################################################################################################
