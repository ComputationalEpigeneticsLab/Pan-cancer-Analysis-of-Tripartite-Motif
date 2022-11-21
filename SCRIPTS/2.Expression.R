remove(list = ls())
source('Code/Defined_Functions.R')
setwd('expression/DE')


#################### Differential expression analysis ####################

cancer_names_analysis=readtxt_1(workdir3,filename = 'cancer_names_analysis.txt',header0 = F)
cancer_names_analysis=cancer_names_analysis[,1]

diff_result_all<-data.frame()
for(i in 1:length(cancer_names_analysis)){
cancer_names_analysis0=cancer_names_analysis[i]
cancerpath=paste0(workdir3,'expression/expression_clear/')

cancerdata=paste0(cancer_names_analysis0,'_cancer.txt')
expdat_cancer=read.csv(paste0(cancerpath,cancerdata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
expdat_cancer=expdat_cancer %>% dplyr::select(-Hugo_Symbol)
expdat_cancer_TRIM=expdat_cancer[intersect(rownames(expdat_cancer),Trimdata$ENSG),]

normaldata=paste0(cancer_names_analysis0,'_normal.txt')
expdat_normal=read.csv(paste0(cancerpath,normaldata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
expdat_normal=expdat_normal %>% dplyr::select(-Hugo_Symbol)
expdat_normal_TRIM=expdat_normal[intersect(rownames(expdat_normal),Trimdata$ENSG),]

geneexpdata=merge(expdat_cancer_TRIM,expdat_normal_TRIM,by='row.names');rownames(geneexpdata)=geneexpdata$Row.names;geneexpdata=geneexpdata[,-1]
cleargene<-rownames(geneexpdata)[which(apply(geneexpdata,1,function(v){return((sum(v==0)/length(v))>=0.7)}))]

remove(geneexpdata)
if(length(cleargene)!=0){
  expdat_cancer_TRIM=expdat_cancer_TRIM[-which(rownames(expdat_cancer_TRIM) %in% cleargene),]
  expdat_normal_TRIM=expdat_normal_TRIM[-which(rownames(expdat_normal_TRIM) %in% cleargene),]
}

expdat_normal_TRIM=expdat_normal_TRIM[rownames(expdat_cancer_TRIM),]

cancer<-rowMeans(expdat_cancer_TRIM)
normal<-rowMeans(expdat_normal_TRIM)
foldchange<- cancer/normal
foldchange2<-as.data.frame(foldchange)
foldchangelog2<-log2(foldchange)

lg2cancerexp<-log2(expdat_cancer_TRIM+1) 
lg2normalexp<-log2(expdat_normal_TRIM+1) 

list0<-cbind(c(1:nrow(lg2cancerexp)),c(1:nrow(lg2normalexp)))
wilcox_test<-function(list0){
  cancer<-as.numeric(lg2cancerexp[list0[1],])
  normal<-as.numeric(lg2normalexp[list0[2],])
  wtestvalue<-wilcox.test(cancer, normal)
  pvalue0<-wtestvalue$p.value
  return(pvalue0)
}
pvalue<-apply(list0,1,wilcox_test)
names(pvalue)<-rownames(lg2cancerexp)
pvalue_adjust<-p.adjust(pvalue,"fdr")
diff_result<-cbind(foldchange2,foldchangelog2,pvalue,pvalue_adjust)
diff_result$ENSG<-rownames(lg2cancerexp)
diff_result$Cancer=cancer_names_analysis0
diff_result_all=rbind(diff_result_all,diff_result)
cat(cancer_names_analysis0,sep = "\n")
}

Trimdata2=Trimdata[,c('Symbol','ENSG')]
diff_result_all2 = merge(Trimdata2,diff_result_all,by='ENSG')
write.csv(diff_result_all2,paste0('Diff_analysis_allcancer.csv'),quote = F,row.names = F)

#---------- heatmap
diff_result_all=readtxt_1(filename = paste0('Diff_analysis_allcancer.csv'),sep0 = ',')
diff_result_all_TRIM=as.data.frame(diff_result_all)
diff_result_all_TRIM$Cancer=factor(diff_result_all_TRIM$Cancer,levels = cancer_names_analysis)
diff_result_all_TRIM=diff_result_all_TRIM[order(diff_result_all_TRIM$Cancer),]

diff_result_all_TRIM2_lgfc=diff_result_all_TRIM %>% dplyr::select(Symbol,Cancer,foldchangelog2)
diff_result_all_TRIM2_lgfc = diff_result_all_TRIM2_lgfc %>% tidyr::spread(Cancer,foldchangelog2)
rownames(diff_result_all_TRIM2_lgfc)=diff_result_all_TRIM2_lgfc$Symbol
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

diff_result_all_TRIM2_p=diff_result_all_TRIM %>% dplyr::select(Symbol,Cancer,pvalue_adjust)
diff_result_all_TRIM2_p = diff_result_all_TRIM2_p %>% spread(Cancer,pvalue_adjust)
rownames(diff_result_all_TRIM2_p)=diff_result_all_TRIM2_p$Symbol
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
bk <- c(seq(-3,-0.1,by=0.06),seq(0,3,by=0.06))
pdf(paste0('heatmap_Log2FC.pdf'),width=10,height=20)
p_heatmap =pheatmap(diff_result_all_TRIM2_lgfc,scale = "none",fontsize = 10,breaks=bk,
                              annotation_row=annotation_row3,annotation_colors=annotation_colors,
                              main='log2Fold Change',cluster_rows=TRUE,cluster_cols = FALSE,
                              display_numbers =data_norm,angle_col = "45",
                              fontsize_number = 10,na_col='grey',
                              color=colorRampPalette(rev(c("#CB181D","white","#1F78B4")))(102),
                              number_color=1,cellheight=15,cellwidth=20)
dev.off()

#---------- barplot
dat=diff_result_all_TRIM
dat$Symbol=factor(dat$Symbol,levels = unique(dat$Symbol))
dat_2=dat %>% filter(pvalue_adjust < 0.05)

dat_H=dat_2 %>% filter(foldchangelog2>0)
hist_H=as.data.frame(table(dat_H$Symbol))
colnames(hist_H)=c('Symbol','Positive')

dat_L=dat_2 %>% filter(foldchangelog2<0)
hist_L=as.data.frame(table(dat_L$Symbol))
colnames(hist_L)=c('Symbol','Negative')

hist=merge(hist_H,hist_L,by='Symbol',all=T)

hist[is.na(hist)]=0

hist2=melt(hist)

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
  )

result_GSEA_bar=result_GSEA_bar+coord_flip()+
  geom_text(                                                  
    aes(label=value,                                          
        hjust = ifelse(variable == "Positive", -0.4, 1.1)   
    ),
    size=2                                                   
    
  )
pdf('Gene_ggplot_bar.pdf',width = 2.5, height = 8)
print(result_GSEA_bar)
dev.off()


#---------- boxplot
diff_result_all_TRIM3=diff_result_all_TRIM %>% dplyr::filter(pvalue_adjust < 0.05)
ENSG=unique(diff_result_all_TRIM3$ENSG)
Symbol=sort(unique(diff_result_all_TRIM3$Symbol))

dat=diff_result_all_TRIM
dat_2=dat %>% filter(pvalue_adjust < 0.05)

dat_H=dat_2 %>% filter(foldchangelog2>0)
hist_H=as.data.frame(table(dat_H$Symbol))
colnames(hist_H)=c('Symbol','Positive')

dat_L=dat_2 %>% filter(foldchangelog2<0)
hist_L=as.data.frame(table(dat_L$Symbol))
colnames(hist_L)=c('Symbol','Negative')

hist=merge(hist_H,hist_L,by='Symbol',all=T)
hist[is.na(hist)]=0

hist2=melt(hist)
hist2$Symbol=factor(hist2$Symbol,levels = rev(Trimdata$Symbol ))
hist2=hist2[order(hist2$Symbol),]


for(k in Symbol){
  up_num=as.numeric(hist2[which(hist2$Symbol == k & hist2$variable == 'Positive'),'value'])
  down_num=as.numeric(hist2[which(hist2$Symbol == k & hist2$variable == 'Negative'),'value'])
  
  
  lg2exp_TRIM_all_1<-as.data.frame(lg2exp_TRIM_all[which(lg2exp_TRIM_all$TRIM==k),])
  Data_summary <- summarySE(lg2exp_TRIM_all_1, measurevar="Count", groupvars=c('Cancer','Sampletype'))
  colnames(Data_summary)[7:8]<-c('quantile_25','quantile_75')
  
  P3 <- ggplot(data=Data_summary, aes(x=Cancer,color=Sampletype)) + 
    geom_errorbar(data = Data_summary,aes(ymin = quantile_25, ymax=quantile_75,group=Sampletype),
                  width=0.2, 
                  position=position_dodge(0.9), 
                  alpha = 1,
                  size=1) + 
    theme_bw()+
    
    geom_point(data = Data_summary,aes(x=Cancer, y=median),pch=19,position=position_dodge(0.9),size=6) +
    scale_color_manual(values = c("#1F78B4","#CB181D"))+ 
    
    labs(title = paste0(k,' UP(',up_num,') Down(',down_num,')'))+
    ylab(paste(k,'log2(count+1)'))+
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
    geom_vline(xintercept=c(seq(1.5,17.5,by=1)), linetype="dotted") +

  pdf(paste0('DE_boxplot/',k,'_boxplot_2.pdf'),width=12,height=5)
  print(P3)
  dev.off()
  
}



###############################################################################################





#################### correlation analysis ####################
setwd('expression/Cor')
library(corrplot)

spearman_result_all<-data.frame()
for(i in 1:length(cancer_names)){
  # i=1
  cancer_names_analysis0=cancer_names[i]
  cancerpath=paste0(workdir3,'expression/expression_clear/')
  
  cancerdata=paste0(cancer_names_analysis0,'_cancer.txt')
  expdat_cancer=read.csv(paste0(cancerpath,cancerdata),sep='\t',header=T,stringsAsFactors=F,row.names = 1,check.names = F)
  expdat_cancer=expdat_cancer %>% dplyr::select(-Hugo_Symbol)

  expdat_cancer_TRIM=expdat_cancer[intersect(rownames(expdat_cancer),Trimdata$ENSG),]
  
  cleargene<-rownames(expdat_cancer_TRIM)[which(apply(expdat_cancer_TRIM,1,function(v){return((sum(v==0)/length(v))==1)}))]
  if(length(cleargene)>0){
  expdat_cancer_TRIM=expdat_cancer_TRIM[-which(rownames(expdat_cancer_TRIM) %in% cleargene),]
  }
  lg2cancerexp<-log2(expdat_cancer_TRIM+1) 
  lg2cancerexp_need=t(lg2cancerexp)
    
  myresult=rcorr(lg2cancerexp_need, type='spearman')
  #extract correlation results
  myresult_r=myresult$r
  myresult_P=myresult$P
  myresult_r2=reshape2::melt(as.matrix(myresult_r))
  myresult_P2=reshape2::melt(as.matrix(myresult_P))
  colnames(myresult_r2)[3]='R'
  colnames(myresult_P2)[3]='Pvalue'
  myresult_all=merge(myresult_r2,myresult_P2,by=c('Var1','Var2'))
  colnames(myresult_all)[1:2]=c('TRIM_1','TRIM_2')
  myresult_all$Padjust=p.adjust(myresult_all$Pvalue,method = 'fdr')
  myresult_all$Cancer=cancer_names_analysis0
  spearman_result_all<-rbind(spearman_result_all,myresult_all)
}

Trimdata2=Trimdata %>% dplyr::select(Symbol,ENSG,Family)
spearman_result_all2=merge(spearman_result_all,Trimdata2,by.x='TRIM_1',by.y='ENSG')
colnames(spearman_result_all2)[(ncol(spearman_result_all2)-1):ncol(spearman_result_all2)]=c('Symbol_1','Family_1')
spearman_result_all2=merge(spearman_result_all2,Trimdata2,by.x='TRIM_2',by.y='ENSG')
colnames(spearman_result_all2)[(ncol(spearman_result_all2)-1):ncol(spearman_result_all2)]=c('Symbol_2','Family_2')

spearman_result_all2$Family_Type='Different'
spearman_result_all2[which(spearman_result_all2$Family_1 == spearman_result_all2$Family_2),'Family_Type']='Same'

write.csv(spearman_result_all2,paste0('Spearman_TRIM.csv'),row.names = F,quote = F)


spearman_result_all2_sig=dplyr::filter(spearman_result_all2,Padjust < 0.05)
write.csv(spearman_result_all2_sig,paste0('Spearman_TRIM_sig.csv'),row.names = F,quote = F)


#---------- boxplot
spearman_result_all2=readtxt_1(filename=paste0('Spearman_TRIM.csv'),sep0=',')
spearman_result_all2=as.data.frame(spearman_result_all2)
spearman_result_all2[,c(3,4,5)]=apply(spearman_result_all2[,c(3,4,5)],2,as.numeric)
spearman_result_all2$Cancer=factor(spearman_result_all2$Cancer,levels = cancer_names)
spearman_result_all2=spearman_result_all2[order(spearman_result_all2$Cancer),]

spearman_result_all2=spearman_result_all2[which(spearman_result_all2$TRIM_2 != spearman_result_all2$TRIM_1),]

spearman_result_all22=spearman_result_all2 %>% dplyr::filter(Cancer %in% cancer_names_analysis)
spearman_result_all22$Cancer=factor(spearman_result_all22$Cancer,levels = cancer_names_analysis)
spearman_result_all22=spearman_result_all22[order(spearman_result_all22$Cancer),]


result_data=spearman_result_all22

compare_pair=list(c('Different','Same'))

p=ggboxplot(result_data,x='Cancer',y='R',
            fill = 'Family_Type',add ='none',outlier.shape=NA,
            bxp.errorbar=FALSE,width=0.5)+
  scale_fill_manual(values = c('#4393C3','#D6604D'))+#c('#56B4E9','#E69F00'))+
  stat_compare_means(aes(group = Family_Type),label = 'p.signif')+
  theme(axis.text.x = element_text(angle=45,hjust = 1,colour = 'black',size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16))+
  geom_vline(xintercept=c(seq(1.5,32.5,by=1)), linetype="dotted")


pdf(paste0('Spearman_boxplot_diff_same.pdf'),
       width = 15,height = 6)
print(p)
dev.off()

#---------- circos

library('circlize')
spearman_result_all2_sig=read.csv(paste0('Spearman_TRIM_sig.csv'),header=T,stringsAsFactors=F)
spearman_result_all2_sig=as.data.frame(spearman_result_all2_sig)

spearman_result_all2_sig[,c(3,4,5)]=apply(spearman_result_all2_sig[,c(3,4,5)],2,as.numeric)
spearman_result_all2_sig=spearman_result_all2_sig[order(spearman_result_all2_sig$Cancer),]
spearman_result_all2_sig=spearman_result_all2_sig[which(spearman_result_all2_sig$TRIM_2 != spearman_result_all2_sig$TRIM_1),]
spearman_result_all2_sig2=dplyr::filter(spearman_result_all2_sig,
                                        R > 0.75)


spearman_result_all2_sig2=spearman_result_all2_sig2 %>% distinct()

df=spearman_result_all2_sig2 %>% dplyr::select(Symbol_1,Family_1,Symbol_2,Family_2,R,Cancer)

TRIM_2=apply(df,1,function(a){
  a=as.character(a)
  b=sort(c(a[1],a[3]))
  b2=paste(b,collapse='_')
  return(b2)
})
df$type=TRIM_2
df[duplicated(df$type),]
df[which(df$type %in% 'TRIM73_TRIM74'),]
df=dplyr::distinct(df,Cancer,type,.keep_all = TRUE)

df$Symbol_1_1=df$Symbol_1
#df[which(df$R<0),'Symbol_1_1']=paste0(df[which(df$R<0),'Symbol_1_1'],'-')
df$Symbol_2_1=df$Symbol_2
#df[which(df$R<0),'Symbol_2_1']=paste0(df[which(df$R<0),'Symbol_2_1'],'-')

pdf(paste0('Spearman_circos0.75.pdf'),width = 15,height = 15)

brand = c(structure(df$Family_1, names=df$Symbol_1_1), 
          structure(df$Family_2,names= df$Symbol_2_1))

brand = brand[!duplicated(names(brand))]
brand = brand[order(brand, names(brand))]
brand=factor(brand,levels = c('C-I','C-II','C-IV','C-V','C-VI','C-VII','C-VIII','C-IX','C-XI','UC'))

brandord=names(brand)
brandord=gsub('\\+','',brandord)
brandord=gsub('-','',brandord)
kk=match(brandord,Trimdata$Symbol)
order(kk)
brand=brand[order(kk)]

mycolor<-brewer.pal(12,"Paired")[c(1,2,4,5,6,7,8,9,11,12)]
brand_color = structure(mycolor, names = levels(brand))
model_color = structure(rep(mycolor,as.numeric(table(brand))), names = names(brand))
model_color[grep('-',names(model_color))]='grey'

circos.clear()

kk=table(brand)
kk=kk[which(kk!=0)]
library(circlize)
gap.after = do.call("c", lapply(kk, function(i) c(rep(2, i-1), 8)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

chordDiagram(df[, c(8,9)], order = names(brand), grid.col = model_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02)
               ),
             link.lty=0.5,
             annotationTrackHeight = mm_h(c(14, 2))
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


###############################################################################################








