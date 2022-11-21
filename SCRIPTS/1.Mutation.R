remove(list = ls())
source('Code/Defined_Functions.R')

############# Calculate mutation frequency of TRIM subfamily ###############
setwd('mutation/mut_fre')
mut_all=readtxt_2(workdir3,'all_cancer_mutdata.txt')

clinicdall<-c()
All_frequence=data.frame()

for(i in 1:length(cancer_names)){

cancer<-as.character(cancer_names[i])
data_TCGA=mut_all%>%filter(Cancer == cancer)

data_TCGA2<-data_TCGA %>% filter(Variant_Type=='SNP') %>% 
  dplyr::select(Hugo_Symbol,Variant_Type,Variant_Classification,
         Tumor_Sample_Barcode,
         HGVSc,HGVSp_Short)
data_TCGA2<-as.data.frame(data_TCGA2)
sample_number<-unique(data_TCGA2$Tumor_Sample_Barcode) %>% length()

data_TCGA4 = data_TCGA2 %>% filter(Hugo_Symbol %in% Trimdata_Symbol) %>% 
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode) %>% distinct()

data_TCGA4=merge(data_TCGA4,Trimdata[,c('Symbol','Family')],by.x='Hugo_Symbol',by.y='Symbol')

data_frequence=data_TCGA4 %>% distinct(Tumor_Sample_Barcode,Family) %>% group_by(Family) %>% dplyr::summarise(n=n()) 
family_num=data_TCGA4 %>% dplyr::select(Hugo_Symbol,Family) %>% distinct() %>% group_by(Family) %>% dplyr::summarise(n=n())
data_frequence=as.data.frame(data_frequence)
colnames(data_frequence)[2]='Family_mut_num'
data_frequence=merge(data_frequence,family_num,by='Family')
data_frequence$sample_number=sample_number
data_frequence$frequence_gene=data_frequence$Family_mut_num / data_frequence$sample_number / data_frequence$n
data_frequence$Cancer=cancer

All_frequence<-rbind(All_frequence,data_frequence)
}
write.csv(All_frequence,paste0('TRIM_subfamily_mut_fre.csv'),quote=F,row.names = F)

Family_num<-as.data.frame(table(Trimdata$Family))
colnames(Family_num)=c('Family','all_gene')
All_frequence=read.csv(paste0('TRIM_subfamily_mut_fre.csv'))
colnames(All_frequence)[3]='mut_gene'

All_frequence=merge(All_frequence,Family_num,by='Family')
All_frequence$New_Mut_Freq=All_frequence$mut_gene * All_frequence$Family_mut_num / All_frequence$all_gene / All_frequence$sample_number

write.table(All_frequence,'New_Mut_Freq_all.txt',sep='\t',quote=F,row.names = F)

#boxplot
result_num_all=read.csv(paste0('New_Mut_Freq_all.txt'),sep = '\t')
library(RColorBrewer)
mycols<-brewer.pal(12,"Paired")
names(mycols)=c('C-I','C-II','C-III',
  'C-IV','C-V','C-VI',
  'C-VII','C-VIII','C-IX',
  'C-X','C-XI','UC')
result_num_all$Family=factor(result_num_all$Family,levels = c('C-I','C-II','C-III',
                                                              'C-IV','C-V','C-VI',
                                                              'C-VII','C-VIII','C-IX',
                                                              'C-X','C-XI','UC'))
result_num_all=result_num_all[order(result_num_all$Family),]

p=ggplot(data=result_num_all, aes(x=Family,y=New_Mut_Freq,color=Family))+
  stat_boxplot(geom='errorbar',width=0.15)+
  geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
  theme_bw()+
  labs(title=paste('Mutation of Family in Cancers'))+
  ylab('Mutation Frequency')+xlab('TRIM Families')+
  geom_jitter(aes(color=Family),shape=16, position=position_jitter(0.2))+
  theme(
    plot.title = element_text(hjust = 0.5, face =  "bold"),
    legend.position = "none")+
  scale_color_manual(values=mycols)

mycompare=list(c('C-IV','C-I'),c('C-IV','C-II'),c('C-IV','C-III'),
               c('C-IV','C-V'),c('C-IV','C-VI'),c('C-IV','C-VII'),
               c('C-IV','C-VIII'),c('C-IV','C-IX'),c('C-IV','C-X'),
               c('C-IV','C-XI'),c('C-IV','UC')
)

pdf('mutfrequency_family_boxplot_P.pdf',width=10,height=6)
p+stat_compare_means(comparisons = mycompare,label = 'p.signif',
                     # label.y = c(0.21,0.18,0.15,seq(0.15,0.36,0.03)))
                     label.y = c(0.29,0.26,0.23,seq(0.23,0.6,0.03)[1:9]))+
  ylim(NA,0.47)
dev.off()

#################################################################



############# Annovar #############
mutdata_all2=readtxt_2(workdir3,'all_cancer_mutdata.txt')

for(i in 1:length(cancer_names)){
   cancer_name0=cancer_names[i]
   testdata<-mutdata_all2[which(mutdata_all2$Cancer == cancer_name0),]
   testdata2<-dplyr::select(testdata,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,
                            Hugo_Symbol,Transcript_ID,Gene,Tumor_Sample_Barcode,Cancer,
                            Variant_Classification,Variant_Type,VARIANT_CLASS,
                            HGVSc,HGVSp_Short,Protein_position,TREMBL)
   write.table(testdata2,paste0("TCGA_mut_for_Annovar2/",cancer_name0,".avinput"),sep='\t',quote = F,col.names = F,row.names = F)
}


```Linux
annovar_path='annovar/'
cancer_annovar_path='mutation/mut_effect/Annovar/'
cancer_names=("ACC" "BLCA" "BRCA" "CESC" "CHOL" "COAD" "DLBC" "ESCA" "GBM" "HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG" "LIHC" "LUAD" "LUSC" "MESO" "OV" "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS" "UVM")
for i in ${cancer_names[@]};do
{
  echo 'perl '${annovar_path}'table_annovar.pl '${cancer_annovar_path}'TCGA_mut_for_Annovar2/'${i}'.avinput '${annovar_path}'humandb/ -buildver hg19 -out '${cancer_annovar_path}'Result2/'${i}' -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -otherinfo' >> ${cancer_annovar_path}Annovar_run2.sh
}
done
#Running in the Linux background
nohup sh Annovar_run2.sh > Annovar_run2.log 2>&1 &
```

fileall=list.files(path=paste0('Result2'),
                   'hg19_multianno.csv',recursive=T,full.names=T)

result_all=data.frame()
for (i in 1:length(fileall)) {
  fileall0=fileall[i]
  fileall1=read.csv(fileall0,sep = ',',header=T)
  fileall1=fileall1[,-c(56:ncol(fileall1))]
  fileall1=separate(fileall1, Otherinfo1, 
                    c('Hugo_Symbol','Transcript_ID','Gene','Tumor_Sample_Barcode','Cancer',
                      'Variant_Classification','Variant_Type','VARIANT_CLASS',
                      'HGVSc','HGVSp_Short','Protein_position','TREMBL'), 
                    sep = "\t")
  result_all=rbind(result_all,fileall1)
}
result_all=distinct(result_all)
write.table(result_all,'Annovar_allcancer_result.txt',sep='\t',quote = F,row.names = F)

#################################################################



############# Calculate mutation frequency of domain ###############
setwd('../mutation/domain/')

#---------------Process domain information

#Protein name and its length information
TRIM_uniprot<-read.xlsx(paste0(workdir1,'TRIM_uniprot.xlsx'))
TRIM_uniprot_need<-TRIM_uniprot %>% dplyr::select(TRIM,Entry,Length)#[,c(1,2,8)]

#Domain information of TRIM
TRIMdat_need=Trimdata[,c(5,10:ncol(Trimdata))]
TRIMdat_need[is.na(TRIMdat_need)]<-''
write.csv(TRIMdat_need,paste0('TRIM_Domain_range.csv'),quote=F,row.names = F)

rownames(TRIMdat_need)<-TRIMdat_need[,1]
TRIMdat_need<-TRIMdat_need[,-1]
TRIMdat_need2<-reshape2::melt(as.matrix(TRIMdat_need))
colnames(TRIMdat_need2)<-c('TRIM','Domain_name','Range')

#Merge the protein information and domain information corresponding to TRIM
TRIMdat_need3<-merge(TRIMdat_need2,TRIM_uniprot_need,by='TRIM')

#Calculate the length of the domain
domain_len<-function(a){
   a1=a[3]
   if(a1!=''){
      a2=unlist(strsplit(a1,';'))
      a3=strsplit(a2,'-')
      b<-c()
      for(i in a3){
         j=as.numeric(i)
         b=c(b,j[1]:j[2])
      }
      return(length(b))
   }else{
      return(0)
   }
}
myddd<-apply(TRIMdat_need3,1,domain_len)
TRIMdat_need3$Domain_len=myddd
write.table(TRIMdat_need3,paste0('TRIM_Uniprot_Domain_info.txt'),sep='\t',quote = F,row.names=F)

#Determine whether the gene has domains and mutation occur in the domain
reset_domain<-function(a){
  a=as.character(a)
      k2<-a[2]
      if(k2!=''){
         b<-c()
         newa<-strsplit(unlist(strsplit(k2,';')),'-')
         for(m in newa){
            m=as.numeric(m)
            b=c(b,m[1]:m[2])
         }
         # print(b)
         if(as.numeric(a[1]) %in% b){
            dd = 'Yes'
         }else{
            dd = 'No'
         }
      }else{
        dd=''
      }
      return(dd)
}
  

#Read result file from Annovar
mutation=readtxt_2(paste0("../mut_effect/Annovar/"),paste0('Annovar_allcancer_result.txt'))
mutation2=mutation %>% filter(Hugo_Symbol %in% Trimdata$Symbol)
rm(mutation)
write.table(mutation2,'mutation_allcancer_onlytrim.txt',sep='\t',quote = F,row.names = F)

#Select SNP for subsequent analysis
mutation_cancer<-filter(mutation2,Variant_Type=='SNP')
mutation_cancer$Protein_position=as.numeric(mutation_cancer$Protein_position)
rm(mutation2)

TRIM_pfam3_4 <- readtxt_1(filename=paste0('TRIM_Uniprot_Domain_info.txt'))
TRIM_pfam4<- read.csv(paste0('TRIM_Domain_range.csv'),header=T)
colnames(TRIM_pfam4)[1]<-'Gene'
domainall<-colnames(TRIM_pfam4)[2:ncol(TRIM_pfam4)]

mutation_cancer2=mutation_cancer[,c('Hugo_Symbol','Protein_position')] %>% distinct()

for(j in 1:length(domainall)){
   domainall0<-domainall[j]
   TRIM_pfam0<-TRIM_pfam4[,c('Gene',domainall0)]
   TRIM_pfam0<-TRIM_pfam0[which(TRIM_pfam0$Gene %in% unique(mutation_cancer$Hugo_Symbol)),]
   TRIM_pfam0<-TRIM_pfam0[which(TRIM_pfam0[,domainall0]!=''),]
   if(nrow(TRIM_pfam0)>0){
      print(domainall0)
      mutation_cancer2<-merge(mutation_cancer2,TRIM_pfam0,by.x='Hugo_Symbol',by.y='Gene',all=TRUE)
      mutation_cancer2[is.na(mutation_cancer2[,domainall0]),domainall0]<-''
      a<-mutation_cancer2[,c('Protein_position',domainall0)]
      b<-apply(a,1,reset_domain) %>% unlist()
      mutation_cancer2[,c(domainall0)]<-b
   }
}

mutation_cancer2<-distinct(mutation_cancer2)
mutation_cancer3=merge(mutation_cancer,mutation_cancer2,by=c('Hugo_Symbol','Protein_position'))
write.table(mutation_cancer3,paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',quote=F,row.names=F)


#---------------Calculate the mutation frequency in domain
mutation_cancer<-read.csv(paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',header=T)
mutation_cancer=mutation_cancer[,-ncol(mutation_cancer)]
mutation=readtxt_2(paste0(workdir3,"mutation/mut_effect/Annovar/"),paste0('Annovar_allcancer_result.txt'))
mutation=filter(mutation,Variant_Type=='SNP')

mutdom_info_all<-c()
for(i in 1:length(cancer_names)){
  cancer0=cancer_names[i]
  mutation_cancer_0<-filter(mutation_cancer,Cancer==cancer0)
  mutation_0<-filter(mutation,Cancer==cancer0)
  aa2=unique(mutation_0$Tumor_Sample_Barcode)
  mutation_N<-length(aa2)
  mutation_cancer_0<-mutation_cancer_0[,c(1,62,59,2,6,7,67:ncol(mutation_cancer_0))]
  domainall<-colnames(mutation_cancer_0)[7:(ncol(mutation_cancer_0))]
  for(j in 1:length(domainall)){
    domainall0<-domainall[j]
    mutation_cancer_0_domain<-mutation_cancer_0[which(mutation_cancer_0[,domainall0] == 'Yes'),]
    aa<-unique(mutation_cancer_0_domain$Tumor_Sample_Barcode)
    mutation_n<-length(aa)
    mutation_info<-c(cancer0,domainall0,mutation_n,mutation_N)
    mutdom_info_all<-rbind(mutdom_info_all,mutation_info)
  }
}

colnames(mutdom_info_all)<-c('Cancer','Domain_name','mutation_n','mutation_N')
mutdom_info_all<-as.data.frame(mutdom_info_all)
mutdom_info_all$mutation_n<-unlist(lapply(mutdom_info_all$mutation_n,as.numeric))
mutdom_info_all$mutation_N<-unlist(lapply(mutdom_info_all$mutation_N,as.numeric))
mutdom_info_all$mutation_freq<-mutdom_info_all$mutation_n / mutdom_info_all$mutation_N
write.csv(mutdom_info_all,paste0('mutdom_info_fre_all_0324修改.csv'),quote=F,row.names=F)

pdf(paste0('alltrim_cancer_domain.pdf'))
myorder<-c('RING','PYRIN',"B.box.type1",'zf.B_box','coiled_coil','COS',
           'PHD','Filamin','MATH','Arf','TM','fn3','Bromodomain','NHL','PRY','SPRY'
)
mutdom_info_all$Domain_level<-factor(mutdom_info_all$Domain_name,levels = myorder)
mutdom_info_all2=mutdom_info_all[which(mutdom_info_all$Domain_name != 'n.a'),]
p=ggboxplot(mutdom_info_all2,x='Domain_level',y='mutation_freq',color='Domain_name',outlier.shape=NA)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")+
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(19)) +ylim(0,0.15)
p+stat_compare_means(method="anova",label.y = 0.15,label.x=2)
dev.off()

#################################################################


############# Comparison of mutation scores whether the mutation occurs in the structural domain #############
setwd('mutation/domain/')

mutation_cancer<-read.csv(paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',header=T)
mutation_cancer=mutation_cancer[,-ncol(mutation_cancer)]
domain_status=apply(mutation_cancer[,67:ncol(mutation_cancer)], 1, function(a){
  a=as.character(a)
  len=length(which(a=="Yes"))
  if(len>0){
    b='Domain'
  }else{
    b='Others'
  }
  return(b)
})

mutation_cancer$domain_status=domain_status

mutation_cancer0=mutation_cancer
scoreorgin<-colnames(mutation_cancer0)[23:50]
scoreorgin<-grep('pred',scoreorgin,invert=T,value=T)
scoreorgin<-scoreorgin[-which(scoreorgin %in% c('CADD_phred','integrated_confidence_value'))]

scoreorgin_select=scoreorgin#[c(1,3,10,15)]

mutation_cancer_1=mutation_cancer0[,c('Hugo_Symbol','domain_status',scoreorgin_select)]
mutation_cancer_1[,3:ncol(mutation_cancer_1)]=apply(mutation_cancer_1[,3:ncol(mutation_cancer_1)],2,function(a){
  a[a=='.']=NA
  a=as.numeric(a)
  return(a)
})

mutation_cancer_2=reshape2::melt(mutation_cancer_1,id.vars=c('Hugo_Symbol','domain_status'))

mutation_cancer_2$variable=lapply(mutation_cancer_2$variable,function(a){
  a=as.character(a)
  a=unlist(strsplit(a,'_'))[1]
  return(a)
}) %>% unlist()

my_comparisons <- list( c("Domain", "Others") )
p<-ggplot(data=mutation_cancer_2, aes(x=variable,y=value,color=domain_status))+
  stat_boxplot(geom='errorbar',width=0.15,position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA,position = position_dodge(0.9))+
  theme_bw()+
  scale_color_brewer(palette='Set1')+
  labs(title=paste0(''))+
  ylab('Score')+
  xlab('')+
  stat_compare_means(aes(group=domain_status),method = 'wilcox.test',label = 'p.signif')+
  facet_wrap(.~variable,scales = 'free',nrow = 1,labeller ='label_value',strip.position='bottom')+
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.placement = "outside",
    strip.text = element_text(size = 5)
  )
  

pdf(paste0('IsOrNotInDomain/allcancer_domain.pdf'),width = 15,height = 5)
print(p)
dev.off()

#################################################################


############# Calculate the cumulative distribution probability of mutation occurring in the domain #############
#---------------Calculate
mutation_cancer<-read.csv(paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',header=T)
mutation_cancer=mutation_cancer%>% filter(Variant_Classification == 'Missense_Mutation')
mutation_cancer=mutation_cancer[,-ncol(mutation_cancer)]

mydata_all<-data.frame(Hugo_Symbol=character(),Domain_name=character(),mutation_n=numeric(),
                       mutation_N=numeric(),mutation_freq=numeric(),domain_mut=character())

mydomain_Ep<-function(IDRData){
  nROI=as.numeric(IDRData[3])#Number of mutations within the domain
  ntot=as.numeric(IDRData[4])#The number of mutations in the gene
  LROI=as.numeric(IDRData[7])#Length of domain
  LGene=as.numeric(IDRData[8])#Length of the whole protein
  p=1-pbinom(nROI,
             ntot,
             LROI/LGene)
  E=nROI/(ntot*LROI/LGene)
  dat=c(E,p)
  return(dat)
}

a=mutation_cancer
Trimall<-unique(a$Hugo_Symbol)

for(j in 1:length(Trimall)){
  Trimall0=Trimall[j]
  aa<-filter(a,Hugo_Symbol == Trimall0)
  Domainall<-colnames(aa)[67:ncol(aa)]
  aa2=unique(aa$Tumor_Sample_Barcode)
  mutation_N<-length(aa2)
  # k=1
  for(k in 1:length(Domainall)){
    Domainall0=Domainall[k]
    aa_domain<-aa[which(aa[,Domainall0] =='Yes'),]
    aa_domain_mut=paste(aa_domain[,'HGVSp_Short'],collapse = ';')
    aaa<-unique(aa_domain$Tumor_Sample_Barcode)
    mutation_n<-length(aaa)
    mutation_freq<-mutation_n/mutation_N
    mydata<-data.frame(Trimall0,Domainall0,mutation_n,mutation_N,mutation_freq,aa_domain_mut)
    mydata_all<-rbind(mydata_all,mydata)
    
  }}
colnames(mydata_all)<-c('Gene','Domain_name','mutation_n','mutation_N','mutation_freq','domain_mut')
mydata_all<-as.data.frame(mydata_all)

write.csv(mydata_all,paste0('mutdom_trim.csv'),quote=F,row.names=F)
#mydata_all<-read.csv(paste0('mutdom_trim.csv'),header=T)

TRIM_pfam3_4<-read.table(paste0('TRIM_Uniprot_Domain_info.txt'),sep='\t',header=T)
TRIM_pfam3_4_Domainlen<-dplyr::select(TRIM_pfam3_4,Domain_name,TRIM,Domain_len)
TRIM_pfam3_4_Proteinlen<-distinct(dplyr::select(TRIM_pfam3_4,TRIM,Length))

mutdom_domain_info_all<-merge(mydata_all,TRIM_pfam3_4_Domainlen,by.x=c('Gene','Domain_name'),by.y=c('TRIM','Domain_name'),all=TRUE)
mutdom_domain_info_all<-merge(mutdom_domain_info_all,TRIM_pfam3_4_Proteinlen,by=c('Gene'),by.y=c('TRIM'),all=TRUE)
write.csv(mutdom_domain_info_all,paste0(workdir3,'mutation/domain/','mutdom_info_fre_len_each_trim.csv'),quote=F,,row.names=F)
# mutdom_domain_info_all<-read.csv('mutdom_info_fre_len_each_trim.csv',header=T)
mutdom_domain_info_all<-as.data.frame(mutdom_domain_info_all)


mut_dom_p<-apply(mutdom_domain_info_all,1,mydomain_Ep) %>% t()
colnames(mut_dom_p)=c('Evalue','Pvalue')
mutdom_domain_info_all2<-cbind(mutdom_domain_info_all,mut_dom_p)
mutdom_domain_info_all2<-filter(mutdom_domain_info_all2,!is.na(mutation_n) & !is.na(Domain_len) & Domain_len != 0 )
mutdom_domain_info_all2$Padjust=p.adjust(mutdom_domain_info_all2$Pvalue,method = 'fdr')
mutdom_domain_info_all3<-filter(mutdom_domain_info_all2,mutation_n >= 3 & Padjust <= 0.05 & Evalue >= 2)

write.csv(mutdom_domain_info_all2,paste0('mutdom_info_fre_each_trim_EP.csv'),quote=F,row.names=F)
write.csv(mutdom_domain_info_all3,paste0('mutdom_info_fre_each_trim_EP_sig.csv'),quote=F,row.names=F)

#---------------boxplot_1
mutation_cancer<-read.csv(paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',header=T)
mutation_cancer<-mutation_cancer%>%filter(Variant_Classification == 'Missense_Mutation')
#Extract Anovar's mutation score results
scoreorgin<-colnames(mutation_cancer)[23:50]
scoreorgin<-grep('pred',scoreorgin,invert=T,value=T)
scoreorgin<-scoreorgin[-which(scoreorgin %in% c('CADD_phred','integrated_confidence_value'))]
scoreorgin_select=scoreorgin

#Extract the optimal gene and its domain information
mutdom_domain_info_all3=read.csv(paste0('mutdom_info_fre_each_trim_EP_sig.csv'),header = T)
mutdom_domain_info_all3=mutdom_domain_info_all3 %>%dplyr::select(Gene,Domain_name,domain_mut)
mutdom_domain_info_all3=separate_rows(mutdom_domain_info_all3,domain_mut, sep = ";")
mutdom_domain_info_all3=as.data.frame(mutdom_domain_info_all3)
mutdom_domain_info_all3$IsPrioritization='Yes'

mutation_cancer2=merge(mutation_cancer,mutdom_domain_info_all3,
                       by.x=c('Hugo_Symbol','HGVSp_Short'),
                       by.y=c('Gene','domain_mut'),
                       all.x=T
                       )

mutation_cancer2[is.na(mutation_cancer2$IsPrioritization),'IsPrioritization']='No'


mutation_cancer_1=mutation_cancer2[,c('Hugo_Symbol','IsPrioritization',scoreorgin_select)]
mutation_cancer_1[,3:ncol(mutation_cancer_1)]=apply(mutation_cancer_1[,3:ncol(mutation_cancer_1)],2,function(a){
  a[a=='.']=NA
  a=as.numeric(a)
  return(a)
})


mutation_cancer_2=reshape2::melt(mutation_cancer_1,id.vars=c('Hugo_Symbol','IsPrioritization'))

mutation_cancer_2$variable=lapply(mutation_cancer_2$variable,function(a){
  a=as.character(a)
  a=unlist(strsplit(a,'_'))[1]
  return(a)
}) %>% unlist()


my_comparisons <- list( c("Yes", "No") )
p<-ggplot(data=mutation_cancer_2, aes(x=variable,y=value,color=IsPrioritization))+
  stat_boxplot(geom='errorbar',width=0.15,position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA,position = position_dodge(0.9))+
  theme_bw()+
  scale_color_brewer(palette='Set1',direction = -1)+
  labs(title=paste0(''))+
  ylab('Score')+
  stat_compare_means(aes(group=IsPrioritization),method = 'wilcox.test',label = 'p.signif')+
  # ylim(-10,10)
  facet_wrap(.~variable,scales = 'free',nrow = 1,labeller ='label_value',strip.position='bottom')+
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.placement = "outside",
    strip.text = element_text(size = 5)
    
  )

pdf(paste0('IsOrNotInDomain/allcancer_domain_IsPrioritization_selected.pdf'),width =15,height = 5)
print(p)
dev.off()

#---------------boxplot_2

mutation_cancer<-read.csv(paste0('mutation_allcancer_onlytrim_domain.txt'),sep='\t',header=T)
mutation_cancer<-mutation_cancer%>%filter(Variant_Classification == 'Missense_Mutation')
scoreorgin<-colnames(mutation_cancer)[23:50]
scoreorgin<-grep('pred',scoreorgin,invert=T,value=T)
scoreorgin<-scoreorgin[-which(scoreorgin %in% c('CADD_phred','integrated_confidence_value'))]
scoreorgin_select=scoreorgin

mutdom_domain_info_all3=read.csv(paste0('mutdom_info_fre_each_trim_EP_sig.csv'),header = T)
mutdom_domain_info_all3=mutdom_domain_info_all3 %>%dplyr::select(Gene,Domain_name,domain_mut)
mutdom_domain_info_all3=separate_rows(mutdom_domain_info_all3,domain_mut, sep = ";")
mutdom_domain_info_all3=as.data.frame(mutdom_domain_info_all3)
mutdom_domain_info_all3$IsPrioritization='Yes'


mutation_cancer=mutation_cancer %>% filter(Hugo_Symbol %in% mutdom_domain_info_all3$Gene)
mutation_cancer2=merge(mutation_cancer,mutdom_domain_info_all3,
                       by.x=c('Hugo_Symbol','HGVSp_Short'),
                       by.y=c('Gene','domain_mut'),
                       all.x=T
)

mutation_cancer2[is.na(mutation_cancer2$IsPrioritization),'IsPrioritization']='No'

mutation_cancer_1=mutation_cancer2[,c('Hugo_Symbol','IsPrioritization','Tumor_Sample_Barcode','HGVSp_Short',scoreorgin_select)]
mutation_cancer_1[,5:ncol(mutation_cancer_1)]=apply(mutation_cancer_1[,5:ncol(mutation_cancer_1)],2,function(a){
  a[a=='.']=NA
  a=as.numeric(a)
  return(a)
})


mutation_cancer_2=reshape2::melt(mutation_cancer_1,id.vars=c('Hugo_Symbol','IsPrioritization','Tumor_Sample_Barcode','HGVSp_Short'))
mutation_cancer_2$variable=lapply(mutation_cancer_2$variable,function(a){
  a=as.character(a)
  a=unlist(strsplit(a,'_'))[1]
  return(a)
}) %>% unlist()

mutation_cancer_3=mutation_cancer_2 %>% dplyr::select(Hugo_Symbol,IsPrioritization,variable,value) %>% distinct()

gene=unique(mutation_cancer_3$Hugo_Symbol)

for (z in seq_along(gene)) {
gene0=gene[z]
  
mutation_cancer_4=mutation_cancer_3 %>% filter(Hugo_Symbol == gene0)

kk=as.data.frame(table(mutation_cancer_4[,1:3]))
toolnamesbad=as.character(unique(kk[which(kk$Freq < 5),'variable']))
mutation_cancer_4 = mutation_cancer_4 %>% filter(!variable %in% toolnamesbad)

my_comparisons <- list( c("Yes", "No") )
p<-ggplot(data=mutation_cancer_4, aes(x=variable,y=value,color=IsPrioritization))+
  stat_boxplot(geom='errorbar',width=0.15,position = position_dodge(0.9))+
  geom_boxplot(outlier.shape=NA,position = position_dodge(0.9))+
  theme_bw()+
  scale_color_brewer(palette='Set1',direction = -1)+
  labs(title=paste0(''))+
  ylab('Score')+
  stat_compare_means(aes(group=IsPrioritization,label = paste0('p=',..p.format..)),method = 'wilcox.test')+
  # ylim(-10,10)
  facet_wrap(.~variable,scales = 'free',nrow = 1,labeller ='label_value',strip.position='bottom')+
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.placement = "outside",
    strip.text = element_text(size = 5)
    
  )

pdf(paste0('IsOrNotInDomain/allcancer_domain_IsPrioritization_selectedgene_',gene0,'.pdf'),width =15,height = 5)
print(p)
dev.off()

try(dev.off())
}


#################################################################



############# Mutation related lollipop diagram #############
mutation<-read.csv(paste0('mutation_allcancer_onlytrim.txt'),sep='\t',header=T)
mutation_cancer<-filter(mutation,Variant_Classification == 'Missense_Mutation')
mutation_cancer<-dplyr::select(mutation_cancer,Hugo_Symbol,
                               Tumor_Sample_Barcode ,HGVSp_Short,
                               Protein_position,Ref,Alt)
mutation_cancer=distinct(mutation_cancer)


goodgene = read.csv('mutdom_info_fre_each_trim_EP_sig.csv',header=T)
goodgene = goodgene %>% dplyr::select(Gene,Domain_name,domain_mut,Evalue,Padjust)
goodgene$IsPrioritization='Yes'
goodgene2 = goodgene %>% dplyr::select(Gene,Domain_name,domain_mut,IsPrioritization)
goodgene2 = goodgene2 %>% separate_rows(domain_mut, sep = ";")
goodgene2=as.data.frame(goodgene2)

mutation_cancer=merge(mutation_cancer,goodgene2,
                      by.x=c('Hugo_Symbol','HGVSp_Short'),
                      by.y=c('Gene','domain_mut'),
                      all.x=T) %>% distinct()

mutation_cancer[is.na(mutation_cancer$IsPrioritization),'IsPrioritization']='No'

TRIM_pfam3_4 <- read.table(paste0('TRIM_Uniprot_Domain_info.txt'),
                           sep='\t',header=T)
TRIM_pfam3_4<-na.omit(TRIM_pfam3_4)
TRIM_pfam3_4=TRIM_pfam3_4[which(TRIM_pfam3_4$Range != ''),]
TRIM_pfam3_4=TRIM_pfam3_4[which(TRIM_pfam3_4$Domain_name !="n/a"),]
library(tidyr)
TRIM_pfam3_4=TRIM_pfam3_4 %>% separate_rows(Range, sep = ";", convert = TRUE) %>% separate(Range, c("begin", "end"), "-")
colnames(TRIM_pfam3_4)=c('Hugo_Symbol','Domain','Begin','End','Protein','Protein_len','Domain_len')

TRIM_pfam3_4=merge(TRIM_pfam3_4,goodgene,
                   by.x=c('Hugo_Symbol','Domain'),
                   by.y=c('Gene','Domain_name'),
                   all.x=T) %>% distinct()
TRIM_pfam3_4[is.na(TRIM_pfam3_4$IsPrioritization),'IsPrioritization']='No'
myorder<-c('RING','PYRIN',"B.box-type1",'zf-B_box','coiled_coil','COS',
           'PHD','Filamin','MATH','Arf','TM','fn3','Bromodomain','NHL','PRY','SPRY'
)
TRIM_pfam3_4$Domain=factor(TRIM_pfam3_4$Domain,levels = myorder)
TRIM_pfam3_4=TRIM_pfam3_4[order(TRIM_pfam3_4$Domain),]


# Basic settings
bgBorderCol = "black"
domainBorderCol = "black"
axisTextSize = c(1, 1)
legendTxtSize = 1;labPosAngle = 90;domainLabelSize = 1
titleSize = c(2, 1); pointSize = 1.5
domainAlpha = 1
labPosSize = 0.5
lim.pos = 2:7
lim.lab = 1:6
lim.dat = data.table::data.table(pos = lim.pos, lab = lim.lab)
lim.dat[,posRounded := round(pos)]
lim.dat = lim.dat[!duplicated(posRounded)]
lim.pos = lim.dat[,pos]
lim.lab = lim.dat[,lab]

library(wesanderson)
names(wes_palettes)
library(paletteer)

kk= TRIM_pfam3_4 %>% filter(Hugo_Symbol %in% goodgene$Gene)
domains = unique(kk$Domain)

domain_cols=c(wes_palette(name="GrandBudapest1")[c(1,2,4)],
              wes_palette(name="Rushmore1")[c(1:2)],
              wes_palette(name="Royal1")[2:3],
              wes_palette(name="Royal2")[4],
              paletteer_d("Redmonder::dPBIPuOr", 10, type = "continuous")[9],
              paletteer_d("nord::frost")[1],
              wes_palette(name="Moonrise3")[4:5],
              wes_palette(name="GrandBudapest2")
)


domain_cols=colorRampPalette(brewer.pal(8, "Set3"))(8)
prismatic::color(domain_cols)

domain_cols=domain_cols[1:8]
names(domain_cols) = as.character(domains)
prismatic::color(domain_cols)

domain_cols2=as.character(domain_cols)
prismatic::color(domain_cols2)

yes='#E41A1C'
no='#377EB8'

geneall=unique(goodgene$Gene)
#m=1#seq_along(geneall)
for(m in seq_along(geneall)){
gene0=geneall[m]
myorder<-c('RING','PYRIN',"B.box-type1",'zf-B_box','coiled_coil','COS',
           'PHD','Filamin','MATH','Arf','TM','fn3','Bromodomain','NHL','PRY','SPRY'
)

Domain_info<-dplyr::filter(TRIM_pfam3_4,Hugo_Symbol==gene0)
Domain_info$Colors=str_replace_all(Domain_info$IsPrioritization,
                                   c('Yes'=yes,
                                     'No'=no
                                   ))
Domain_info$IsPrioritization=str_replace_all(Domain_info$IsPrioritization,
                                             c('Yes'='2',
                                               'No'='1'
                                             ))
Domain_info$Domain_len2=as.numeric(Domain_info$End) - as.numeric(Domain_info$Begin) +1
Domain_info$Domain=factor(Domain_info$Domain,levels = myorder)
Domain_info=Domain_info[order(Domain_info$Domain),]


Mutation_info<-dplyr::filter(mutation_cancer,Hugo_Symbol==gene0)
Mutation_info$Color=str_replace_all(Mutation_info$IsPrioritization,
                                    c('Yes'=yes,
                                      'No'=no
                                    ))

library(trackViewer)
features <- GRanges(gene0, IRanges(c(1,as.numeric(Domain_info$Begin)), 
                                   width=c(unique(Domain_info$Protein_len),as.numeric(Domain_info$Domain_len2)),
                                   names=c('',as.character(Domain_info$Domain))))
features$fill <- c("white", Domain_info$Colors)
features$color <-c('black',rep('white',nrow(Domain_info)))

features$height <-c(0.053,rep(0.05,nrow(Domain_info)))
features$cex=rep(0.5,nrow(Domain_info)+1)

Mutation_info2=Mutation_info %>% group_by(Hugo_Symbol,HGVSp_Short,Protein_position,Domain_name,IsPrioritization,Color) %>% summarise(n=n())
Mutation_info2=as.data.frame(Mutation_info2)
SNP <- as.numeric(Mutation_info2$Protein_position)
sample.gr <- GRanges(gene0, IRanges(SNP, width=1, names=Mutation_info2$HGVSp_Short))

sample.gr$color <- Mutation_info2$Color
sample.gr$border <-Mutation_info2$Color

sample.gr$label <- as.character(Mutation_info2$n)
sample.gr$label.col <- "white"

sample.gr$score <- Mutation_info2$n
sample.gr$cex <- Mutation_info2$n/2

legend <- list(labels=c('Prioritization','Other'),fill=c(yes,no)) ## legend fill color
sample.gr$SNPsideID=str_replace_all(Mutation_info2$IsPrioritization,
                                    c(
                                      'Yes'='top',
                                      'No'='bottom'
                                    ))

pdf(paste0(gene0,'_lolliplot.pdf'),width = 10,height = 6)

heightmy=lollipop_track(sample.gr, features,yaxis=F,xaxis = F,legend=legend,jitter='node',label_on_feature=T,cex=0.5)

Evalue=unique(Domain_info$Evalue[!is.na(Domain_info$Evalue)])
Padjust=unique(Domain_info$Padjust[!is.na(Domain_info$Padjust)])
cbioSubTitle=paste0(gene0,' in Pan-cancer')
cbio_EP=paste0('E = ',signif(Evalue,3),'\n','adj.P = ',signif(Padjust,digits=3),'')
grid.text(cbio_EP, x=0.1, y=heightmy-heightmy/10, just="top",
          gp=gpar(cex=1, fontface="bold"))

test_Prioritization=Domain_info[which(Domain_info$IsPrioritization=='2'),'Domain']
mut_number=Mutation_info[which(Mutation_info$IsPrioritization=='Yes'),]
grid.text(paste0('Prioritization', " [",paste(test_Prioritization,collapse = ',') , "]"), 
          x=.1, y=heightmy-heightmy*3/10, just="top",
          gp=gpar(cex=1, fontface="bold"))
mutation_num=nrow(Mutation_info[which(Mutation_info$IsPrioritization=='Yes'),])
patient_num=length(unique(Mutation_info[which(Mutation_info$IsPrioritization=='Yes'),'Tumor_Sample_Barcode']))
grid.text(paste0(mutation_num,' mutations\nin ',patient_num,' samples'), x=.1, y=heightmy-heightmy*4/10, just="top",
          gp=gpar(cex=1, fontface="bold"))

grid.text(paste0('Other'), 
          x=.1, y=heightmy-heightmy*7/10, just="bottom",
          gp=gpar(cex=1, fontface="bold"))
mutation_num=nrow(Mutation_info[which(Mutation_info$IsPrioritization=='No'),])
patient_num=length(unique(Mutation_info[which(Mutation_info$IsPrioritization=='No'),'Tumor_Sample_Barcode']))
grid.text(paste0(mutation_num,' mutations\nin ',patient_num,' samples'), x=.1, y=heightmy-heightmy*9/10, just="top",
          gp=gpar(cex=1, fontface="bold"))

dev.off()
}
#################################################################


