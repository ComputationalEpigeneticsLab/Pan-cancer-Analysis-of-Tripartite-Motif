# work path：
##################### work path #########################################
#work_path='D:/Document/Postgraduate_study/Subject/one/'
work_path = '/boot3/bio_gaoyueying/Study/one2020/'
mutpath=paste0(work_path,'Data/TCGA_mut/')
exppath=paste0(work_path,'Data/cancer_expression/')
workdir1=paste0(work_path,'Data/')
workdir3=paste0(work_path,'Data5/')
Work_mut=paste0(work_path,'Data/TCGA_mut/')
path_mutscore<-paste0(workdir3,'mutation/Annovar/')
setwd(workdir3)
###############################################################



###################### load packages #####################################
# library(tidyr)
library(tidyverse)
library(vroom)
library(openxlsx)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(ggplotify)
library(aplot)
library(ggsignif)
library(ggpubr)
require(ggsci)
library(circlize)
library(ComplexHeatmap)
library("survival")
#install.packages("survminer")
library("survminer")
library(clusterProfiler)#读取基因集文件的read.gmt是该包中的函数
# display.brewer.all()
library(NMF) 
library(corrplot)
library(GSVA)
library(GSEABase)
library(genefilter)
library(Biobase)
library(stringr)
library(utils)
library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos = rforge, dependencies = TRUE)
library(estimate)
library(enrichplot)###3******
library(ggalluvial)
library(networkD3)
library(riverplot)
library(randomcoloR)
library(diffdf)
library(Hmisc)
library(devtools)
#BiocManager::install('simplifyEnrichment')
#install_github("jokergoo/simplifyEnrichment")
library(simplifyEnrichment)
library(GOSemSim)
# BiocManager::install("rrvgo")
library(rrvgo)
library(GEOquery)
library(Seurat)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
#BiocManager::install('ggsci')
library(ggsci)
#install.packages('tidyverse')
library(devtools)
library(patchwork)
library(circlize)
library(wesanderson)
library(ggsci)
# cl=pal_lancet("lanonc",alpha=0.6)(9)
library(scales)
# show_col(cl)
# prismatic::color(cl)
library(trackViewer)

# add Cancers-related colors
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
colorsmy=colors

tumorOrder=read.xlsx(paste0(workdir3,'tumorOrder.xlsx'),sheet=2,colNames =F)
tumorOrder=tumorOrder[,1]
cancer_names=tumorOrder

library(openxlsx)
Trimdata=read.xlsx(paste(workdir1,"TRIM.xlsx",sep = ""))
# head(Trimdata)
Trimdata_Symbol<-Trimdata$Symbol
Family<-unique(Trimdata$Family)
Symbol<-Trimdata$ENSG

Trimdata_Family=Trimdata %>% dplyr::select(Symbol,Family)
Trimdata2=Trimdata[,c('Symbol','ENSG')]

#####################################################################



###################### define functions##########################################################
readtxt_1<-function(filepath='./',filename,sep0='\t',header0=T){
   data<-read.csv(paste(filepath,filename,sep = ""),
                  header=header0,sep=sep0,fill=T,stringsAsFactors=F,check.names = F)
   return(data)
}

readtxt_2<-function(filepath='./',filename,col_names=TRUE){
   data<-vroom::vroom(paste(filepath,filename,sep = ""),col_names=col_names,)
   data=as.data.frame(data)
   return(data)
}

splitsample=function(a){
  a=as.character(a)
  a1=unlist(strsplit(a,'-'))[1:3]
  b=paste(a1,collapse = '-')
  return(b)
}

cutmydata_Ref<-function(a){
  library(stringr)
  str1<-strsplit(a, split = ">")[[1]]
  Ref0<-str1[1]
  Ref0_1<-str_sub(Ref0,nchar(Ref0))
  return(Ref0_1)
}
cutmydata_Alt<-function(a){
  str1<-strsplit(a, split = ">")[[1]]
  Alt0<-str1[2]
  return(Alt0)
}

removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}


transid<-function(){
  ENSG=readtxt(workdir3,'ENSGall.txt',sep0 = '\t',header0 = F)
  ENSG=ENSG$V1
  library(org.Hs.eg.db)
  # keytypes(org.Hs.eg.db)
  # columns(org.Hs.eg.db)
  ensembls <- mapIds(org.Hs.eg.db, keys = ENSG, keytype = "ENSEMBL", column="SYMBOL")
  ensembls=as.matrix(ensembls)
  ensembls=cbind(ensembls,rownames(ensembls))
  colnames(ensembls)=c('Hugo_Symbol','ENSG')
  ensembls=apply(ensembls, 2, as.character)
  ensembls=as.data.frame(ensembls)
  # ensembls2=ensembls[!is.na(ensembls$Hugo_Symbol),]
  write.table(ensembls,'ENSGall_Symbol.txt',sep='\t',quote=F,row.names = F)
  return(ensembls)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
   library(plyr)
   
   # New version of length which can handle NA's: if na.rm==T, don't count them
   length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
   }
   # This does the summary. For each group's data frame, return a vector with
   # N, mean, and sd
   datac <- ddply(data, groupvars, .drop=.drop,
                  .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),          
                       median = median   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm),
                       quantile_25 = quantile(xx[[col]],0.25),
                       quantile_75 = quantile(xx[[col]],0.75)
                     )
                  },
                  measurevar
   )
   
   # Rename the "mean" column    
   #datac <- rename(datac, c("median" = measurevar))
   
   
   datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
   
   # Confidence interval multiplier for standard error
   # Calculate t-statistic for confidence interval: 
   # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
   ciMult <- qt(conf.interval/2 + .5, datac$N-1)
   datac$ci <- datac$se * ciMult
   
   return(datac)
}



#------------------- GSEA plot

#Functions extracted from R package: enrichplot and DOS
gsInfo <- getFromNamespace("gsInfo", "enrichplot")
ggtable <- getFromNamespace("ggtable", "enrichplot")
tableGrob2 <- getFromNamespace("tableGrob2", "enrichplot")
gseaScores <- getFromNamespace("gseaScores", "DOSE")



MyGSEAplot=function(dat,top1,condition){
  #Parameter setting
  ES_geom='line'
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  pvalue_table=F
  
  #dat is a matrix about pathway
  #top1 is the pathway you want to draw
  title=top1#paste(top1,condition,sep='<-')
  dat_top1=dat %>% filter(ID == top1)
  needcancer=unique(dat_top1$Cancer)
  needcancer=as.character(needcancer)
  
  gsdata=data.frame()
  x=data.frame()
  for(i in 1:length(needcancer)){
    needcancer0=needcancer[i]
    geneSetID=top1
	
    x0=readRDS(paste0(needcancer0,'_GSEA_all.Rds'))
    x0@result$ID=gsub('HALLMARK_','',x0@result$ID)
    x0@result$ID=gsub('_',' ',x0@result$ID)
    x0@result$ID=str_to_title(x0@result$ID)
  
    x1=x0@result %>% filter(ID == geneSetID)
    x1$Description=needcancer0
    x1$Color=colors[as.character(needcancer0)] %>% as.character()
    x=rbind(x,x1)
    
    gsdata0 <- gsInfo(x0, geneSetID)
    gsdata0$Description=needcancer0
    gsdata = rbind(gsdata,gsdata0)
  }
  

  library(ggplot2)
  base_size=10
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                        size=1)
  p.res <- p + es_layer +
    theme( legend.title = element_blank(),#legend.position = c(.8, .8),
           legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
  
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(unique(gsdata$Description)) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    library(RColorBrewer)
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
    # pvalue_table=T
  }
  
  # p2 <- p2 +
  # geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  #           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  # theme(legend.position="none") +
  # scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab("Rank in Ordered Dataset") +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  p.res <- p.res + ylab('Enrichment score') 
  p.pos <- p.pos + xlab('Rank of genes') + ylab('SCC') 
  
  
  # color=x$Color
  colorsmy0=colorsmy[intersect(names(colorsmy),unique(gsdata$Description))]
  # if (length(color) == length(geneSetID)) {
  if(length(colorsmy0) == length(unique(gsdata$Description))) {
    p.res <- p.res + scale_color_manual(values=colorsmy0)
    if (length(colorsmy0) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = colorsmy0)
    }
  }
  

  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    pd <- round(pd, 4)
    library(gridExtra)
    library(gtable)
    library(grid)
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  subplots=1:3
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))
  rel_heights=c(1.5, 0.5, 1)
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  aa=aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
  # print(aa)
  return(aa)
}





#------------------- plotLollipops
lollipop_track <- getFromNamespace("lolliplot", "trackViewer")
lollipop <- getFromNamespace("lollipopPlot", "maftools")
lollipop2 <- getFromNamespace("lollipopPlot2", "maftools")
lollipop_track <- getFromNamespace("lolliplot", "trackViewer")
convertHeight2NPCnum <- getFromNamespace("convertHeight2NPCnum", "trackViewer")
setFeatureLayerID <- getFromNamespace("setFeatureLayerID", "trackViewer")
plotFeatureLegend <- getFromNamespace("plotFeatureLegend", "trackViewer")
plot.grid.xaxis <- getFromNamespace("plot.grid.xaxis", "trackViewer")
getHeight <- getFromNamespace("getHeight", "trackViewer")
plotLegend <- getFromNamespace("plotLegend", "trackViewer")
plotFeatures <- getFromNamespace("plotFeatures", "trackViewer")
plotLollipops <- getFromNamespace("plotLollipops", "trackViewer")
handleLegend=getFromNamespace('handleLegend','trackViewer')
handleRanges=getFromNamespace('handleRanges','trackViewer')
cutSNP=getFromNamespace('cutSNP','trackViewer')
lolliplot=getFromNamespace('lolliplot','trackViewer')
jitterLables=getFromNamespace('jitterLables','trackViewer')
plotFeatures=getFromNamespace('plotFeatures','trackViewer')
plotLollipops=getFromNamespace('plotLollipops','trackViewer')
handleRanges=getFromNamespace('handleRanges','trackViewer')
cutSNP=getFromNamespace('cutSNP','trackViewer')
lolliplot=getFromNamespace('lolliplot','trackViewer')
cleanDataMcols=getFromNamespace('cleanDataMcols','trackViewer')
grid.lollipop=getFromNamespace('grid.lollipop','trackViewer')
plotFeatures=getFromNamespace('plotFeatures','trackViewer')
plotLollipops=getFromNamespace('plotLollipops','trackViewer')
reAdjustLabels <- getFromNamespace("reAdjustLabels", "trackViewer")



plotLollipops=function (SNPs, feature.height, bottomHeight, baseline, type, 
                        ranges, yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType, 
                        LINEW, cex, ratio.yx, GAP, pin, dashline.col, side = c("top","bottom"), 
                        jitter = c("node", "label"), main = TRUE) 
{
  side <- match.arg(side)
  jitter <- match.arg(jitter)
  if (side == "top") {
    pushViewport(viewport(y = bottomHeight, height = 1, just = "bottom", 
                          xscale = c(start(ranges), end(ranges)), clip = "off"))
  }
  else {
    pushViewport(viewport(y = bottomHeight + feature.height, 
                          height = 1, just = "top", xscale = c(start(ranges), 
                                                               end(ranges)), yscale = c(1, 0), clip = "off"))
  }
  if (type == "pie.stack" && length(SNPs$stack.factor) > 0) {
    stopifnot(is.vector(SNPs$stack.factor, mode = "character"))
    if (length(SNPs$stack.factor.order) > 0 || length(SNPs$stack.factor.first) > 
        0) {
      warning("stack.factor.order and stack.factor.first are used by this function!", 
              "The values in these column will be removed.")
    }
    stack.factors <- unique(as.character(SNPs$stack.factor))
    stack.factors <- sort(stack.factors)
    stack.factors.order <- seq_along(stack.factors)
    names(stack.factors.order) <- stack.factors
    SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs), 
                       as.character(SNPs$stack.factor))]
    SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
    SNPs$stack.factor.first <- !duplicated(SNPs)
    SNPs.condense <- SNPs
    SNPs.condense$oid <- seq_along(SNPs)
    SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), 
                                  start(SNPs), end(SNPs))
    SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
    SNPs.condense <- lapply(SNPs.condense, function(.ele) {
      .oid <- .ele$oid
      .gr <- .ele[1]
      mcols(.gr) <- NULL
      .gr$oid <- NumericList(.oid)
      .gr
    })
    SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
    SNPs.condense <- sort(SNPs.condense)
    #change
    lab.pos.condense <- jitterLables(start(SNPs.condense), 
                                     xscale = c(start(ranges), end(ranges)), lineW = LINEW * 
                                       cex * 2)
    lab.pos.condense <- reAdjustLabels(lab.pos.condense, 
                                       lineW = LINEW * cex * 2)
    
    condense.ids <- SNPs.condense$oid
    lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
    lab.pos <- lab.pos[order(unlist(condense.ids))]
  }
  else {
    #change
    lab.pos <- jitterLables(start(SNPs), xscale = c(start(ranges), 
                                                    end(ranges)), lineW = LINEW * cex * 2)
    lab.pos <- reAdjustLabels(lab.pos, lineW = LINEW * cex * 2)
  }
  if (length(SNPs) > 0) {
    yaxisat <- NULL
    yaxisLabel <- TRUE
    if (length(yaxis) > 1 && is.numeric(yaxis)) {
      yaxisat <- yaxis
      if (length(names(yaxis)) == length(yaxis)) 
        #change:Add y-axis label to display gene name
        yaxisLabel <- names(yaxis)
      yaxis <- TRUE
    }
    if (yaxis && scoreMax > 1 && !type %in% c("pie", "pie.stack")) {
      if (side == "top") {
        
        grid.yaxis(at = yaxisat, label = yaxisLabel, 
                   main = main, gp = yaxis.gp, vp = viewport(x = 0.5 + 
                                                               ifelse(main, -1, 1) * LINEW, y = feature.height + 
                                                               5.25 * GAP * cex + scoreMax * LINEW * ratio.yx/2 * 
                                                               cex, width = 1, height = scoreMax * LINEW * 
                                                               ratio.yx * cex, yscale = c(0, scoreMax0 + 
                                                                                            0.5)))
      }
      else {
        grid.yaxis(at = yaxisat, label = yaxisLabel, 
                   main = main, gp = yaxis.gp, vp = viewport(x = 0.5 + 
                                                               ifelse(main, -1, 1) * LINEW, y = 1 - (feature.height + 
                                                                                                       5.25 * GAP * cex + scoreMax * LINEW * ratio.yx/2 * 
                                                                                                       cex), width = 1, height = scoreMax * LINEW * 
                                                               ratio.yx * cex, yscale = c(scoreMax0 + 0.5, 
                                                                                          0)))
      }
    }
    if (length(SNPs$alpha) == length(SNPs)) {
      SNPs$alpha[is.na(SNPs$alpha)] <- 0
      if (all(is.numeric(SNPs$alpha))) {
        if (any(SNPs$alpha > 1)) {
          SNPs$alpha <- SNPs$alpha/max(SNPs$alpha)
        }
      }
      else {
        SNPs$alpha <- as.numeric(factor(as.character(SNPs$alpha)))
        SNPs$alpha <- (SNPs$alpha + max(SNPs$alpha))/max(SNPs$alpha)/2
      }
    }
    else {
      SNPs$alpha <- NULL
    }
    if (type == "circle") {
      if (length(SNPs$shape) == length(SNPs)) {
        if (!all(SNPs$shape %in% c("circle", "square", 
                                   "diamond", "triangle_point_up", "triangle_point_down"))) {
          message("shape must be \"circle\", \"square\", \"diamond\", \"triangle_point_up\", or \"triangle_point_down\"")
          SNPs$shape <- as.numeric(factor(SNPs$shape))
          SNPs$shape <- rep(c("circle", "square", "diamond", 
                              "triangle_point_up", "triangle_point_down"), 
                            max(SNPs$shape))[SNPs$shape]
        }
      }
      else {
        SNPs$shape <- NULL
      }
    }
    for (m in seq_along(SNPs)) {
      this.dat <- SNPs[m]
      color <- if (is.list(this.dat$color)) 
        this.dat$color[[1]]
      else this.dat$color
      border <- if (is.list(this.dat$border)) 
        this.dat$border[[1]]
      else this.dat$border
      fill <- if (is.list(this.dat$fill)) 
        this.dat$fill[[1]]
      else this.dat$fill
      alpha <- if (length(this.dat$alpha) > 0) 
        this.dat$alpha[[1]]
      else 1
      lwd <- if (is.list(this.dat$lwd)) 
        this.dat$lwd[[1]]
      else this.dat$lwd
      id <- if (is.character(this.dat$label)) 
        this.dat$label
      else NA
      id.col <- if (length(this.dat$label.col) > 0) 
        this.dat$label.col
      else "black"
      shape <- if (length(this.dat$shape) > 0) 
        this.dat$shape[[1]]
      else "circle"
      rot <- if (length(this.dat$label.rot) > 0) 
        this.dat$label.rot
      else 15
      this.cex <- if (length(this.dat$cex) > 0) 
        this.dat$cex[[1]][1]
      else 1
      this.dashline.col <- if (length(this.dat$dashline.col) > 
                               0) 
        this.dat$dashline.col[[1]][1]
      else dashline.col
      if (length(names(this.dat)) < 1) 
        this.dashline.col <- NA
      this.dat.mcols <- mcols(this.dat)
      this.dat.mcols <- cleanDataMcols(this.dat.mcols, 
                                       type)
      grid.lollipop(x1 = convertX(unit(start(this.dat), 
                                       "native"), "npc", valueOnly = TRUE), 
                    y1 = baseline, 
                    x2 = convertX(unit(ifelse(jitter == "node", lab.pos[m], 
                                              start(this.dat)), "native"), "npc", valueOnly = TRUE), 
                    y2 = feature.height, 
                    y3 = 4 * GAP * cex * 2, 
                    y4 = 2.5 * GAP * cex * 2, 
                    radius = LINEW * cex/2, 
                    col = color, 
                    border = border, 
                    percent = this.dat.mcols, 
                    edges = 100, 
                    type = type, ratio.yx = ratio.yx, pin = pin, 
                    scoreMax = scoreMax * LINEW * cex / 2, 
                    scoreType = scoreType, 
                    id = id, id.col = id.col, 
                    cex = this.cex, lwd = lwd, 
                    dashline.col = this.dashline.col, side = side, 
                    rot = rot, alpha = alpha, shape = shape)
    }
    this.height <- getHeight(SNPs, ratio.yx, LINEW, GAP, 
                             cex, type, scoreMax = scoreMax, level = "data")
    labels.rot <- 60#90
    if (length(names(SNPs)) > 0) {
      if (type == "pie.stack") {
        idx <- !duplicated(names(SNPs))
        lab.pos <- lab.pos[idx]
        SNPs <- SNPs[idx]
      }
      labels.x <- lab.pos
      labels.text <- names(SNPs)
      labels.just <- ifelse(side == "top", "left", "right")
      labels.hjust <- NULL
      labels.vjust <- NULL
      labels.check.overlap <- FALSE
      labels.default.units <- "native"
      labels.gp <- gpar(cex = cex)
      for (label.parameter in c("x", "y", "just", "hjust", 
                                "vjust", "rot", "check.overlap", "default.units", 
                                "gp")) {
        label.para <- paste0("label.parameter.", label.parameter)
        if (label.para %in% colnames(mcols(SNPs))) {
          assign(paste0("labels.", label.parameter), 
                 mcols(SNPs)[, label.para])
        }
      }
      if (!"cex" %in% names(labels.gp)) {
        labels.gp <- c(labels.gp, cex = cex)
      }
      mergeList <- function(.ele) {
        .n <- unique(unlist(lapply(.ele, names)))
        .out <- list()
        if (length(.n) > 0) {
          for (.name in .n) {
            .out[[.name]] <- sapply(.ele, function(.e) {
              if (.name %in% names(.e)) {
                .e[[.name]][1]
              }
              else {
                NA
              }
            })
          }
        }
        else {
          .n <- unique(names(.ele))
          for (.name in .n) {
            .out[[.name]] <- unlist(.ele[names(.ele) %in% 
                                           .name])
          }
        }
        .out
      }
      labels.gp <- mergeList(labels.gp)
      labels.gp[duplicated(names(labels.gp))] <- NULL
      labels.gp <- do.call(gpar, labels.gp)
      if (jitter == "label") {
        rased.height <- 4 * GAP * cex
        guide.height <- 2.5 * GAP * cex
        for (i in seq_along(SNPs)) {
          this.dashline.col <- if (length(SNPs[i]$dashline.col) > 
                                   0) 
            SNPs[i]$dashline.col[[1]][1]
          else dashline.col
          if (length(names(SNPs[i])) < 1) 
            this.dashline.col <- NA
          grid.lines(x = c(start(SNPs[i]), labels.x[i]), 
                     y = c(this.height + feature.height - cex * 
                             LINEW, (this.height + feature.height + rased.height)/1.5), 
                     default.units = labels.default.units, gp = gpar(col = this.dashline.col, 
                                                                     lty = 3))
          grid.lines(x = c(labels.x[i], labels.x[i]), 
                     y = c(this.height + rased.height + feature.height, 
                           (this.height + rased.height + guide.height + 
                             feature.height)/1.5), default.units = labels.default.units, 
                     gp = gpar(col = this.dashline.col, lty = 3))
        }
        this.height <- this.height + rased.height + guide.height
      }
      grid.text(x = labels.x, y = (this.height + feature.height)/1.5, 
                label = labels.text, just = labels.just, hjust = labels.hjust, 
                vjust = labels.vjust, rot = labels.rot, check.overlap = labels.check.overlap, 
                default.units = labels.default.units, gp = labels.gp)
    }
  }
  popViewport()
}



lollipop_track=function (SNP.gr, features = NULL, ranges = NULL, type = "circle", 
          newpage = TRUE, ylab = TRUE, ylab.gp = gpar(cex=1,col = "black",fontface="bold"), 
          yaxis = F, yaxis.gp = gpar(col = "black"), xaxis = TRUE, 
          xaxis.gp = gpar(col = "black"), legend = NULL, cex = 1, dashline.col = "gray80", 
          jitter = c("node", "label"), rescale = FALSE, label_on_feature = T, 
          ...) 
{
  stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList", "list")))
  stopifnot(inherits(features, c("GRanges", "GRangesList", 
                                 "list")))
  jitter <- match.arg(jitter)
  rescale.old <- rescale
  xaxis.old <- xaxis
  if (any(type != "circle" & jitter == "label")) {
    jitter[which(type != "circle" & jitter == "label")] <- "node"
    warning("if jitter set to label, type must be cirle.")
    message("jitter is set to node.")
  }
  SNP.gr.name <- unique(as.character(SNP.gr@seqnames))#deparse(substitute(SNP.gr))
  if (is(SNP.gr, "GRanges")) {
    SNP.gr <- list(SNP.gr)
    if (length(SNP.gr.name) == length(SNP.gr)) {
      names(SNP.gr) <- SNP.gr.name
    }
  }
  len <- length(SNP.gr)
  for (i in seq.int(len)) {
    stopifnot(is(SNP.gr[[i]], "GRanges"))
  }
  if (inherits(features, c("GRangesList", "list"))) {
    for (i in seq_along(features)) {
      stopifnot(`features must be a GRanges or GRangesList object` = is(features[[i]], 
                                                                        "GRanges"))
    }
    features <- features[seq.int(len)]
  }
  else {
    stopifnot(`features must be a GRanges or GRangesList object` = is(features, 
                                                                      "GRanges"))
    features <- list(features)[seq.int(len)]
  }
  TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
  if (any(!type %in% TYPES)) {
    stop("Error in match argument: ", paste0("'type' should be one of '", 
                                             paste(TYPES, collapse = "', '"), "'."))
  }
  types <- rep(type, length = len)[seq.int(len)]
  rm(type)
  legend <- handleLegend(legend, len, SNP.gr)
  ranges <- handleRanges(ranges, SNP.gr, features, len)
  SNP.gr <- cutSNP(SNP.gr, ranges, len)
  height <- 1/sum(lengths(ranges))
  args <- as.list(match.call())
  if (length(args$height0) == 0) {
    height0 <- 0
  }
  else {
    height0 <- args$height0
  }
  if (newpage) 
    grid.newpage()
  for (i in seq.int(len)) {
    if (length(ranges[[i]]) > 1) {
      args$newpage <- FALSE
      for (j in rev(seq_along(ranges[[i]]))) {
        args$ranges <- ranges[[i]][j]
        args$SNP.gr <- SNP.gr[i]
        args$features <- features[[i]]
        args$type <- types[i]
        args$legend <- legend[[i]]
        args$height0 <- height0
        height0 <- do.call(what = lolliplot, args = args)
      }
    }
    else {
      type <- match.arg(types[i], TYPES)
      if (type == "pin") {
        pinpath <- system.file("extdata", "map-pin-red.xml", 
                               package = "trackViewer")
        pin <- readPicture(pinpath)
      }
      else {
        pin <- NULL
      }
      vp <- viewport(x = 0.5, y = height0 + height * 0.5, 
                     width = 1, height = height)
      pushViewport(vp)
      LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
      LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
      totalH <- as.numeric(unit(1, "npc"))
      if (LINEH > totalH/20) {
        LINEH <- totalH/20
      }
      GAP <- 0.2 * LINEH
      ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), 
                                        "npc"))
      SNPs <- SNP.gr[[i]]
      strand(SNPs) <- "*"
      SNPs <- sort(SNPs)
      feature <- features[[i]]
      rescale <- rescale.old
      xaxis <- xaxis.old
      if (is.logical(rescale)[1]) {
        if (rescale[1]) {
          range.tile <- tile(ranges[[i]], n = 5)[[1]]
          if (all(width(range.tile) > 2)) {
            range.tile.cnt <- countOverlaps(range.tile, 
                                            SNPs)
            feature.start <- feature.end <- feature
            end(feature.start) <- start(feature.start)
            start(feature.end) <- end(feature.end)
            range.tile.cnt2 <- countOverlaps(range.tile, 
                                             unique(c(feature.start, feature.end)))
            range.tile.cnt <- range.tile.cnt + range.tile.cnt2
            range.width <- width(ranges[[i]])
            range.tile.width <- log2(range.tile.cnt + 
                                       1)
            range.tile.width <- range.tile.width/sum(range.tile.width)
            range.tile.width <- range.width * range.tile.width
            range.tile.width <- cumsum(range.tile.width)
            range.tile.width <- start(ranges[[i]]) + 
              c(0, round(range.tile.width) - 1)
            rescale <- data.frame(from.start = start(range.tile), 
                                  from.end = end(range.tile), to.start = range.tile.width[-length(range.tile.width)], 
                                  to.end = range.tile.width[-1])
            rescale$to.start[-1] <- rescale$to.start[-1] + 
              1
          }
        }
      }
      else {
        if (is.numeric(rescale)) {
          feature.rd <- disjoin(c(feature, ranges[[i]]))
          feature.segment.points <- sort(unique(c(start(feature.rd), 
                                                  end(feature.rd))))
          feature.segment.points <- feature.segment.points[feature.segment.points >= 
                                                             start(ranges[[i]]) & feature.segment.points <= 
                                                             end(ranges[[i]])]
          rescale <- rescale/sum(rescale, na.rm = TRUE)
          rescale <- rescale[!is.na(rescale)]
          if (length(rescale) == length(feature.segment.points) - 
              1) {
            rescale.ir <- IRanges(feature.segment.points[-length(feature.segment.points)] + 
                                    1, feature.segment.points[-1])
            start(rescale.ir)[1] <- start(rescale.ir)[1] - 
              1
            rescale.ir.width <- sum(width(rescale.ir))
            rescale.ir.new.width <- cumsum(round(rescale.ir.width * 
                                                   rescale, digits = 0))
            rescale <- data.frame(from.start = start(rescale.ir), 
                                  from.end = end(rescale.ir), to.start = feature.segment.points[1] + 
                                    c(0, rescale.ir.new.width[-length(rescale.ir.new.width)]), 
                                  to.end = feature.segment.points[1] + rescale.ir.new.width)
          }
          else {
            stop("The length of rescale is not as same as the number of segments (including features and non-features).")
          }
        }
      }
      if (is.data.frame(rescale)) {
        if (all(c("from.start", "from.end", "to.start", 
                  "to.end") %in% colnames(rescale))) {
          rescale.gr <- function(x) {
            if (is(x, "GRanges")) {
              x.start <- start(x)
              x.end <- end(x)
              y <- c(x.start, x.end)
              x.cut <- cut(y, breaks = c(rescale$from.start[1], 
                                         rescale$from.end + 1), labels = seq.int(nrow(rescale)), 
                           right = FALSE)
              y <- mapply(function(a, b) {
                if (!is.na(b)) {
                  rescale(a, to = c(rescale$to.start[b], 
                                    rescale$to.end[b]), from = c(rescale$from.start[b], 
                                                                 rescale$from.end[b]))
                }
                else {
                  a
                }
              }, y, as.numeric(as.character(x.cut)))
              y <- round(y)
              start(x) <- 1
              end(x) <- y[seq_along(x) + length(x)]
              start(x) <- y[seq_along(x)]
              x
            }
            else {
              x.cut <- cut(x, breaks = c(rescale$from.start[1], 
                                         rescale$from.end + 1), labels = seq.int(nrow(rescale)), 
                           right = FALSE)
              y <- mapply(function(a, b) {
                if (!is.na(b)) {
                  rescale(a, to = c(rescale$to.start[b], 
                                    rescale$to.end[b]), from = c(rescale$from.start[b], 
                                                                 rescale$from.end[b]))
                }
                else {
                  a
                }
              }, x, as.numeric(as.character(x.cut)))
              y <- round(y)
              y
            }
          }
          feature <- rescale.gr(feature)
          SNPs <- rescale.gr(SNPs)
          if (is.logical(xaxis)[1]) {
            xaxis <- c(rescale$to.start[1], rescale$to.end)
            names(xaxis) <- c(rescale$from.start[1], 
                              rescale$from.end)
          }
          else {
            xaxis.names <- names(xaxis)
            if (length(xaxis.names) != length(xaxis)) {
              xaxis.names <- as.character(xaxis)
            }
            xaxis <- rescale.gr(xaxis)
            names(xaxis) <- xaxis.names
          }
        }
      }
      feature$height <- convertHeight2NPCnum(feature$height)
      feature <- setFeatureLayerID(feature, ranges[[i]])
      feature.splited <- split(feature, feature$featureLayerID)
      bottomblank <- plotFeatureLegend(feature, as.numeric(convertY(unit(1, 
                                                                         "line"), "npc")), ranges[[i]], xaxis, xaxis.gp, 
                                       label_on_feature)
      if (length(SNPs$score) > 0) {
        SNPs$score <- sapply(SNPs$score, mean)
      }
      scoreMax0 <- scoreMax <- if (length(SNPs$score) > 
                                   0) 
        ceiling(max(c(SNPs$score, 1), na.rm = TRUE))
      else 1
      if (type == "pie.stack") 
        scoreMax <- length(unique(SNPs$stack.factor))
      if (!type %in% c("pie", "pie.stack")) {
        scoreType <- if (length(SNPs$score) > 0) 
          all(floor(SNPs$score) == SNPs$score)
        else FALSE
        if (length(yaxis) > 1 && is.numeric(yaxis)) {
          if (length(names(yaxis)) != length(yaxis)) {
            names(yaxis) <- yaxis
          }
          scoreMax0 <- max(yaxis, scoreMax0)
          scoreMax <- max(yaxis, scoreMax)
        }
        if (scoreMax > 10) {
          SNPs$score <- 10 * SNPs$score/scoreMax
          scoreMax <- 10 * scoreMax0/scoreMax
          scoreType <- FALSE
        }
        else {
          scoreMax <- scoreMax0*3
        }
      }
      else {
        scoreType <- FALSE
      }
      IsCaterpillar <- length(SNPs$SNPsideID) > 0
      if (IsCaterpillar) {
        if (any(is.na(SNPs$SNPsideID)) || !all(SNPs$SNPsideID %in% 
                                               c("top", "bottom"))) {
          warning("Not all SNPsideID is top or bottom")
          IsCaterpillar <- FALSE
        }
      }
      if (IsCaterpillar) {
        SNPs.top <- SNPs[SNPs$SNPsideID == "top"]
        SNPs.bottom <- SNPs[SNPs$SNPsideID == "bottom"]
      }
      else {
        SNPs.top <- SNPs
        SNPs.bottom <- GRanges()
      }
      if (length(SNPs.bottom) < 1) 
        IsCaterpillar <- FALSE
      if (!IsCaterpillar) {
        bottomblank <- bottomblank
      }
      pushViewport(viewport(x = LINEW + 0.5, y = bottomblank/2 + 
                              0.5, width = 1 - 7 * LINEW, height = 1 - bottomblank, 
                            xscale = c(start(ranges[[i]]), end(ranges[[i]])), 
                            clip = "off"))
      bottomHeight <- 0
      if (IsCaterpillar) {
        bottomHeight <- getHeight(SNPs = SNPs.bottom, 
                                  ratio.yx = ratio.yx, LINEW = LINEW, GAP = GAP, 
                                  cex = cex, type = type, scoreMax = scoreMax, 
                                  level = "data&labels")
        vp <- viewport(y = bottomHeight, just = "bottom", 
                       xscale = c(start(ranges[[i]]), end(ranges[[i]])))
        pushViewport(vp)
        xaxis.gp$col <- "gray"
        plot.grid.xaxis(xaxis, gp = xaxis.gp)
        popViewport()
      }
      else {
        plot.grid.xaxis(xaxis, gp = xaxis.gp)
      }
      baseline <- max(c(feature.splited[[1]]$height/2, 
                        1e-04)) + 0.2 * LINEH
      baselineN <- max(c(feature.splited[[length(feature.splited)]]$height/2, 
                         1e-04)) + 0.2 * LINEH
      feature.height <- plotFeatures(feature.splited, LINEH, 
                                     bottomHeight, label_on_feature)
      #change
      bottomHeight2=bottomHeight-baselineN
      yaxis=F#whether show y axis
      LINEW2=LINEW*1.2
      
      if (length(SNPs.bottom) > 0) {
        plotLollipops(SNPs.bottom, feature.height, bottomHeight2, 
                      baselineN, type, ranges[[i]], yaxis, yaxis.gp, 
                      scoreMax, scoreMax0, scoreType, LINEW2, cex, 
                      ratio.yx, GAP, pin, dashline.col, side = "bottom", 
                      jitter = jitter)
      }
      feature.height <- feature.height + 2 * GAP
      
      #change
      bottomHeight3=bottomHeight+baselineN
      
      if (length(SNPs.top) > 0) {
        plotLollipops(SNPs.top, feature.height, bottomHeight3, 
                      baseline, type, ranges[[i]], yaxis, yaxis.gp, 
                      scoreMax, scoreMax0, scoreType, LINEW2, cex, 
                      ratio.yx, GAP, pin, dashline.col, side = "top", 
                      jitter = jitter)
      }
      this.height <- getHeight(SNPs.top, ratio.yx, LINEW, 
                               GAP, cex, type, scoreMax = scoreMax, level = "data&labels")
      this.height <- this.height + bottomHeight + feature.height
      this.height <- plotLegend(legend[[i]], this.height, 
                                LINEH)
      popViewport()
      this.height <- bottomblank + this.height * (1 - bottomblank)
      if (length(yaxis) > 1 && is.numeric(yaxis)) {
        x <- LINEW
      }
      else {
        x <- unit(3, "lines")
        if (yaxis) {
          x <- LINEW
        }
      }
      vp <- viewport(x = 0.5, y = this.height * 0.5, width = 1, 
                     height = this.height)
      pushViewport(vp)
      if (is.logical(ylab)) {
        if (ylab && length(names(SNP.gr)) > 0) {
          grid.text(names(SNP.gr)[i], x = x, y = bottomHeight+baseline, 
                    rot = 0,#90, 
                    gp = ylab.gp)
        }
      }
      if (is.character(ylab)) {
        if (length(ylab) == 1) 
          ylab <- rep(ylab, len)
        grid.text(ylab[i], x = x, y = bottomHeight+baseline, rot = 0,#90, baselineN
                  gp = ylab.gp)
      }
      popViewport()
      popViewport()
      height0 <- height0 + this.height * height
    }
  }
  return(invisible(height0))
}

####################################################################################################
