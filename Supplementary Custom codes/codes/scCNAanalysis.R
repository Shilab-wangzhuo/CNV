#===================================================
rm(list = ls())
args=commandArgs(TRUE)

library(GenomicRanges)
library(HMMcopy)
library(DNAcopy)


#set working direct
setwd("/pathway/")
getwd()

mfile='MapabilityAlign40mer_win500kb_hg19.wig'
gfile='gc_win500kb_hg19.wig'
wigData=list.files(pattern='.wig$')

sample.name=data.frame(matrix(unlist(strsplit(wigData,".wig")), nrow=length(wigData), byrow=T),stringsAsFactors=FALSE)[,1]

alpha=0.02;
nperm=1000;
undo.SD=1;
min.width=5;

for(i in 1:length(wigData)){
  #read data
  print(i)
  rfile=wigData[i]
  readcounts <- wigsToRangedData(rfile, gfile, mfile)
  
  # Correct reads into copy number
  corrected_copy <- correctReadcount(readcounts);
  #  copy is cor.map values after log base 2
  write.table(corrected_copy,paste(sample.name[i],'corrected_copy_500k.txt',sep='_'),sep='\t',quote=FALSE, row.names=FALSE)
  # Segment copy number profile
  #segmented_copy <- HMMsegment(corrected_copy)
  
  corrected_copy_dnacopy=as.data.frame(corrected_copy);
  corrected_copy_dnacopy$chr=as.character(corrected_copy_dnacopy$chr)
  corrected_copy_dnacopy[which(corrected_copy_dnacopy[,1]=='chrX'),1]='chr23';
  corrected_copy_dnacopy[which(corrected_copy_dnacopy[,1]=='chrY'),1]='chr24';
  corrected_copy_dnacopy$chr=as.integer(substr(corrected_copy_dnacopy$chr,4,nchar(corrected_copy_dnacopy$chr)))
  corrected_copy_dnacopy <- corrected_copy_dnacopy[sort(corrected_copy_dnacopy$chr,index.return=TRUE)$ix,]
  corrected_copy_dnacopy=na.omit(corrected_copy_dnacopy);
  
  library(DNAcopy)
  CNA.object <- CNA(corrected_copy_dnacopy$copy, corrected_copy_dnacopy$chr, corrected_copy_dnacopy[,2], data.type="logratio", sampleid=sample.name[i]) 
  smoothed.CNA.object <- smooth.CNA(CNA.object) 
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width,verbose=2) 
  #alpha=0.01, min.width=5, undo.splits="prune", undo.prune=0.05
  ##get short segment lists
  thisShort <- segment.smoothed.CNA.object[[2]]
  ## get segments long data
  seg.mean <- rep(2^thisShort$seg.mean, thisShort$num.mark)
  corrected_copy_dnacopy$seg.mean <- seg.mean;
  write.table(thisShort, file=paste(sample.name[i],'corrected_DNAcopySegment_500k.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE) 
  write.table(corrected_copy_dnacopy, file=paste(sample.name[i],'corrected_DNAcopySegWithRatio_500k.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE);
  
}

filePattern='corrected_DNAcopySegWithRatio_500k.txt$'
copyfile=list.files(pattern=filePattern);
varbinfile=list.files(pattern=filePattern);

copyfileSUbnames=data.frame(matrix(unlist(strsplit(copyfile,"_")), nrow=length(copyfile), byrow=T),stringsAsFactors=FALSE)[,1]

MAPD_value=matrix(NA,length(varbinfile),2);
row.names(MAPD_value)=data.frame(matrix(unlist(strsplit(varbinfile,"_")), nrow=length(varbinfile), byrow=T),stringsAsFactors=FALSE)[,1]

for(i in 1:length(varbinfile)){
  tempcopyfile=read.table(varbinfile[i],header=T,sep='\t');
  tempBincount=tempcopyfile$reads+1; 
  tempBincount_log2Ratio=log2(tempBincount / median(tempBincount));
  log2Ratio1=tempBincount_log2Ratio[c(1:(length(tempBincount_log2Ratio)-1))]
  log2Ratio2=tempBincount_log2Ratio[c(2:length(tempBincount_log2Ratio))]
  MAPD_value[i,1]=median(abs(log2Ratio2-log2Ratio1))
}
write.table(MAPD_value,'a_MAPD_value.txt',sep='\t')

tiff("@MAPD_samples.tiff")
barplot(MAPD_value[,1],las=2,main='MAPD',cex.names=1,ylab='MAPD')
abline(h = 0.45, col = "red", lty = 1)
dev.off()

nsample=length(sample.name)
ggplot_cnv <- function(copyfile,copyfileSUbnames){
  
  samples=NULL;
  library("ggplot2")
  i=1
  for(i in 1:length(copyfile)){
    tempcopyfile=read.table(copyfile[i],header=T,sep='\t');
    tempcopyfile[,1]=as.numeric(tempcopyfile[,1])
    tempcopyfilenew=cbind(tempcopyfile[,c(1,2)],2^tempcopyfile$copy
                          *2,round(tempcopyfile$seg.mean*2),rep(copyfileSUbnames[i],dim(tempcopyfile)[1]));
    indexXY=which(tempcopyfilenew[,1]==23 | tempcopyfilenew[,1]==24);
    tempcopyfilenew[indexXY,3]=tempcopyfilenew[indexXY,3]/2
    tempcopyfilenew[indexXY,4]=round(tempcopyfile[indexXY,12])	
    colnames(tempcopyfilenew)=c('chrom','chrompos','cn.ratio','copy.number','Sample')
    write.table(tempcopyfilenew,paste(copyfileSUbnames[i],'HmmRatio_DNAcopySegWithRatio_new_copy.txt',sep='_'),row.names=F,col.names=T,sep='\t')
    samples=rbind(samples,c(i,paste(copyfileSUbnames[i],'HmmRatio_DNAcopySegWithRatio_new_copy.txt',sep='_')));
    
  }
  write.table(samples,'AllSamples.txt',row.names=F,col.names=F,sep='\t')
  
  samples_subForppt=data.frame(matrix(unlist(strsplit(samples[,2],"_")), nrow=dim(samples)[1], byrow=T),stringsAsFactors=FALSE)
  write.table(samples_subForppt[,1],'AllSamples_list.txt',row.names=F,col.names=F,sep='\t')
  
  samp <- read.table('AllSamples.txt', header=FALSE)
  num <- length(samp$V1)
  
  for (j in c(1:num)){
    file=as.character(samp$V2[j])
    if(length(grep(".gz", file))>0){dataTable<-read.table(gzfile(file), header=TRUE)
    }else{dataTable<-read.table(file, header=TRUE)}
    
    dataTable<-cbind(dataTable, Type=NA)
    
    tt <- which(dataTable$copy.number == 2 & dataTable$chrom <= 22)
    dataTable$Type[tt]="Normal"
    tt <- which(dataTable$copy.number > 2 & dataTable$chrom <= 22)
    dataTable$Type[tt]="Amplification"
    tt <- which(dataTable$copy.number < 2 & dataTable$chrom <= 22)
    dataTable$Type[tt]="Deletion"
    
    tt <- which(dataTable$copy.number == 1 & dataTable$chrom > 22)
    dataTable$Type[tt]="Normal"
    tt <- which(dataTable$copy.number > 1 & dataTable$chrom > 22)
    dataTable$Type[tt]="Amplification"
    tt <- which(dataTable$copy.number < 1 & dataTable$chrom > 22)
    dataTable$Type[tt]="Deletion"
    
    if(j==1){Merge<-rbind(dataTable)
    }else{Merge<-rbind(Merge, dataTable)}
  }
  
  cols <- c("Normal" = "#5DF300","Amplification" = "#FF0033","Deletion" = "#0020E4")
  #cols <- c("Normal" = "#5DF300","Amplification" = "#FF0033","Deletion" = "#0020E4")
  Merge$chrom[which(Merge$chrom == 23)] <- "X"
  Merge$chrom[which(Merge$chrom == 24)] <- "Y"
  Merge$chrom <- factor(Merge$chrom, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
  Merge$Sample <- factor(Merge$Sample,levels=as.vector(copyfileSUbnames))
  
  theme_set(theme_bw())
  p1=ggplot(Merge, aes(group=Sample)) +
    geom_point(aes(x=chrompos,y=cn.ratio, colour=factor(Type)), size=1, alpha=2/3) +
    geom_point(aes(x=chrompos,y=copy.number), size=1, alpha=2/3) +
    theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), legend.position="bottom", 
          legend.background=element_blank(), legend.key=element_blank(), legend.title=element_blank(), legend.text=element_text(size=12), 
          strip.text=element_text(size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14), 
          panel.grid=element_blank()) +
    scale_colour_manual(values = cols,breaks= c("Amplification", "Normal", "Deletion")) +
    labs(x="Chromosome", y="Copy Number", size=12) +
    facet_grid(Sample ~ chrom, space="free_x", scales="free_x") +
    guides(colour = guide_legend(override.aes=list(size=3))) +
    ylim(0,8)
  
  ggsave('@tumorCNA_samples.tiff', p1, width = 15,height = nsample+1, dpi=300,limitsize = FALSE)
}
ggplot_cnv(copyfile,copyfileSUbnames)

