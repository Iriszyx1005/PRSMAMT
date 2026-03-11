library(data.table)
library(dplyr)
library(openxlsx)
library(R.utils)
## map to cytoband
cytobandmap = function(chr,pos,hg = "hg38"){
  if(hg == "hg19") cytoband = read.delim("/home/sshen/software/cytoBand.txt",header = F)
  if(hg == "hg38") cytoband = read.delim("/home/sshen/software/cytoBand_hg38.txt",header = F)
  cytoband$CHR = unlist(lapply(cytoband$V1,function(x) strsplit(x,"chr")[[1]][2]))
  
  band = NULL
  for (i in 1:length(chr)){
    cytoband1 = subset(cytoband,CHR == chr[i])
    band[i]  = subset(cytoband1,pos[i]>=V2 & pos[i]<V3)$V4
  }
  
  band = paste0(chr,band)
  return(band)
}
#read GWAS result
gwas_data =fread("/home/sshen/EPYC/DiskA/LC_GWAS/LC_META_shen/01META_All.txt.gz")
gwas_data$Allele1=toupper(gwas_data$Allele1)
gwas_data$Allele2=toupper(gwas_data$Allele2)

data=gwas_data[,c("MarkerName","Allele1","Allele2","Freq1","EffectARE","StdErrARE","PvalueARE","chr","pos")]
data$N=1963325
colnames(data)=c("SNP","A1","A2","freq","b","se","p","chr","pos","N")

#data=data[order(data$chr,-data$pos),]
#data=data[,c("SNP","A1","A2","freq","b","se","p")]

for(i in 1:22){
  x=data[data$chr==i,]
  y=x[,c("SNP","A1","A2","freq","b","se","p","N")]
  write.table(y,file = paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/meta",i,".txt"),row.names = FALSE, quote = FALSE)
  write.table(y,file = paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/meta",i,".ma"),row.names = FALSE, quote = FALSE,sep = "\t")
  print(i)
}

###NG
load("/home/sshen/public/resource/G1000_s_hg38.RData")

G1000_s = subset(G1000_s,!is.na(ID))

ng=read.xlsx("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/41588_2022_1115_MOESM4_ESM.xlsx",sheet = 6,startRow = 1)
colnames(ng)=ng[1,]
ng=ng[-1,]
ng$P_BE1=as.numeric(ng$P_BE1)
ng=ng[ng$P_BE1<5e-8,]

temp=ng
out="/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp"
setwd(out)
temp$id_hg37=paste0("chr",temp$SNP) 
temp$BP=as.numeric(temp$BP)
lift = data.frame(V1 = paste0("chr",temp$CHR),V2 = as.integer(as.numeric(temp$BP-1)),V3 =temp$BP,V4=temp$id_hg37)
write.table(lift,paste0(out,"/input.bed"),row.names = F,col.names = F,quote=F,sep= "\t")

system("/home/sshen/miniconda3/bin/CrossMap bed /home/sshen/software/hg19ToHg38.over.chain.gz input.bed output_ng_GRCh38.txt")

lift1 = fread("output_ng_GRCh38.txt")

lift2=data.frame(pos_hg38=lift1$V3,id_hg37=lift1$V4)

data=merge(temp,lift2,by="id_hg37")
data$ID=paste0(data$CHR,":",data$pos_hg38) 

y=data[,"ID"]
write.table(y,paste0(out,"/ng.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

ng <- fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/ng.txt"),header=F)
colnames(ng) <- c("SNP")
#ng$chr=sapply(ng$SNP,function(x) strsplit(x,":")[[1]][1])
#ng$pos=sapply(ng$SNP,function(x) strsplit(x,":")[[1]][2]) 
#ng$cytoband = cytobandmap(ng$chr,ng$pos)

###+GWAS Catalog
catalog <- read.table("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/gwas-association-downloaded_2025-10-20-MONDO_0008903-withChildTraits.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
catalog=catalog[catalog$'P.VALUE'<5e-8,]             
catalog=catalog[,c("CHR_ID","CHR_POS","SNPS","REGION")]
catalog$SNP=paste0(catalog$CHR_ID,":",catalog$CHR_POS) #685

c1=catalog[catalog$CHR_ID=="",]
c1$SNP=gsub("chr", "",c1$SNPS)
c1 <- c1[!grepl("^rs", c1$SNP), ]

catalog=catalog[catalog$CHR_ID %in% 1:22,]
catalog=rbind(catalog,c1)
#catalog$CHR_ID=as.numeric(catalog$CHR_ID)
#catalog$CHR_POS=as.numeric(catalog$CHR_POS)
#catalog$cytoband = cytobandmap(catalog$CHR_ID,catalog$CHR_POS)

x=data.frame(SNP=catalog[,"SNP"])
x=unique(x) #434
catalog_ng_snp=rbind(ng,x)
catalog_ng_snp=unique(catalog_ng_snp)  #1373
write.table(catalog_ng_snp,paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/catalog_ng_snp.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

catalog_ng_snp$chr=sapply(catalog_ng_snp$SNP,function(x) strsplit(x,":")[[1]][1])         
for(i in 1:22){
  x=catalog_ng_snp[catalog_ng_snp$chr==i,]
  x=x[,"SNP"]
  write.table(x,file = paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/catalog_ng_snp",i,".txt"),row.names = FALSE, col.names = FALSE,quote = FALSE)
  
  y=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/meta",i,".txt"))
  y=y[,"SNP"]
  tmp=rbind(x,y)
  snplist=unique(tmp)
  write.table(snplist,file = paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/snplist1_ukb",i,".txt"),row.names = FALSE, col.names = FALSE,quote = FALSE)
  print(i)
}

###筛选UKB 
cd /home/sshen/Disk_m2/PRS_yxzhang/indepent_snp


for chr in $(seq 1 22); do
nohup /home/sshen/Disk_m2/PRS_yxzhang/GCTA/gcta64 --bfile /home/sshen/DiskA/WGS_10k/chr${chr} \
--chr ${chr} --maf 0.001 \
--cojo-file meta${chr}.ma \
--cojo-cond catalog_ng_snp.txt \
--cojo-slct \
--out new/new_snp${chr} &
  done 

#independent loci
data=NULL
for(i in 1:22){
  file_path <-paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/new/new_snp",i,".jma.cojo")
  if (file.exists(file_path)==F ) next
  a=fread(file_path)
  a=a[a$p<5e-8 & a$pJ<5e-8,]
  data=rbind(data,a)
  print(i)
}

data$Cytoband = cytobandmap(data$Chr,data$bp)

####+rsid                   
gwas_data =fread("/home/sshen/EPYC/DiskA/LC_GWAS/LC_META_shen/01META_All.txt.gz")
gwas_data$Allele1=toupper(gwas_data$Allele1)
gwas_data$Allele2=toupper(gwas_data$Allele2)

data_id=data$SNP
c=gwas_data[gwas_data$MarkerName %in% data_id,]
df=merge(data,c,by.x="SNP",by.y="MarkerName")
final=df[,c("chr","pos","rsid","SNP","Cytoband","Allele1","Allele2","Freq1","Effect","StdErr","Pvalue","Direction","HetISq","HetChiSq","HetDf","HetPVal","EffectARE","StdErrARE","PvalueARE","pJ","tausq","StdErrMRE","PvalueMRE")]                              
colnames(final)=c("CHROM","GENPOS","rsid","SNP","Cytoband","effect_allele","other_allele","Freq1","Effect","StdErr","Pvalue","Direction","HetISq","HetChiSq","HetDf","HetPVal","EffectARE","StdErrARE","PvalueARE","P_conditional","tausq","StdErrMRE","PvalueMRE")
final=final[order(final$CHROM,final$GENPOS),]   
write.xlsx(final,paste0("meta_sig_indep.xlsx"),rowNames = F,quote=F,sep="\t")



##Annotation
load("/home/sshen/public/resource/G1000_s_hg38.RData")
G1000_s$label = paste(G1000_s$chr,G1000_s$hg38,sep=":")
snplist=final$SNP 
a=G1000_s[G1000_s$label %in% snplist,]
final$gene_name <- a$gene_name[match(final$SNP, a$label)]
final$Consequence <- a$Consequence[match(final$SNP, a$label)]
#final=merge(final,a,by.x="SNP",by.y="label",all.x = TRUE)  
final$rsid=ifelse(final$rsid=="",final$SNP,final$rsid)

#new cytoband
ng <- fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/ng.txt"),header=F)
colnames(ng) <- c("SNP")

catalog <- read.table("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/gwas-association-downloaded_2025-10-20-MONDO_0008903-withChildTraits.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#catalog=catalog[catalog$'P.VALUE'<5e-8,]             
catalog=catalog[,c("CHR_ID","CHR_POS","SNPS","REGION")]
catalog$SNP=paste0(catalog$CHR_ID,":",catalog$CHR_POS)

c1=catalog[catalog$CHR_ID=="",]
c1$SNP=gsub("chr", "",c1$SNPS)
c1 <- c1[!grepl("^rs", c1$SNP), ]

catalog=catalog[catalog$CHR_ID %in% 1:22,]
catalog=rbind(catalog,c1)
x=data.frame(SNP=catalog[,"SNP"])
x=unique(x) #1412
catalog_ng_snp=rbind(ng,x)
catalog_ng_snp=unique(catalog_ng_snp)  #2333

catalog_ng_snp$chr=sapply(catalog_ng_snp$SNP,function(x) strsplit(x,":")[[1]][1])
catalog_ng_snp$pos=sapply(catalog_ng_snp$SNP,function(x) strsplit(x,":")[[1]][2])
catalog_ng_snp$chr=as.numeric(catalog_ng_snp$chr)   
catalog_ng_snp$pos=as.numeric(catalog_ng_snp$pos)   
catalog_ng_snp$cytoband = cytobandmap(catalog_ng_snp$chr,catalog_ng_snp$pos)

##region
region=setdiff(final$Cytoband,catalog_ng_snp$cytoband)
region1=final[final$Cytoband %in% region,]

write.xlsx(region1,paste0("unreported_region.xlsx"),rowNames = F,quote=F,sep="\t")
write.xlsx(final,paste0("cond_snp_ref_ng_catalog.xlsx"),rowNames = F,quote=F,sep="\t")

