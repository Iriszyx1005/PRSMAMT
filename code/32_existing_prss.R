suppressMessages(library(data.table))
suppressMessages(library(dplyr));require(BEDMatrix)

args = commandArgs(trailingOnly = TRUE)
plink_dir = args[1]   #get the first parameter
trait_dir = args[2]    #get the second parameter
trait = args[3]  
pop = args[4]    

setwd(paste0(trait_dir,"/",trait,"/",pop))

a=dir(paste0(trait_dir,"/",trait))
dta=a[grepl("_hmPOS_GRCh38",a)]

x=fread(paste0(trait_dir,"/",trait,"/",dta))
x=x[x$hm_chr %in% c(1:22),]
x$hm_chr=as.numeric(x$hm_chr)
x=x[order(x$hm_chr,x$hm_pos),]

prs=NULL;beta_all = NULL
for(j in 1:22){
  pgs=x[x$hm_chr==j,]
  if(nrow(pgs) > 0) {
    bim=fread(paste0(plink_dir,"/chr",j,".bim"))
    pos = intersect(pgs$hm_pos,bim$V4)
    if(length(pos) > 0){
    bim2 = bim[match(pos,bim$V4),]
    pgs2 = pgs[match(pos,pgs$hm_pos),]

    for(i in 1:nrow(bim2)){
        if(pgs2$effect_allele[i] == bim2$V5[i]) next
        if(pgs2$effect_allele[i] != bim2$V5[i]){
          if(pgs2$effect_allele[i] == "A") pgs2$effect_allele[i] = "T"
          else if(pgs2$effect_allele[i] == "T") pgs2$effect_allele[i] = "A"
          else if(pgs2$effect_allele[i] == "C") pgs2$effect_allele[i] = "G"
          else if(pgs2$effect_allele[i] == "G") pgs2$effect_allele[i] = "C"
        }
        if(pgs2$effect_allele[i] != bim2$V5[i]) pgs2$effect_weight[i] = -pgs2$effect_weight[i]
    } 
    score=data.frame(ID=bim2$V2,chr = bim2$V1,pos = bim2$V4,effect_allele=bim2$V5,other_allele=bim2$V6,beta=pgs2$effect_weight)
    write.table(bim2$V2,paste0("snplist.txt"),row.names = F,col.names = F,quote=F,sep="\t")
    system(paste0("plink2 --bfile ",plink_dir,"/chr",j," --extract snplist.txt --make-bed --out temp"))
    bed = as.matrix(BEDMatrix("temp.bed",simple_names=T))
    score1 = bed%*%pgs2$effect_weight
    prs=cbind(prs,score1)
    beta_all = rbind(beta_all,score)
    print(j)
  }
}
}

prs = as.data.frame(prs)
prs$prs = rowSums(prs)
write.table(prs,paste0(trait,"_prs.txt"),quote=F,sep="\t")
write.table(beta_all,paste0(trait,"_beta.txt"),row.names = F,quote=F,sep="\t")
