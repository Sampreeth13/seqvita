suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rfPred))
arg <- commandArgs(trailingOnly = TRUE)
txdb <- suppressWarnings(makeTxDbFromGFF(file=arg[1],format="gtf"))
vcf <- suppressWarnings(readVcf(arg[2], "hg19"))
if(identical(intersect(seqlevels(vcf),seqlevels(txdb)),character(0))){
  if(identical(intersect(seqlevels(vcf),"1"),character(0))){
    txdb <- renameSeqlevels(txdb, paste0("chr", seqlevels(txdb)))
  }
  else{
    vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))
  }
}
fa <- open(FaFile(arg[3]))
coding1 <- suppressWarnings(predictCoding(vcf,txdb,fa))
final <- data.frame(seqnames=seqnames(coding1),start=start(coding1),ref=coding1$REF,alt=coding1$ALT,type=coding1$CONSEQUENCE,gene_id=coding1$GENEID)
if(nrow(final) == 0){
  final$type <- numeric(nrow(final))
}
test <- final
final <- test[,c("seqnames","start","ref","alt.value","type","gene_id")]
ans <- unique(final)
genes = fread("../Liver_Cancer/ENS_ID_to_GENE")
colnames(genes) = c("ENSID","gene_name")
ans = merge(x=ans,y=genes,by.x=c("gene_id"),by.y=c("ENSID"),all.x=TRUE)
ans = ans[,c(2:ncol(ans))]
ans= ans[order(ans$seqnames,ans$start),]
ans$seqnames<-gsub("chr(*)","\\1",ans$seqnames)
if(nrow(ans) > 0){
  test <-  rfPred_scores(variant_list = ans[1:4],data=arg[4],index=arg[5],all.col = TRUE)
  p <- merge(x=ans,y=test,by.x=c("seqnames", "start","gene_name"),by.y=c("chromosome","position_hg19","genename"),all.x=TRUE)
  q<-p[with(p,order(start)), ]
  r<-q[,c("seqnames","start","ref","alt.value","type","gene_name","SIFT_score","Polyphen2_score","MutationTaster_score","PhyloP_score","LRT_score")]
  r$seqnames <- paste("chr",r$seqnames,sep="")
}else{
  r = data.frame(seqnames=character(),start=character(),ref=character(),alt.value=character(),type=character(),gene_name=character(),SIFT_score=character(),Polyphen2_score=character(),MutationTaster_score=character(),PhyloP_score=character(),"LRT_score"=character())
}
loc_all <- suppressWarnings(locateVariants(vcf, txdb, AllVariants()))
final1 <- data.frame(seqnames=seqnames(loc_all),start=start(loc_all),location=loc_all$LOCATION,gene_id=loc_all$GENEID,QUERYID=loc_all$QUERYID)
final1<- unique(final1)
test1 = merge(x=final1,y=genes,by.x=c("gene_id"),by.y=c("ENSID"),all.x=TRUE)
final1 = test1[,c("seqnames","start","location","gene_name","QUERYID")]
final1=final1[order(final1$seqnames,final1$start),]
#final1 <- final1 %>% group_by(seqnames,start,gene_id,QUERYID) %>% summarise(location = toString(location)) %>% as.data.frame
region <- IntergenicVariants(upstream=70000, downstream=70000)
loc_int <- suppressWarnings(locateVariants(vcf, txdb, region))
test <- mcols(loc_int)[c("LOCATION", "PRECEDEID", "FOLLOWID","QUERYID")]
xxx = CharacterList(character(0))
if(nrow(test)>0){
  for(i in 1:nrow(test))
  {
    if(identical(test[i,"PRECEDEID"],xxx)){
      test[i,"PRECEDEID"] = "NA"
    }
    else{
      test[i,"PRECEDEID"]=paste(test[i,"PRECEDEID"],collapse=",")
    }
    if(identical(test[i,"FOLLOWID"],xxx)){
      test[i,"FOLLOWID"] = "NA"
    }
    else{
      test[i,"FOLLOWID"]=paste(test[i,"FOLLOWID"],collapse=",")
    }
  }
}
indx <- setNames(as.character(genes$gene_name), genes$ENSID)
test$PRECEDEID <- vapply(strsplit(as.character(test$PRECEDEID),','),function(x) paste(indx[x],collapse=','),character(1L))
test$FOLLOWID <- vapply(strsplit(as.character(test$FOLLOWID),','),function(x) paste(indx[x],collapse=','),character(1L))
a <- merge(x=final1,y=test,by = "QUERYID",all.x = TRUE)
v <- rowRanges(vcf)
var <- data.frame(seqnames=seqnames(v),start=start(v),REF=v$REF,ALT=v$ALT)
var1 <- var[,c("seqnames","start","REF","ALT.value")]
b <- merge(x=a,y=var1,by.x=c("seqnames", "start"),by.y=c("seqnames","start"),all.x = TRUE)
new <- b[order(b$QUERYID),]
c <- new[,c("seqnames","start","REF","ALT.value","location","gene_name","PRECEDEID","FOLLOWID")]
d <- merge(x=c,y=r,by.x=c("seqnames","start","gene_name"),by.y=c("seqnames","start","gene_name"),all.x=TRUE)
e <- d[,c("seqnames","start","REF","ALT.value","location","gene_name","PRECEDEID","FOLLOWID","type","SIFT_score","Polyphen2_score","MutationTaster_score","PhyloP_score","LRT_score")]
e <- unique(e)
suppressPackageStartupMessages(library(tidyr))
kgp = fread("../Liver_Cancer/clinVar_variants_hg19")
test = kgp %>% separate(name, into = c("ref","alt"), sep='>')
clinvar = test[,c("#chrom","chromEnd","ref","alt","clinSign","phenotypeList")]
e = merge(x=e,y=clinvar,by.x=c("seqnames","start","REF","ALT.value"),by.y=c("#chrom","chromEnd","ref","alt"),all.x=TRUE)
cosmic = fread("../Liver_Cancer/cosmic_with_chr.vcf")
cosmic = cosmic[,c("#CHROM","POS","ID","REF","ALT")]
colnames(cosmic) = c("CHROM","POS","Cosmic_Ids","REF","ALT")
map = e[,c("seqnames","start","REF","ALT.value")]
map = unique(map)
map = merge(x=map,y=cosmic,by.x=c("seqnames","start","REF","ALT.value"),by.y=c("CHROM","POS","REF","ALT"),all.x=TRUE)
map = as.data.frame(map)
map <- map %>% group_by(seqnames,start,REF,ALT.value) %>% summarise(Cosmic_Ids = toString(Cosmic_Ids)) %>% as.data.frame
e = merge(x=e,y=map,by.x=c("seqnames","start","REF","ALT.value"),by.y=c("seqnames","start","REF","ALT.value"),all.x=TRUE)
omim  = read.csv("../Liver_Cancer/hg19_OMIM",sep="\t",header=FALSE)
colnames(omim) = c("CHROM","start","end","OMIM_ID")
omim = omim[,c("CHROM","start","OMIM_ID")]
e = merge(x=e,y=omim,by.x = c("seqnames","start"),by.y=c("CHROM","start"),all.x=TRUE)
decipher  = fread("../Liver_Cancer/hg19_DECIPHER")
colnames(decipher) = c("CHROM","start","end","Decipher_values")
decipher = decipher %>% separate(Decipher_values, into = c("Gene_name","Decipher_val_1","Decipher_val_2"), sep='[|]')
map_decipher = e[,c("seqnames","start","REF","ALT.value")]
map_decipher = unique(map_decipher)
ivar = with(map_decipher, GRanges(seqnames, IRanges(start,width=1)))
ideci = with(decipher, GRanges(CHROM, IRanges(start,end)))
olaps = suppressWarnings(findOverlaps(ivar,ideci))
decipher_results = cbind(map_decipher[queryHits(olaps),],decipher[subjectHits(olaps),])
map_decipher = as.data.frame(decipher_results[,c("seqnames","start","REF","ALT.value","Gene_name","Decipher_val_1","Decipher_val_2")])
e = merge(x=e,y=map_decipher,by.x=c("seqnames","start","REF","ALT.value","gene_name"),by.y=c("seqnames","start","REF","ALT.value","Gene_name"),all.x=TRUE)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
map_dbsnp =  e[,c("seqnames","start","REF","ALT.value")]
map_dbsnp = unique(map_dbsnp)
map_dbsnp$seqnames<-gsub("chr(*)","\\1",map_dbsnp$seqnames)
ivar = with(map_dbsnp, GRanges(seqnames, IRanges(start,width=1)))
test=suppressWarnings(snpsByOverlaps(snps,ivar,drop.rs.prefix=FALSE))
dbsnp = data.frame(seqnames=seqnames(test),start=pos(test),dbsnp_id=test$RefSNP_id)
dbsnp$seqnames <- paste("chr",dbsnp$seqnames,sep="")
e = merge(x=e,y=dbsnp,by.x=c("seqnames","start"),by.y=c("seqnames","start"),all.x=TRUE)
pgkb = fread("../Liver_Cancer/pharmGKB.tsv")
pgkb = pgkb[,c(2:5)]
colnames(pgkb) = c("gene","pgkb_type","pgkb_level","pgkb_chemicals")
e = merge(x=e,y=pgkb,by.x=c("gene_name"),by.y=c("gene"),all.x=TRUE)
e = e[,c(2:5,1,6:ncol(e))]
e=e[order(e$seqnames,e$start),]
f = as.data.frame(e)
f %>% mutate_if(is.factor,as.character) -> f
f[is.na(f)]<- "."
if(nrow(f)>0){
  for(i in 1:nrow(f))
  { 
    if(identical(f[i,"PRECEDEID"],"NA")){
      f[i,"PRECEDEID"] = "."
    }
    if(identical(f[i,"FOLLOWID"],"NA")){
      f[i,"FOLLOWID"] = "."
    }
    functional_score = suppressWarnings(as.numeric(f[i,"SIFT_score"]) + as.numeric(f[i,"Polyphen2_score"]) + as.numeric(f[i,"MutationTaster_score"]) + as.numeric(f[i,"PhyloP_score"]) + as.numeric(f[i,"LRT_score"]))
    if(identical(f[i,"SIFT_score"],"NA")){
      f[i,"SIFT_score"] = "."
    }
    if(identical(f[i,"Polyphen2_score"],"NA")){
      f[i,"Polyphen2_score"] = "."
    }
    if(identical(f[i,"MutationTaster_score"],"NA")){
      f[i,"MutationTaster_score"] = "."
    }
    if(identical(f[i,"PhyloP_score"],"NA")){
      f[i,"PhyloP_score"] = "."
    }
    if(identical(f[i,"LRT_score"],"NA")){
      f[i,"LRT_score"] = "."
    }
    functional_priority = 0
    if(is.na(functional_score)){
      functional_score = 0
    }
    if(functional_score < 1.5){
      functional_priority = 1
    }
    else if(functional_score < 3.5){
      functional_priority = 2
    }
    else{
      functional_priority = 3
    }
    clinical_priority = 0
    if(!identical(f[i,"clinSign"],".")){
      clinical_priority = clinical_priority + 1
    }
    if(!identical(f[i,"OMIM_ID"],".")){
      clinical_priority = clinical_priority + 1
    }
    if(!identical(f[i,"Cosmic_Ids"],"NA")){
      clinical_priority = clinical_priority + 1
    }
    else{
      f[i,"Cosmic_Ids"] = "."
    }
    if(!identical(f[i,"Decipher_val_1"],".")){
      if(clinical_priority < 3){
        clinical_priority = clinical_priority + 1
      }
    }
    if(identical(clinical_priority,0)){
      clinical_priority = 1
    } 
    drug_priority = 1  
    if(!identical(f[i,"pgkb_type"],".")){
      drug_priority = 3
    }
    total_score = functional_priority + clinical_priority + drug_priority
    priority = ""
    if(total_score < 4){
      priority = "Low"
    }
    else if(total_score < 7){
      priority = "Medium"
    }
    else{
      priority = "High"
    }
    f[i,"Priority"] = priority
  }
}
write.table(f, file=arg[6], quote=F, sep="\t", row.names=F, col.names=T)