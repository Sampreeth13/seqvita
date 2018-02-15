suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
fl <- ("hg19_refGene.gtf")
txdb <- makeTxDbFromGFF(file=fl,format="gtf")
fl <- system.file("extdata", "SRR1518358_default.vcf",package="VariantAnnotation")
vcf <- readVcf('Variants_20_30.vcf', "hg19")

fa <- open(FaFile('human_gik_v37.fasta'))
coding1 <- predictCoding(vcf,txdb,fa)
final <- data.frame(seqnames=seqnames(coding1),start=start(coding1),ref=coding1$REF,alt=coding1$ALT,type=coding1$CONSEQUENCE,gene_id=coding1$GENEID)
test <- final
final <- test[,c("seqnames","start","ref","alt.value","type","gene_id")]
ans <- final %>% group_by(seqnames,start,ref,alt.value,type) %>% summarise(gene_id = toString(gene_id)) %>% as.data.frame
test <-  rfPred_scores(variant_list = ans[1:4],data="all_chr_rfPred.txtz",index="all_chr_rfPred.txtz.tbi",all.col = TRUE)
write.table(ans, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)

loc_all <- locateVariants(vcf, txdb, AllVariants())
tets <- unique(loc_all)
final1 <- data.frame(seqnames=seqnames(tets),start=start(tets),location=tets$LOCATION,gene_id=tets$GENEID,QUERYID=tets$QUERYID)
region <- IntergenicVariants(upstream=70000, downstream=70000)
loc_int <- locateVariants(vcf, txdb, region)
test <- mcols(loc_int)[c("LOCATION", "PRECEDEID", "FOLLOWID","QUERYID")]
xxx = CharacterList(character(0))
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
a <- merge(x=final1,y=test,by = "QUERYID",all.x = TRUE)
v <- rowRanges(vcf)
var <- data.frame(seqnames=seqnames(v),start=start(v),REF=v$REF,ALT=v$ALT)
var1 <- var[,c("seqnames","start","REF","ALT.value")]
b <- merge(x=a,y=var1,by.x=c("seqnames", "start"),by.y=c("seqnames","start"),all.x = TRUE)
new <- b[order(b$QUERYID),]
c <- new[,c("seqnames","start","REF","ALT.value","location","gene_id","PRECEDEID","FOLLOWID")]
d <- merge(x=b,y=ans,by.x=c("seqnames","start"),by.y=c("seqnames","start"),all.x=TRUE)
e <- d[,c("seqnames","start","REF","ALT.value","location","gene_id.x","PRECEDEID","FOLLOWID","type")]
write.table(e, file="Ann.txt", quote=F, sep="\t", row.names=F, col.names=F)
