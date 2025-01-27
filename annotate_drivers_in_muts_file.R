patient="PD48372i"

snvs<-read.table(paste0(patient,"_ndp_assigned_snvs.csv"),sep=",",header=T)
indels<-read.table(paste0(patient,"_ndp_assigned_indels.csv"),sep=",",header=T)

snvs$Gene <- as.character(snvs$Gene)
indels$Gene <- as.character(indels$Gene)
snvs$Driver <- ""

indels[indels$Pos>24305187 & indels$Pos<24311430 & indels$Chrom=="chr14","Gene"] <- "CIDEB"
indels$Driver <- ""
indels[indels$Gene %in% c("FOXO1","GPAM","TNRC6B","CIDEB","INSR","FASN","KLF15","CDKN1B","ALB","ACVR2A","CYP2E1") & indels$Effect %in% c("frameshift","5prime_UTR_ess_splice","-") & indels$vaf>0,"Driver"] <- "DRIVER"
driver_indels <- indels[indels$Driver=="DRIVER",]

snvs[snvs$Gene  %in% c("FOXO1","GPAM","TNRC6B","CIDEB","INSR","FASN","KLF15","CDKN1B","ALB","ACVR2A","CYP2E1") & snvs$Effect %in% c("ess_splice","missense","nonsense","stop_loss","start_loss") & snvs$vaf>0,"Driver"] <- "DRIVER"

muts <- rbind(snvs,driver_indels)

# ! FOXO1 hitspots: manually assign cluster based on the heatmap

write.table(muts,file=paste0(patient,"_ndp_assigned_muts.csv"),sep=",",col.names=T,row.names=F,quote=F)

