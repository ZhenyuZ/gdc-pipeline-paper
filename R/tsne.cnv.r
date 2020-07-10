library(data.table)


source("./tsne.util.r")

cnv.file = "TCGA.cnv.window_1M.txt"
cnv = read.table(cnv.file, h=T, stringsAsFactors=F, row.names=1, check.names=F)

# aliquot selection
aliquot = colnames(cnv)
aliquot2 = SelectAliquot(aliquot, my.type = c("01", "03"))
patient = substr(aliquot2, 1, 12)
m = match(patient, substr(aliquot, 1, 12))
cnv = cnv[, m]
colnames(cnv) = patient

# replace sample names to barcode
meta = fread("../meta/snp6.sample.txt", h=T)
meta = meta[, c("patient_barcode", "disease"), with=F]
meta$disease = toupper(meta$disease)
colnames(meta) = c("patient", "project")
m = match(colnames(cnv2), meta$patient)
meta = meta[m, ]


# QC
# rnaseq data do not have NA. Because this is accross multiple project, zero.ratio set to 0.99
cnv2 = CompletenessFilter(cnv, feature.missing = 0.1, sample.mising = 1, dup.removal = T)
if(sum(is.na(cnv2)) > 0) {
	cnv2 = ImputeNA(cnv2)
}
save(cnv2, file="cnv2.1M.rda")

# pca
load("cnv2.1M.rda")
ptm <- proc.time()
pca2 = pca(cnv2)
proc.time() - ptm
save(pca2, file="TCGA.cnv.1m.pca.rda")

# t-SNE
library(Rtsne)
ptm <- proc.time()
tsne2 = tsne(pca2$x, permute=1, use=50)
proc.time() - ptm
# save(tsne2, file="TCGA.rna.allgene.tsne.rda")

pdata2 = tsne2
pdata2$project = meta$project
write.table(pdata2, "TCGA.cnv.10M.per10.pdata.txt", col.names=T, row.names=F, sep="\t", quote=F)



png("TCGA.cnv.10M.per10.png", width=1200, height=800)
p = ggplot(pdata2, aes(x=feature_1, y=feature_2, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


