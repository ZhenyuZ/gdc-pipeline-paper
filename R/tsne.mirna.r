library(data.table)

source("./tsne.util.r")

load("mirna.count.rda")
# variable name "count"

# replace sample names to barcode
meta = fread("../meta/mirna.meta.txt")
m = match(colnames(count), meta$aligned_file_gdcid)
aliquot = meta$barcode[m]
colnames(count) = aliquot

# aliquot selection
aliquot2 = SelectAliquot(aliquot, my.type = c("01", "03"))
m = match(aliquot2, colnames(count))
count = count[, m]

# QC
# rnaseq data do not have NA. Because this is accross multiple project, zero.ratio set to 0.99
count2 = CompletenessFilter(count, feature.missing = 1, sample.mising = 1, zero.ratio = 0.99, dup.removal = T)
# low count expression filter
low.expr = apply(count2, 1, function(x) sum(x <= 1) / ncol(count2))
count2 = count2[which(low.expr < 0.99, ), ]
save(count2, file="mirna.count2.rda")

# make log(tpm)
total = apply(count2, 2, sum)
count2[which(count2 == 0)] = 1 
tpm2 = t(t(count2)/total)
tpm2 = log2(tpm2)



# pca
ptm <- proc.time()
pca2 = pca(tpm2)
proc.time() - ptm
save(pca2, file="TCGA.mirna.filter1.allgene.pca.rda")

# t-SNE
library(Rtsne)
ptm <- proc.time()
tsne2 = tsne(pca2$x, permute=100, use=50)
proc.time() - ptm
save(tsne2, file="TCGA.mirna.allgene.tsne.rda")


aliquot = fread("disease_aliquot.txt")

meta$aliquot = meta$barcode
pdata2 = annotate.tsne(tsne2, meta)
pdata2 = pdata2[, c("aliquot", "feature_1", "feature_2", "disease", "center")]
colnames(pdata2)[4] = "project"
write.table(pdata2, "TCGA.mirna.allgene.pdata.txt", col.names=T, row.names=F, sep="\t", quote=F)



png("TCGA.mirna.allgene.pdata.perm100.png", width=1200, height=800)
p = ggplot(pdata2, aes(x=feature_1, y=feature_2, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


