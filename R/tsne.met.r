library(data.table)
source("./tsne.util.r")

load("TCGA.methylation.beta.sd0.2.rda")

# pca
ptm <- proc.time()
pca = pca(met2)
proc.time() - ptm
save(pca, file="TCGA.met.pca.perm100.rda")

# t-SNE
library(Rtsne)
ptm <- proc.time()
tsne = tsne(pca$x, permute=100, use=50)
proc.time() - ptm
# about 6 hours seed = 25
save(tsne, file="TCGA.met.tsne.perm100.rda")


pdata = annotate.tsne(tsne, aliquot)
write.table(pdata, "TCGA.met.tsne.pdata.perm100.txt", col.names=T, row.names=F, sep="\t", quote=F)

# t-SNE of randome sampleing
# pca3 = pca(met3)
# tsne3 = tsne(pca3$x)
# tsne3.pdata = annotate.tsne(tsne3, aliquot)
# tsne2.pdata = annotate.tsne(tsne2, aliquot)
# p3 = ggplot(tsne3.pdata, aes(x=feature_1, y=feature_2, col=project)) + geom_point()

save(tsne2.pdata, tsne3.pdata, file="tsne.pdata.rda")

pdata$sh = rep(c(0, 2, 4), 11)

png("TCGA.met.tsne.pdata.perm100.png", width=1200, height=800)
p = ggplot(pdata, aes(x=feature_1, y=feature_2, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()

