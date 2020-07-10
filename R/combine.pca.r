library(data.table)
library(ggplot2)

source("./tsne.util.r")

load("./rna/TCGA.rna.filter1.allgene.pca.rda")
pca.rna = pca
rm(list=c("pca"))
patient.rna = substr(rownames(pca.rna$x), 1, 12)

load("./mirna/TCGA.mirna.filter1.allgene.pca.rda")
pca.mirna = pca2
rm(list=c("pca2"))
patient.mirna = substr(rownames(pca.mirna$x), 1, 12)

load("./meth/TCGA.met.pca.perm100.rda")
pca.met = pca
rm(list=c("pca"))
patient.met = substr(rownames(pca.met$x), 1, 12)

load("./cnv/TCGA.cnv.1m.pca.rda")
pca.cnv = pca2
rm(list=c("pca2"))
patient.cnv = substr(rownames(pca.cnv$x), 1, 12)


# Merge PCs
use = 200
types = c("rna", "met", "cnv", "mirna")
weight = data.frame(type = types)
weight$weight = c(3, 3, 1, 1)

varexp = rbind( cbind(var=pca.met$varexp, type="met", index=1:length(pca.met$varexp)), 
				cbind(var=pca.rna$varexp, type="rna", index=1:length(pca.rna$varexp)), 
				cbind(var=pca.mirna$varexp, type="mirna", index=1:length(pca.mirna$varexp)), 
				cbind(var=pca.cnv$varexp, type="cnv", index=1:length(pca.cnv$varexp)))
varexp = data.table(varexp)
varexp$var = as.numeric(varexp$var)
varexp$index = as.integer(varexp$index)
# put on weight
varexp = merge(varexp, weight, by="type")
varexp$adj.var = varexp$var * varexp$weight

# order by variation explained
varexp = varexp[order(-varexp$adj.var), ]
# pick top ones
varexp = varexp[1:use, ]
temp =  aggregate(varexp$var, by=list(varexp$type), FUN=length)
names(temp) = c("type", "numPC")
weight = merge(weight, temp, by="type")
temp =  aggregate(varexp$var, by=list(varexp$type), FUN=sum)
names(temp) = c("type", "varexp")
weight = merge(weight, temp, by="type")
rownames(weight) = weight$type

# common patient
patients = Reduce(intersect, list(patient.met, patient.rna, patient.mirna, patient.cnv))
weight = weight[match(types, weight$type), ]
weight$patient = list(patient.rna, patient.met, patient.cnv, patient.mirna)
weight$patient.index = lapply(weight$patient, function(x) match(patients, x))


weight$sumvar = c(	sum(pca.rna$sdev[1: weight$numPC[weight$type=="rna"]] ^2), 
					sum(pca.met$sdev[1: weight$numPC[weight$type=="met"]] ^2), 
					sum(pca.cnv$sdev[1: weight$numPC[weight$type=="cnv"]] ^2), 
					sum(pca.mirna$sdev[1: weight$numPC[weight$type=="mirna"]] ^2)
					)

weight$adj.factor = sqrt(weight$weight / weight$sumvar)

PC = cbind(	pca.rna$x[unlist(weight$patient.index[weight$type=="rna"]), 1:weight$numPC[weight$type=="rna"]] * weight$adj.factor[weight$type=="rna"], 
			pca.met$x[unlist(weight$patient.index[weight$type=="met"]), 1:weight$numPC[weight$type=="met"]] * weight$adj.factor[weight$type=="met"], 
			pca.cnv$x[unlist(weight$patient.index[weight$type=="cnv"]), 1:weight$numPC[weight$type=="cnv"]] * weight$adj.factor[weight$type=="cnv"], 
			pca.mirna$x[unlist(weight$patient.index[weight$type=="mirna"]), 1:weight$numPC[weight$type=="mirna"]] * weight$adj.factor[weight$type=="mirna"]
			)

rownames(PC) = patients
save(PC, file="combined.PC.20161108_3_3_1_1.rda")

library(Rtsne)
ptm <- proc.time()
tsne3 = tsne(PC, permute=1000, use=use, dims=3)
proc.time() - ptm
# set expected contribution 


# variable name "count"

# replace sample names to barcode
meta = fread("./meta/mirna.meta.txt")
m = match(rownames(tsne3), substr(meta$barcode, 1, 12))
tsne3$project = meta$disease[m]
save(tsne3, file="combined.tsne3.PC200.dim3.1_1_1_1.perm1000.rda")


png("combined.tsne2.3311.PC200.perm1000.png", width=1280, height=800)
p = ggplot(tsne2, aes(x=feature_1, y=feature_2, col=project, shape = project)) + 
	geom_point(size = 2) + 
	scale_shape_manual(values=c(rep(c(0, 2, 8, 16, 17), 6), 0, 2, 8) )
p
dev.off()


