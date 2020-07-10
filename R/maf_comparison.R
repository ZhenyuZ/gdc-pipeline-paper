# Code to generate varies mutation plots in the paper


# function to read MAF file
ReadMAF <- function(maffile, small=T, short_key=T, full_key=F, gz=T) {
  # support gz file
  if (gz) {
    maffile <- paste("zcat", maffile)
  }
  select <- c("Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", 
        "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
        "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "Allele", "FILTER", 
        "GDC_Validation_Status", "MC3_Overlap", "GDC_FILTER")
# class <- c("factor", "factor", "factor", integer", "integer", "factor", "factor", "character", "character",
#     "character", "integer", "integer", "integer", "integer", "integer", "integer", "character", "character", 
#     "factor", "factor", "character")
  # read maf
  # if small is True, only read key columns, and specific colClasses; otherwise, full information with automatic classes
  if(small) {
    maf <- fread(maffile, h=T, select=select, colClasses="character")
    maf$t_depth <- as.integer(maf$t_depth)
    maf$t_ref_count <- as.integer(maf$t_ref_count)
    maf$t_alt_count <- as.integer(maf$t_alt_count)
    maf$n_depth <- as.integer(maf$n_depth)
    maf$n_ref_count <- as.integer(maf$n_ref_count)
    maf$n_alt_count <- as.integer(maf$n_alt_count)
    maf$Start_Position <- as.integer(maf$Start_Position)
    maf$End_Position <- as.integer(maf$End_Position)
    maf$Tumor_Sample_Barcode <- factor(maf$Tumor_Sample_Barcode)
    maf$Matched_Norm_Sample_Barcode <- factor(maf$Matched_Norm_Sample_Barcode)
    maf$Chromosome <- factor(maf$Chromosome)
    maf$Variant_Classification <- factor(maf$Variant_Classification)
    maf$Variant_Type <- factor(maf$Variant_Type)
    maf$GDC_Validation_Status <- factor(maf$GDC_Validation_Status)
    maf$MC3_Overlap <- factor(maf$MC3_Overlap)
  } else {
    maf <- fread(maffile, h=T)
  }
  # add keys
  # short_key: tumor barcode and variant locus only
  # full_key: short_key with tumor sequencing allleles
  if(short_key) {
    maf$short_key <- with(maf, paste(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, sep="|")) 
  }
  if(full_key) {
    maf$full_key <- with(maf, paste(Tumor_Sample_Barcode, Chromosome, Start_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep="|"))
  } 
  return(maf)
}

# function to extract overlapping metrics
getOverlapStat <- function(overlap) {
  pipelines <- colnames(overlap)
  t <- table(apply(overlap, 1, function(x) interaction(x[1], x[2], x[3], x[4])))
  stat <- data.frame(names(t))
  colnames(stat) <- interaction(pipelines[1], pipelines[2], pipelines[3], pipelines[4])
  stat$count <- t
  stat$percente <- t/ nrow(overlap)* 100
  cat("by number of pipelines: \n")
  t <- table(apply(overlap, 1, sum)) / nrow(overlap) * 100
  cat("4 pipelines:", t[names(t)==4], "%\n")
  cat("At least 3 pipelines:", t[names(t)==4] + t[names(t)==3], "%\n")
  cat("At least 2 pipelines:", t[names(t)==4] + t[names(t)==3] + t[names(t)==2], "%\n")
  cat("At least 1 pipelines:", t[names(t)==4] + t[names(t)==3] + t[names(t)==2] + t[names(t)==1], "%\n")
  cat("no pipelines:", t[names(t)==0], "%\n")
  return(stat)
}


library(data.table)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(gridExtra)


##############################################################
# Read All TCGA MAFs from GDC DR 10
##############################################################

# read all MAFs manifest
mafs <- fread("gdc_maf_20170309.txt")

protected.mafs <- mafs[type=="protected"]
public.mafs <- mafs[type=="public"]

# read all public MAF
public <- NULL
for(i in 1: nrow(public.mafs)) {
  my.entry <- public.mafs[i]
  my.maf <- ReadMAF(my.entry$filename)
  my.maf$program <- factor(my.entry$program)
  my.maf$project <- factor(my.entry$project)
  my.maf$pipeline <- factor(my.entry$pipeline)
  public <- rbind(public, my.maf)
}

# read all protected MAF
protected <- NULL
for(i in 1: nrow(protected.mafs)) {
  my.entry <- protected.mafs[i]
  my.maf <- ReadMAF(my.entry$filename)
  my.maf$program <- factor(my.entry$program)
  my.maf$project <- factor(my.entry$project)
  my.maf$pipeline <- factor(my.entry$pipeline)
  protected <- rbind(protected, my.maf)
}
save(protected, file="protected.maf.short.rda")

# fill in depth information in public maf
public.key <- with(public, paste(pipeline, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Chromosome, Start_Position, End_Position, sep="|"))
protected.key <- with(protected, paste(pipeline, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Chromosome, Start_Position, End_Position, sep="|"))
m <- match(public.key, protected.key)
rm(list = c("public.key", "protected.key"))
gc()
public$n_ref_count <- protected$n_ref_count[m]
public$n_alt_count <- protected$n_alt_count[m]
save(public, file="public.maf.short.rda")
# 1.8G memory

# add column in protected to indicate if the pair is in public, and whether the variant is in public 
protected <- protected[program=="TCGA"]
public.pair <- unique(paste(public$pipeline, public$Tumor_Sample_Barcode, public$Matched_Norm_Sample_Barcode, sep="|"))
protected.pair <- paste(protected$pipeline, protected$Tumor_Sample_Barcode, protected$Matched_Norm_Sample_Barcode, sep="|")
in_public <- array("Unknown", nrow(protected))
w <- which(protected.pair %in% public.pair)
in_public[w] <- "False"
in_public[m] <- "True"
in_public <- factor(in_public)
protected$in_public <- in_public
save(protected, file="tcga.protected.maf.short.rda")
# 32G memory


##############################################################
# Comparing to TCGA Validation and collect caller overlaps
##############################################################
# TCGA validation sequences are extracted from old TCGA MAFs from DCC
validation <- fread("TCGA_Validation_Sequence.w.aliquot.GRCh38.txt")
# 127910 validation entry

# sample level validation
mapping <- fread("barcode.project.mapping.tsv")
mapping$sample <- substr(mapping$tumor_barcode, 1, 15)

protected$sample <- substr(protected$Tumor_Sample_Barcode, 1, 15)
common.sample <- Reduce(intersect, list( unique(protected[pipeline=="muse"]$sample),
                    unique(protected[pipeline=="mutect"]$sample),
                    unique(protected[pipeline=="varscan"]$sample),
                    unique(protected[pipeline=="somaticsniper"]$sample), 
                    validation$sample))
protected <- protected[sample %in% common.sample]
validation <- validation[sample %in% common.sample]
validation$project <- mapping$project[match(validation$sample, mapping$sample)]
validation.gr <- with(validation, GRanges(seqnames <- paste(sample, chr, sep="|"), ranges <- IRanges(start <- start, end <- end)))

overlap <- data.table(matrix(0, nrow(validation), 4))
pipelines <- c("muse", "mutect", "somaticsniper", "varscan")
colnames(overlap) <- pipelines
for(i in 1:4) {
  my.pipeline <- pipelines[i]
  maf.gr <- with(protected[pipeline==my.pipeline], GRanges( seqnames <- paste(sample, Chromosome, sep="|"), ranges <- IRanges(start <- Start_Position, end <- End_Position)))
  overlap[, i] <- countOverlaps(validation.gr, maf.gr)
}
overlap <- overlap > 0
validation <- cbind(validation, overlap)
write.table(validation, "sample.validation.txt", col.names=T, row.names=F, sep="\t", quote=F)

stat <- getOverlapStat(overlap)
write.table(stat, "validation.stat", col.names=T, row.names=F, sep="\t", quote=F)

# by number of pipelines:
# 4 pipelines: 71.57072 %
# At least 3 pipelines: 86.18154 %
# At least 2 pipelines: 90.00658 %
# At least 1 pipelines: 96.77509 %
# no pipelines: 3.224913 %
# data used in figure 1


# Collecting similar metrics that also requires tumor normal pair are the same as in the validation data
protected$pair <- paste(protected$Tumor_Sample_Barcode, protected$Matched_Norm_Sample_Barcode, sep="|")
validation$pair <- paste(validation$tumor_sample_barcode, validation$matched_norm_sample_barcode, sep="|")
common.pair <- Reduce(intersect, list( unique(protected[pipeline=="muse"]$pair),
                    unique(protected[pipeline=="mutect"]$pair),
                    unique(protected[pipeline=="varscan"]$pair),
                    unique(protected[pipeline=="somaticsniper"]$pair), 
                    validation$pair))
protected2 <- protected[pair %in% common.pair]

# filter for variants with PASS FILTER and blank GDC_FILTER
protected3 <- protected2[FILTER=="PASS" & GDC_FILTER==""]
protected3.gr <- with(protected3, GRanges(seqnames <- paste(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Chromosome, sep="|"), ranges <- IRanges(start <- Start_Position, end <- End_Position)))
mutect.gr <- protected3.gr[which(protected3$pipeline=="mutect")]
other.gr <- protected3.gr[which(protected3$pipeline!="mutect")]
other.overlap <- countOverlaps(mutect.gr, other.gr)
unique.gr <- c(mutect.gr, other.gr[which(other.overlap==0)])
protected.overlap <- data.table(short_key <- c(protected3[pipeline=="mutect"]$short_key, protected3[pipeline!="mutect"]$short_key[which(other.overlap==0)]))
protected.overlap$mutect <- countOverlaps(unique.gr, protected3.gr[which(protected3$pipeline=="mutect")]) > 0
protected.overlap$muse <- countOverlaps(unique.gr, protected3.gr[which(protected3$pipeline=="muse")]) > 0
protected.overlap$varscan <- countOverlaps(unique.gr, protected3.gr[which(protected3$pipeline=="varscan")]) > 0 
protected.overlap$somaticsniper <- countOverlaps(unique.gr, protected3.gr[which(protected3$pipeline=="somaticsniper")]) > 0
protected.overlap$n.caller <- protected.overlap$muse + protected.overlap$mutect + protected.overlap$somaticsniper + protected.overlap$varscan


short_key <- unique(protected2$short_key)
protected.variants <- data.table(short_key)
protected.variants$mutect <- protected.variants$short_key %in% protected2[pipeline=="mutect"]$short_key


validation2 <- validation[pair %in% common.pair]
validation2.gr <- with(validation2, GRanges(seqnames <- paste(tumor_sample_barcode, matched_norm_sample_barcode, chr, sep="|"), ranges <- IRanges(start <- start, end <- end)))
overlap2 <- data.table(matrix(0, nrow(validation2), 4))
colnames(overlap2) <- pipelines
for(i in 1:4) {
  my.pipeline <- pipelines[i]
  maf.gr <- with(protected2[pipeline==my.pipeline], GRanges( seqnames <- paste(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Chromosome, sep="|"), ranges <- IRanges(start <- Start_Position, end <- End_Position)))
  overlap2[, i] <- countOverlaps(validation2.gr, maf.gr)
}
overlap2 <- overlap2 > 0
validation2 <- cbind(validation2, overlap2)
write.table(validation2, "pair.validation.txt", col.names=T, row.names=F, sep="\t", quote=F)

stat2 <- getOverlapStat(overlap2)
write.table(stat2, "sample.validation.stat2", col.names=T, row.names=F, sep="\t", quote=F)

# by number of pipelines:
# 4 pipelines: 70.56742 %
# At least 3 pipelines: 86.73151 %
# At least 2 pipelines: 91.02027 %
# At least 1 pipelines: 98.44827 %
# no pipelines: 1.551729 %


##############################################################
# calculate recall of TCGA validated variants by caller
##############################################################
recall <-  melt (
        validation %>%
          group_by(sample) %>%
          summarise(
            project <- first(project), 
            MuSE <- sum(muse) /n(),
            MuTect2 <- sum(mutect) /n(), 
            VarScan2 <- sum(varscan) /n(), 
            SomaticSniper <- sum(somaticsniper) /n()
          ),
        measure.vars <- c("MuSE", "VarScan2", "SomaticSniper", "MuTect2"), 
        value.name <- "Recall_Rate", 
        variable.name <- "Pipeline"
      )
# order by average recall rate
average.recall <- recall %>% group_by(project) %>% summarise(rate <- mean(Recall_Rate))
project.order <- average.recall[order(-average.recall$rate), ]$project
pipeline.order <- c("MuTect2", "MuSE", "VarScan2", "SomaticSniper")
recall$project <- factor(recall$project, levels <- project.order)
recall$Pipeline <- factor(recall$Pipeline, levels <- pipeline.order)

# png("called.validation.freq.png", width=1024, height=800)
theme_set(theme_gray(base_size <- 18))
p1 <- ggplot(recall, aes(x=project, y=Recall_Rate, fill=Pipeline)) +
  geom_boxplot(alpha=0.8) +
  theme(legend.position <- c(0.07, 0.12))

#p <- ggplot(recall, aes(x=project, y=Recall_Rate, fill=Pipeline)) +
# geom_boxplot(alpha=0.8) +
# geom_line(data <- recall, aes(x=project, y=Recall_Rate, group <- Pipeline, colour <- Pipeline), size=3, alpha=0.5) +
# theme(legend.position <- c(0.07, 0.12))
# p1
# dev.off()

##############################################################
# calculate recall of TCGA validated variants by caller combinations
##############################################################
validation$n.caller <- validation$muse + validation$mutect + validation$somaticsniper + validation$varscan
sample.stat <-   validation %>% 
          group_by(sample, project) %>% 
            summarise(count=n(), At_Least_1_Caller=sum(n.caller>=1)/count, At_Least_2_Callers=sum(n.caller>=2)/count, At_Least_3_Callers=sum(n.caller>=3)/count, All_4_Callers=sum(n.caller>=4)/count)
write.table(sample.stat, "sample.validation.stat", col.names=T, row.names=F, sep="\t", quote=F)

sample.stat.melt <- melt(select(sample.stat, -count), variable.name="Num_Callers", value.name="Recall_Rate")
sample.stat.melt$project <- factor(sample.stat.melt$project, levels <- project.order)
# png("n.caller.validation.freq.png", width=1024, height=800)
# theme_set(theme_gray(base_size <- 18))
p2 <- ggplot(sample.stat.melt, aes(x=project, y=Recall_Rate, fill=Num_Callers)) +
  geom_boxplot(alpha=0.8) +
  theme(legend.position <- c(0.08, 0.12), axis.title.x=element_blank(), axis.text.x=element_blank()) + scale_x_discrete(position <- "top") 
# p2
# dev.off()
png("validation.freq.png", width=1280, height=1024)
grid.arrange(p1, p2, ncol=1)
dev.off()

################################################
# Make called.validation.freq.png
################################################
# load("validation0613.rda")
# project <- validation[, unique(project), by = sample]
# names(project) <- c("sample", "Project")
# total <- validation[, length(variant), by = sample]
# names(total) <- c("sample", "total")
# muse <- validation[, sum(muse.pass), by = sample]
# names(muse) <- c("sample", "MuSE")
# mutect <- validation[, sum(mutect.pass), by = sample]
# names(mutect) <- c("sample", "MuTect2")
# varscan <- validation[, sum(varscan.pass), by = sample]
# names(varscan) <- c("sample", "VarScan2")
# somaticsniper <- validation[, sum(somaticsniper.pass), by <- sample]
# names(somaticsniper) <- c("sample", "SomaticSniper")

# validation2 <- Reduce(function(...) merge(..., by="sample"), list(project, total, muse, mutect, varscan, somaticsniper))
# validation2$MuSE <- validation2$MuSE / validation2$total
# validation2$MuTect2 <- validation2$MuTect2 / validation2$total
# validation2$VarScan2 <- validation2$VarScan2 / validation2$total
# validation2$SomaticSniper <- validation2$SomaticSniper / validation2$total
# validation3= melt(validation2, measure.vars <- c("MuSE", "VarScan2", "SomaticSniper", "MuTect2"), value.name ="Proportion_Recalled")

# save(validation3, file="validation3.0613.rda")
# colnames(validation3) <- c("sample", "Project", "total", "Pipeline", "Proportion_Recalled")
validation3 < fread("sample.validation3.tsv")

validation4 <- read.table("called.validation.freq.txt")
validation5 <- data.frame(t(validation4[5:8, ]))
colnames(validation5) <- c("MuSE", "MuTect2", "VarScan2", "SomaticSniper")
validation5$Project= rownames(validation5)
validation6 <- melt(validation5, measure.vars <- c("MuSE", "VarScan2", "SomaticSniper", "MuTect2"), value.name ="Proportion_Recalled", variable.name <- "Pipeline")

# sort by order
validation6 <- data.table(validation6)
temp <- validation6[, mean(Proportion_Recalled), by=Project]
levels <- temp$Project[order(-temp$V1)]

validation3$Project <- factor(validation3$Project, levels <- levels)
validation6$Project <- factor(validation6$Project, levels <- levels)


png("called.validation.freq.png", width=1024, height=800)
theme_set(theme_gray(base_size <- 18))
p <- ggplot(validation3, aes(x=Project, y=Proportion_Recalled, fill=Pipeline)) +
  geom_boxplot(alpha=0.8) +
  geom_line(data <- validation6, aes(x=Project, y=Proportion_Recalled, group <- Pipeline, colour <- Pipeline), size=3, alpha=0.5) +
  theme(legend.position <- c(0.07, 0.12))
p
dev.off()

################################################
# Make called.validation.by.num.caller.png
################################################
validation7 <- validation[, c("project", "sample", "pair", "variant"), with=F]
validation7$pipeline.called <- with(validation, mutect.pass + somaticsniper.pass + somaticsniper.pass + muse.pass)

project <- validation7[, unique(project), by <- sample]
names(project) <- c("sample", "Project")
temp0 <- validation7[, sum(pipeline.called >= 0), by="sample"]
names(temp0) <- c("sample", "total")
temp1 <- validation7[, sum(pipeline.called >= 1), by="sample"]
names(temp1) <- c("sample", "at_least_1_pipeline")
temp2 <- validation7[, sum(pipeline.called >= 2), by="sample"]
names(temp2) <- c("sample", "at_least_2_pipelines")
temp3 <- validation7[, sum(pipeline.called >= 3), by="sample"]
names(temp3) <- c("sample", "at_least_3_pipelines")
temp4 <- validation7[, sum(pipeline.called == 4), by="sample"]
names(temp4) <- c("sample", "all_4_pipelines")
validation8 <- Reduce(function(...) merge(..., by="sample"), list(project, temp0, temp1, temp2, temp3, temp4))


validation8$at_least_1_pipeline <-  validation8$at_least_1_pipeline / validation8$total
validation8$at_least_2_pipelines <-  validation8$at_least_2_pipelines / validation8$total
validation8$at_least_3_pipelines <-  validation8$at_least_3_pipelines / validation8$total
validation8$all_4_pipelines <-  validation8$all_4_pipelines / validation8$total

validation8$Project <- factor(validation8$Project, levels <- c("UCS", "BRCA", "READ", "OV", "COAD", "SARC", "CESC", "BLCA", "KIRC", "THYM", "PAAD", "LAML"))
validation8 <- validation8[, -"total", with=F]
validation9 <- melt(validation8, value.name ="Proportion_Recalled", variable.name <- "Category")


temp0 <- validation7[, sum(pipeline.called >= 0), by="project"]
names(temp0) <- c("project", "total")
temp1 <- validation7[, sum(pipeline.called >= 1), by="project"]
names(temp1) <- c("project", "at_least_1_pipeline")
temp2 <- validation7[, sum(pipeline.called >= 2), by="project"]
names(temp2) <- c("project", "at_least_2_pipelines")
temp3 <- validation7[, sum(pipeline.called >= 3), by="project"]
names(temp3) <- c("project", "at_least_3_pipelines")
temp4 <- validation7[, sum(pipeline.called == 4), by="project"]
names(temp4) <- c("project", "all_4_pipelines")
validation10 <- Reduce(function(...) merge(..., by="project"), list(temp0, temp1, temp2, temp3, temp4))

validation10$at_least_1_pipeline <-  validation10$at_least_1_pipeline / validation10$total
validation10$at_least_2_pipelines <-  validation10$at_least_2_pipelines / validation10$total
validation10$at_least_3_pipelines <-  validation10$at_least_3_pipelines / validation10$total
validation10$all_4_pipelines <-  validation10$all_4_pipelines / validation10$total
validation10$Project <- factor(validation10$Project, levels <- c("UCS", "BRCA", "READ", "OV", "COAD", "SARC", "CESC", "BLCA", "KIRC", "THYM", "PAAD", "LAML"))
validation10 <- validation10[, -"total", with=F]
validation11 <- melt(validation10, value.name ="Proportion_Recalled", variable.name <- "Category")
colnames(validation11) <- c("Project", "Category", "Proportion_Recalled")

# save(validation9, file="validation9.0614.rda")
png("called.validation.by.num.caller.png", width=1024, height=800)
theme_set(theme_gray(base_size <- 18))
p2 <- ggplot(validation9, aes(x=Project, y=Proportion_Recalled, fill=Category)) +
  geom_boxplot(alpha=0.8) +
  geom_line(data <- validation11, aes(x=Project, y=Proportion_Recalled, group <- Category, colour <- Category), size=3, alpha=0.5) +
  theme(legend.position <- c(0.09, 0.12))
p2
dev.off()



