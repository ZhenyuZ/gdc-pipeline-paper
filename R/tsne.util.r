############################
# utility functions to run tSNE
##############################

# annotate tsne output with metadata for plot
# tsne is the tsne output, with aliquot barcode as rownames
# metadata is a dataframe with at least one columne "aliquot", and other metadata columns
annotate.tsne = function(tsne, metadata) {
	tsne$aliquot = rownames(tsne)
	data = merge(tsne, metadata, by="aliquot")
	if (nrow(data) != nrow(tsne)) {
		cat("\nWarning!", nrow(tsne)-nrow(data), "aliquots are removed due to lack of metadata")
	}
	return(data)
}


# perform t-SNE algorithm with the first "use" number of columns from the input data
# output is either Y matrix or cost value
# preprocessing with PCA is required to normalize data and reduce dimension for t-SNE input
tsne.single = function(data, use = 200, dims = 2, seed = 1, perplexity = 30, return_cost=F){
	library("Rtsne")
	# Filter by number of columns
	use = min(use, ncol(data))
	data = as.matrix(data)[, 1:use]
	# run
	set.seed(seed)
	tsne = Rtsne(data, dims = dims, pca = F, perplexity = perplexity)
	# select to return Y matrix or cost
	if(return_cost) {
		out = sum(tsne$costs)
	} else {
		out = data.frame(tsne$Y)
		rownames(out) = rownames(data)
		colnames(out) = paste0("feature_", 1:dims)
	}
	return(out)
}

# perform t-SNE algorithm with the first "use" number of columns from the input data
# this function try to pick the best seed that minimizes cost by looping
# output is the best seed among 1:permute
# preprocessing with PCA is required to normalize data and reduce dimension for t-SNE input
tsne.find.seed.by.loop = function(data, use = 200, dims = 2, permute = 1, perplexity = 30){
	# repeated run with different seed to pick the best one
	cat(permute, "rounds of permutation selected:\n\n")
	costs = array(0, permute)
	for(i in 1: permute) {
		cat("  round", i, "out of", permute, "\n")
		costs[i] = tsne.single(data=data, use=use, dims=dims, seed=i, perplexity=perplexity, return_cost=T)
	}
	return(which.min(costs))
}


# perform t-SNE algorithm with the first "use" number of columns from the input data
# this function try to pick the best seed that minimizes cost by parallel
# output is the best seed among 1:permute
# preprocessing with PCA is required to normalize data and reduce dimension for t-SNE input
tsne.find.seed.by.parallel = function(data, use = 200, dims = 2, permute = 1, perplexity = 30){
	library("parallel")
	# Calculate the number of cores
	num_cores <- detectCores() - 1
	# Initiate cluster
	cl <- makeCluster(num_cores, type="FORK")
	clusterExport(cl, c("dims", "perplexity", "use"), envir = .GlobalEnv)
	costs = parLapply(cl, 1:permute, function(x) tsne.single(data=data, use=use, dims=dims, seed=x, perplexity=perplexity, return_cost=T))
	stopCluster(cl)
	return(which.min(unlist(costs)))
}


# perform t-SNE algorithm with the first "use" number of columns from the input data
# this function will try "permute" number of seeds and pick the best one that minimizes cost
# output is the Y matrix
# preprocessing with PCA is required to normalize data and reduce dimension for t-SNE input
tsne = function(data, use = 200, dims = 2, permute = 1, seed = 1, perplexity = 30, parallel = T){
	# Filter by number of columns
	use = min(use, ncol(data))
	data = as.matrix(data)[, 1:use]
	cat("\nUse the first", use, "features to compute t-SNE\n\n")
	cat("t-SNE feature reduction to", dims, "dimensional space\n\n")

	if(permute > 1 & parallel) {
		# repeated run with different seed to pick the best one
		library("parallel")
		# Calculate the number of cores
		num_cores <- detectCores() - 1
		# Initiate cluster
		cl <- makeCluster(num_cores, type="FORK")
		clusterExport(cl, c("dims", "perplexity", "use"), envir = environment())
		costs = parLapply(cl, 1:permute, function(x) tsne.single(data=data, use=use, dims=dims, seed=x, perplexity=perplexity, return_cost=T))
		stopCluster(cl)
		best.seed = which.min(unlist(costs))
	} else if(permute > 1 ) {
		costs = array(0, permute)
		for(i in 1: permute) {
			cat("  round", i, "out of", permute, "\n")
			costs[i] = tsne.single(data=data, use=use, dims=dims, seed=i, perplexity=perplexity, return_cost=T)
		}
		best.seed = which.min(costs)
	} else {
		best.seed = seed
	}
	cat("Best seed selected among", permute, "seeds tested:", best.seed, "\n\n")
	out = tsne.single(data=data, use=use, dims=dims, seed=best.seed, perplexity=perplexity, return_cost=F)
	return(out)
}


# Caculate PCA
# Input: data matrix, feature in rows, and sample in columns
# Output: pca object and variation explained
pca = function(data){
	pca = prcomp(t(data), center = TRUE, scale. = TRUE) 
	pc = pca$x
	sdev = pca$sdev
	total_variation = sum(sdev^2)
	varexp = sdev^2 / total_variation
	names(varexp) = colnames(pc)
	pca$varexp = varexp
	return(pca)
}


# function to select aliquot barcode with the provided sample_type, 
# and only pick one per sample
SelectAliquot = function(my.barcode, my.type = "01"){
	# select sample type code
	my.barcode = my.barcode[substr(my.barcode, 14,15) %in% my.type]
	# sort and remove duplicates
	my.barcode = sort(my.barcode)
	my.barcode = my.barcode[!duplicated(substr(my.barcode, 1, 15), fromLast=T)]
	# return value
	return(my.barcode)
}


# Filter for Completeness
# my.data as a data matrix of numeric values, rows as features and columns as samples
# This function filters features by max missing data rate, max proportaion of zero values, and 0 variance
# 				filters samples by max missing data rate
# benchmark information: about half an hour on a 31 Gb, 8551 column, 485577 row matrix 
CompletenessFilter = function(my.data, feature.missing = 0.05, sample.mising = 0.05, zero.ratio = 1, dup.removal = T) {
	# print out summary
	cat("\nInput data size:", format(object.size(my.data), units = "auto"), "\n")
	cat("Number of samples:", ncol(my.data), "\n")
	cat("Number of features:", nrow(my.data), "\n\n")

	# filter features with NA in more than feature.missing samples
	if (feature.missing < 1 & feature.missing >= 0) {
		cat("Inspecting missing value per feature\n")
		row.filter = which(apply(my.data, 1, function(x) sum(is.na(x))) / ncol(my.data) > feature.missing)
		cat("  ", length(row.filter), "Features contain more than", round(feature.missing * 1000)/10, "% of missing data, and removed\n\n")
		if(length(row.filter) > 0) {
			my.data = my.data[-row.filter, ]
		}
	}

	# filter samples with NA in more than sample.missing features
	if (sample.mising < 1 & sample.mising >= 0) {
		cat("Inspecting missing value per sample\n")
		col.filter = which(apply(my.data, 2, function(x) sum(is.na(x))) / nrow(my.data) > sample.mising)
		cat("  ", length(col.filter), "Samples contain more than", round(sample.mising * 1000)/10, "% of missing data\n")
		if(length(col.filter) > 0) {
			cat("  The following samples are removed\n")
			cat(colnames(my.data)[col.filter], "\n")
			my.data = my.data[, -col.filter]
		}
		cat("\n")
	}

	# remove features that have identical values
	if(dup.removal) {
		cat("Inspecting features with the same values in all samples\n")
		row.filter = which(apply(my.data, 1, function(x) diff(range(x, na.rm=T)) < .Machine$double.eps ^ 0.5))
		cat("  ", length(row.filter), "Features contain only unique values, and removed\n\n")
		if(length(row.filter) > 0) {
			my.data = my.data[-row.filter, ]
		}
	}

	# remove features with too many zero values
	if(zero.ratio < 1 & zero.ratio >= 0) {
		cat("Inspecting features with too many zero values\n")
		row.filter = which(apply(my.data, 1, function(x) sum(x==0) / sum(!is.na(x))) > zero.ratio)
		cat("  ", length(row.filter), "Features contain more than", round(zero.ratio * 1000)/10, "% of zeros, and removed\n\n")
		if(length(row.filter) > 0) {
			my.data = my.data[-row.filter, ]
		}
	}

	# print out summary after filter
	cat("Output data size:", format(object.size(my.data), units = "auto"), "\n")
	cat("Number of samples:", ncol(my.data), "\n")
	cat("Number of features:", nrow(my.data), "\n\n")

	return(my.data)
}



# Missing data imputation by mean
# Need CompletenessFilter to run before ImputeNA
ImputeNA = function(my.data, large.matrix=F){
	# impute NA with row means
	rowmean = apply(my.data, 1, mean, na.rm=T)
	if(large.matrix) {
		# this is a special procedure created to get NA index if 
		# the matrix size is too larget: row x col > 2.147e9 
		row.index = apply(my.data, 2, function(x) which(is.na(x)))
		na.per.col = unlist(lapply(row.index, length))
		na.per.col = cbind(na.per.col, 1:length(na.per.col))
		col = unlist(apply(na.per.col, 1, function(x) rep(x[2], x[1])))
		row = unlist(row.index)
		na.index = as.matrix(cbind(row, col))
	} else {
		na.index = which(is.na(my.data), arr.ind=T)
	}
	my.data[na.index] = rowmean[na.index[,1]]
	return(my.data)
}
