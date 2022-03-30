#!/usr/bin/env Rscript

library(ggplot2)
library(gggenes)
library(data.table)
library(dplyr)

# the main dataframe must have the following columns:
# numeric:
# 	start
# 	end
# 	substart
# 	subend
# 	pos
# 	feature
# 	value
# 	variable
# 	orientation
# 	gene
# 	molecule
# 	mutation_type

# add_extra_col = function(frame, colname) {
# 	if (!colname %in% colnames(frame)) {
# 		frame[,"colname"] = NA
# 	}
# }
# 
# add_extra_cols = function(frame) {
# 	toadd = c("start", "end", "substart", "subend", "pos", "feature", "value", "variable", "orientation"

splitattrs = function(attributes) {
	return sapply(str_split(attributes, ";"), function(x){trimws(x, which = "both")})
}

parseattr = function(attribute) {
	return sapply(str_split(attribute, "="), function(x){trimws(x, which = "both")})
}

attrnames = function(split_attributes) {
	return sapply(attributes, function(x) {parseattr[1]})
}

attrvals = function(split_attributes) {
	return sapply(attributes, function(x) {parseattr[2]})
}

attrs = function(split_attributes) {
	anames = attrnames(split_attributes)
	avals = attrvals(split_attributes)
	names(avals) = anames
	return avals
}

attribute = function(attributes, name) {
	return attrs(splitattrs(attributes))[name]
}

get_genes = function(gff) {
	genes = gff[gff["feature"] == "gene",]
	genes$gene = attribute(genes$attributes, "ID")
	return genes
}

get_transcripts = function(gff) {
	transcripts = gff[gff["feature"] == "transcript",]
	transcripts$transcript = attribute(transcripts$attributes, "Parent")
	return transcripts
}

read_gff = function(path) {
	raw = as.data.frame(fread(path))
	colnames(raw) = c("chrom", "source", "feature", "start", "end", "score", "orientation", "frame", "attributes", "comments")
	genes = get_genes(raw)
	transcripts = get_transcripts(raw)
	data = bind_rows(genes, transcripts)
	data = data[,c("chrom", "feature", "start", "end", "orientation", "gene")]
	return data
}

read_bedgraph_noheader = function(path) {
	raw = as.data.frame(fread(path))
	colnames(raw) = c("chrom", "start", "end", "val", "var")
	starts = raw[, c("chrom", "start", "val", "var")]
	ends = raw[, c("chrom", "end", "val", "var")]
	data = bind_rows(starts, ends)
	colnames(data) = c("chrom", "pos", "value", "variable")
	data$start = 0
	data$end = 0
	data$molecule = 0
	data$gene = ""
	return data
}

read_snps = function(path) {
	data = as.data.frame(fread(path))
	# colnames of data are c("chrom", "start", "end", "mutation_type")
	data$pos = (data$start + data$end) / 2
	data$molecule = 0
	data$gene = ""
	return data
}

# 		geom_gene_arrow(
# 			data = data[data$type=="genes",],
# 			aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)
# 		) +
ggplot_full = function(genome_wide_data, chrom, start, end) {
	data = genome_wide_data[genome_wide_data$chrom == chrom,]
	return ggplot(data=data) +
		geom_gene_arrow(data = data[data$type=="genes",], aes(fill = "white", label = gene)) +
		geom_gene_subarrow(data = data[data$type=="genes",],
			aes(xmin = start, xmax = end, y = molecule, xsubmin = substart, xsubmax = submax),
			color = "black", fill = "dark grey", alpha = .7) +
		geom_point(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_line(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_vline(data = data[data$graph="var",], aes(pos, color = mutation_type)) +
		facet_grid(graph~., scales="free_y")
}

main = function() {
	args = commandArgs(trailingOnly = TRUE)
	features = read_gff(args[1])
	data = read_bedgraph_noheader(args[2])
	vars = read_snps(args[3])

	chrompos = args[4]
	chrompos = str_split_n(chrompos, ":", 3)
	chrom = chrompos[1]
	start = as.numeric(chrompos[2])
	end = as.numeric(chrompos[3])

	outpath = args[5]

	fulldata = bind_rows(data, vars, features)
	p = ggplot_full(fulldata, chrom, start, end)

	pdf(outpath, width = 4, height = 3)
		print(p)
	dev.off()
}

main()
