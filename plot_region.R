#/usr/bin/env Rscript

library(ggplot2)
library(gggenes)
library(data.table)
library(dplyr)

splitattrs = function(attributes) {
	return(sapply(
		strsplit(
			trimws(attributes, which="both", whitespace=";"),
			';'
		),
		function(x){trimws(x, which = "both")}
	))
}

parseattr = function(attribute) {
	return(sapply(strsplit(attribute, '='), function(x){trimws(x, which = "both")}))
}

attrnames = function(split_attributes) {
	return(sapply(split_attributes, function(x) {parseattr(x)[1]}))
}

attrvals = function(split_attributes) {
	return(sapply(split_attributes, function(x) {parseattr(x)[2]}))
}

attrs = function(split_attributes) {
	anames = attrnames(split_attributes)
	avals = attrvals(split_attributes)
	names(avals) = anames
	avals = unlist(avals)
	return(avals)
}

attribute = function(attributes, name) {
	return(attrs(splitattrs(attributes))[name])
}

get_genes = function(gff) {
	genes = gff[gff["feature"] == "gene",]
	genes$gene = attribute(genes$attributes, "ID")
	return(genes)
}

get_transcript_gene = function(gene_name, genes) {
	return(genes[genes["gene"] == gene_name,][1,])
}

get_transcripts = function(gff, genes) {
	transcripts = gff[gff["feature"] == "mRNA",]
	transcripts$gene = attribute(transcripts$attributes, "Parent")
	transcripts$substart = transcripts$start
	transcripts$subend = transcripts$end
	transcripts$start = sapply(transcripts$gene, function(x){get_transcript_gene(x, genes)$start})
	transcripts$end = sapply(transcripts$gene, function(x){get_transcript_gene(x, genes)$end})
	return(transcripts)
}

orientone = function(strandstr) {
	if (strandstr == "-") {
		return(-1)
	}
	return(1)
}

orient = function(strandvec) {
	return(sapply(strandvec, orientone))
}

read_gff = function(path) {
	raw = as.data.frame(fread(path))
	colnames(raw) = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
	raw$orientation = orient(raw$strand)
	genes = get_genes(raw)
	transcripts = get_transcripts(raw, genes)
	data = bind_rows(genes, transcripts)
	data = data[,c("chrom", "feature", "start", "end", "substart", "subend", "orientation", "gene")]
	data$type = "genes"
	data$graph = "genes"
	data$molecule = 1
	return(data)
}

read_bedgraph_noheader = function(path) {
	raw = as.data.frame(fread(path))
	colnames(raw) = c("chrom", "start", "end", "val", "var")
	starts = raw[, c("chrom", "start", "val", "var")]
	ends = raw[, c("chrom", "end", "val", "var")]
	data = bind_rows(starts, ends)
	colnames(data) = c("chrom", "pos", "value", "variable")
	data$type = "data"
	data$graph = "data"
	return(data)
}

read_snps = function(path) {
	data = as.data.frame(fread(path))
	colnames(data) = c("chrom", "start", "end", "mutation_type")
	data$pos = (data$start + data$end) / 2
	data$type = "var"
	data$graph = "genes"
	return(data)
}

ggplot_full = function(genome_wide_data, chrom, start, end) {
	data = genome_wide_data[genome_wide_data$chrom == chrom,]
	p = ggplot(data=data) +
		geom_gene_arrow(data = data[data$type=="genes",], aes(y = molecule, fill = "white", label = gene, xmin = start, xmax = end)) +
		geom_subgene_arrow(data = data[data$type=="genes",],
			aes(y = molecule, xmin = start, xmax = end, xsubmin = substart, xsubmax = subend),
			color = "black", fill = "dark grey", alpha = .7) +
		geom_point(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_line(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_vline(data = data[data$graph=="var",], aes(pos, color = mutation_type)) +
		facet_grid(graph~., scales="free_y") +
		scale_x_continuous(lim = c(start, end))
	return(p)
}

main = function() {
	args = commandArgs(trailingOnly = TRUE)
	features = read_gff(args[1])
	data = read_bedgraph_noheader(args[2])
	vars = read_snps(args[3])

	chrompos = args[4]
	chrompos = strsplit(chrompos, ':')
	chrom = chrompos[[1]][1]
	start = as.numeric(chrompos[[1]][2])
	end = as.numeric(chrompos[[1]][3])

	outpath = args[5]

	fulldata = bind_rows(data, vars, features)
	p = ggplot_full(fulldata, chrom, start, end)

	pdf(outpath, width = 4, height = 3)
		print(p)
	dev.off()
}

main()
