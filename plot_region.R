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
	# print("attributes:")
	# print(attributes)
	# print("name:")
	# print(name)
	# sa = splitattrs(attributes)
	# a = attrs(sa)
	# out = a[name]
	# print("sa")
	# print(sa)
	# print("a")
	# print(a)
	# print("out")
	# print(out)
	return(attrs(splitattrs(attributes))[name])
}

oneparent = function(string) {
	split = strsplit(string, ",")
	return(split[[1]][[1]])
}

get_genes = function(gff) {
	genes = gff[gff["feature"] == "gene",]
	genes$gene = sapply(genes$attributes, function(x) {attribute(x, "ID")})
	return(genes)
}

get_transcript_gene = function(gene_name, genes) {
	return(genes[genes["gene"] == gene_name,][1,])
}

get_transcripts = function(gff, genes) {
	transcripts = gff[gff["feature"] == "mRNA",]
	transcripts$gene = sapply(transcripts$attributes, function(x){unname(unlist(oneparent(attribute(x, "Parent"))))})
	transcripts$substart = transcripts$start
	transcripts$subend = transcripts$end
	transcripts$start = sapply(transcripts$gene, function(x){get_transcript_gene(x, genes)$start})
	transcripts$end = sapply(transcripts$gene, function(x){get_transcript_gene(x, genes)$end})
	transcripts$transcript = sapply(transcripts$attributes, function(x) {attribute(x, "ID")})
	# print("transcripts$transcript:")
	# print(transcripts$transcript)
	return(transcripts)
}

get_exon_transcript = function(transcript_name, transcripts) {
	# print("exon transcript name:")
	# print(transcript_name)
	
	out = transcripts[transcripts["transcript"] == transcript_name,][1,]
	# print("exon transcript:")
	# print(out)
	return(out)
}

get_exon_gene = function(transcript_name, transcripts, genes) {
	transcript = get_exon_transcript(transcript_name, transcripts)
	# print("transcript:")
	# print(transcript)
	gene_name = transcript$gene
	# print("gene_name:")
	# print(gene_name)
	gene = genes[genes["gene"] == gene_name,][1,]
	# print("gene:")
	# print(gene)
	# return(genes[genes["gene"] == gene_name,][1,])
	return(gene)
}

get_exons = function(gff, genes, transcripts) {
	exons = gff[gff["feature"] == "exon",]
	exons$transcript = sapply(exons$attributes, function(x){unname(unlist(oneparent(attribute(x, "Parent"))))})
	exons$gene = sapply(exons$transcript, function(x){get_exon_gene(x, transcripts, genes)$gene})
	# print("genes:")
	# print(exons$gene)
	exons$substart = exons$start
	exons$subend = exons$end
	exons$start = sapply(exons$gene, function(x){get_transcript_gene(x, genes)$start})
	exons$end = sapply(exons$gene, function(x){get_transcript_gene(x, genes)$end})
	return(exons)
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
	exons = get_exons(raw, genes, transcripts)
	data = bind_rows(genes, transcripts, exons)
	data = data[!(data$feature %in% c("mRNA")),]
	# data = exons

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
		geom_gene_arrow(data = data[data$type=="genes",], aes(y = molecule, xmin = start, xmax = end), fill = "white") +
		geom_subgene_arrow(data = data[data$type=="genes",],
			aes(y = molecule, xmin = start, xmax = end, xsubmin = substart, xsubmax = subend),
			color = "black", fill = "dark grey", alpha = .7) +
		geom_point(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_line(data = data[data$type=="data",], aes(pos, value, color = variable)) +
		geom_vline(data = data[data$graph=="var",], aes(pos, color = mutation_type)) +
		geom_gene_label(aes(y = molecule, xmin = start, xmax = end, label = gene)) +
		facet_grid(graph~., scales="free_y") +
		scale_y_log10() +
		scale_x_continuous(lim = c(start, end)) +
		theme_classic()
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

	xydims = args[6]
	xydims = strsplit(xydims, ':')
	xdim = as.numeric(xydims[[1]][1])
	ydim = as.numeric(xydims[[1]][2])

	fulldata = bind_rows(data, vars, features)
	write.table(fulldata, "fulldata.txt")
	p = ggplot_full(fulldata, chrom, start, end)

	pdf(outpath, width = xdim, height = ydim)
		print(p)
	dev.off()
}

main()
