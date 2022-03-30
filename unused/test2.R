#!/usr/bin/env Rscript

library(ggplot2)
library(gggenes)

print(1)
data = data.frame(y = 1:100, x = as.integer(rnorm(100, 50, 10)), graph = "data", molecule = 0,
	gene = "none", start = 0, end = 0, strand = "forward", orientation = 0)
print(2)

genes = data.frame(
	molecule = c(1, 2),
	gene = c("genA", "genB"),
	start = c(25, 50),
	end = c(75, 101),
	strand = c("forward", "reverse"),
	orientation = c(-1, 1),
	graph = "genes",
	y = 0,
	x=0
)

print(data)
print(genes)
# full = rbind(data, genes)
full = rbind(genes, data)

p = ggplot(data=full) +
	geom_gene_arrow(data = full[full$graph=="genes",], aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
	geom_point(     data = full[full$graph=="data", ], aes(x, y)) +
	facet_grid(graph~., scales = "free_y")

pdf("temp.pdf", width = 4, height = 3)
	print(p)
dev.off()

