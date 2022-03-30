#!/usr/bin/env Rscript

library(ggplot2)
library(gggenes)

print(1)
data = data.frame(y = 1:100, start = as.integer(rnorm(100, 50, 10)), graph = "data",
	gene = "none", end = 0, strand = "forward", orientation = 0)
# molecule = 0, 
print(2)

#	molecule = "Potato",
#	molecule = 1,
genes = data.frame(
	gene = c("genA", "genB"),
	start = c(25, 50),
	end = c(75, 101),
	strand = c("forward", "reverse"),
	orientation = c(-1, 1),
	graph = "genes",
	y = 0
)

print(3)
full = rbind(data, genes)
print(full)
print(str(full))

print(4)
p = ggplot(data=full) +
	geom_gene_arrow(data = full[full$graph=="genes",], aes(xmin = start, xmax = end, y = y, fill = gene)) +
	geom_point(     data = full[full$graph=="data", ], aes(start, y)) +
	facet_grid(graph~., scales = "free_y")

print("str(p):")
print(str(p))


print(5)
pdf("temp.pdf", width = 4, height = 3)
print(6)
	print(p)
print(7)
dev.off()
print(8)

