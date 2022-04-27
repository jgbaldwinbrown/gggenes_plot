#!/bin/bash
set -e

if [ ! -f setup.done ] ; then
	bash setup.sh
	touch setup.done
fi

Rscript plot_region.R \
	louse_annotation_0.1.1_chr1.gff.gz \
	data.bed \
	snps.bed \
	"chr1:1000000:1200000" \
	potato.pdf \
	8:3
