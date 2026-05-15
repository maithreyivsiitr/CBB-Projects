#!/bin/bash

# Create genome size file
cut -f1,2 hg38.fa.fai | sort -k1,1 > hg38.genome

# Generate promoter regions
bedtools slop \
  -i genes_tss.clean.bed \
  -g hg38.genome \
  -l 500 -r 500 -s > promoters.bed

# Extract promoter sequences
bedtools getfasta \
  -fi hg38.fa \
  -bed promoters.bed \
  -s -name \
  -fo promoters.fa

# Find motif-containing promoters
grep -B1 "GCGC" promoters.fa > motif_hits.txt

# Extract gene names
grep ">" motif_hits.txt | \
sed 's/>//' | \
cut -d"|" -f2 | \
sort | uniq > genes_with_motif.txt
