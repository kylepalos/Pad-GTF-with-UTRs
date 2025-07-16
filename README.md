# Pad-GTF-with-UTRs
Adding fixed length UTR features to a GTF if they are missing.

I've come across some genome annotations that aren't annotated very comprehensively. One example was a GTF that lacked UTRs for basically all protein-coding genes. It's hard to find a tool that will take RNA-seq data and add UTR features or extend terminal exons appropriately (if you know of a straight-forward solution, please let me know.) Therefore, I am using a somewhat crude approach to add on UTR features that should appropriately handle gene orientation (strand) and chromsome boundaries (if you provide a [chrom.sizes file]([url](https://www.biostars.org/p/173963/#9600912)). **Note, this won't discriminate non-coding RNA entries vs protein-coding entries so you may want to split the annotation ahead of time and merge after padding.**

**Example usage:**

```
python gtf_add_utrs.py input.gtf output.gtf -5 100 -3 150 --chrom-sizes chrom.sizes

```
