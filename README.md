# PNDS RNA-seq/RNA-Rx pipelines
[![DOI](https://zenodo.org/badge/1016758818.svg)](https://doi.org/10.5281/zenodo.17397631)

This is a public collection of documentation about RNA-seq and RNA-Rx pipelines used by the PNDS (Plant Nuclear Dynamics & Signaling) Team led by Clara Bourbousse & Fredy Barneche.

Scripts and documentation by Adrien Vidal.

## Pipeline guides
Step-by-step guides to realize analyses.

[RNA-Rx pipeline](rnarx.md)

[Selection of stable genes](rnaRxStableGenes.md)

## Custom tools
Scripts used in the above guides and other custom tools.

[bedFromFasta.pl](bedFromFasta.md): Perl script. Creates a `.bed` table of the full length of the sequences from a `.fasta` file.

[bedFromGff.pl](bedFromGff.md): Perl Script. Creates a `.bed` table of the regions from a `.gff` file. With the possibility to specify which tag contains the ID, to filter by feature type and to enforce ID uniqueness.

[tx2GeneFromGff.pl](tx2GeneFromGff.md) Perl Script. Extracts the transcript to gene correspondance from a `.gff` annotation file.
## Genomics resources:

**Reference genomes:**
*  [⇗TAIR10_chr_all.fas.gz](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz): *Arabidopsis thaliana* TAIR10 genome assembly.
*  [⇗Download page](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/) for *Drosophila melanogaster* release 6 genome assembly.

**Annotation:**
*  [Araport11_GFF3.gene.201606.bed](resources/Araport11_GFF3.gene.201606.bed): Araport 11 annotation for genes on *Arabidopsis thaliana* TAIR10 genome assembly as a `.bed` file. Converted from the [june 2016 annotation](https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/archived/Araport11_GFF3_genes_transposons.Jun2016.gff.gz).

**Blacklist:**
*  [TAIR10_blacklist.bed](resources/TAIR10_blacklist.bed): A blacklist of aberrant regions for the *Arabidopsis thaliana* TAIR10 genome.
