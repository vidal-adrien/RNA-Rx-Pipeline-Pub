# PNDS Genomics pipelines

## Pipeline guides
Step-by-step guides to realize analyses.

[RNA-Rx pipeline](rnarx.md)

[RNA-Rx stable genes](rnaRxStableGenes.md)

## Custom tools
Scripts used in the above guides and other custom tools.

[bedFromFasta.pl](bedFromFasta.md): Perl script. Creates a `.bed` table of the full length of the sequences from a `.fasta` file.

[bedFromGff.pl](bedFromGff.md): Perl Script. Creates a `.bed` table of the regions from a `.gff` file. With the possibility to specify which tag contains the ID, to filter by feature type and to enforce ID uniqueness.

[mergeOverlappingRegions.sh](mergeOverlappingRegions.md): Bash script. Uses a comination of bedtools functions to search for overlaps of the genomic regions of the between two files and generate merged regions bed file out of selected regions.

[validPairsForJuicer.sh](validPairsForJuicer.md): Bash script. Converts a validPair format from HiC-pro into a valid format for juicer_tools pre.

[Genomic tools wrappers and HTCondor cluster launcher](wrappers.md): A set of wrapper scripts for many steps of the genomics pipelines and a system for launching scripts and HTCondor submissions through R.

## Genomics resources:

**Reference genomes:**
*  [⇗TAIR10_chr_all.fas.gz](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz): *Arabidopsis thaliana* TAIR10 genome assembly.
*  [⇗Download page](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/) for *Drosophila melanogaster* release 6 genome assembly.

**Annotation:**
*  [Araport11_GFF3.gene.201606.bed](resources/Araport11_GFF3.gene.201606.bed): Araport 11 annotation for genes on *Arabidopsis thaliana* TAIR10 genome assembly as a `.bed` file. Converted from the [june 2016 annotation](https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/archived/Araport11_GFF3_genes_transposons.Jun2016.gff.gz).

**Blacklist:**
*  [TAIR10_blacklist.bed](resources/TAIR10_blacklist.bed): A blacklist of aberrant regions for the *Arabidopsis thaliana* TAIR10 genome.