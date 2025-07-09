# tx2GeneFromGff.pl

Creates a transcript to gene (tx2gene) two column table from a `.gff` file.

Usage:`tx2GeneFromGff [options] -i <in.gff>`

Options:
*  `-i FILE`        Path of the input file (expects gff format.).
*  `-o FILE`        Path of the output file.
* `-ida STR`        Name of the attribute containing the identifier of the transcript. \"transcript_id\" by default.
* `-gida STR`       Name of the attribute containing the identifier of the transcript's gene. \"gene_id\" by default.
* `-feat STR`       Feature type to select using the third column of the file. Use to select specific types of transcripts (e.g. mRNA). This argument may be repeated multipe times to select more than one feature.
* `-h`              Print this help."