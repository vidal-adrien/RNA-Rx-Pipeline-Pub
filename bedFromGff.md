# bedFromGff.pl

Creates a `.bed` table of the regions from a `.gff` file.

Usage:`bedFromGff [options] -i <in.gff>`

Options:
*  `-i FILE`         Path of the input file (expects gff format.).
*  `-o FILE`         Path of the output file.
*  `-ida STR`        Name of the attribute containing the identifier of a region. "ID" by default.
*  `-feat STR`       Feature type to select using the third column of the file. This argument may be repeated multipe times to select more than one feature.
*  `-nona`           Remove any region which ID could not be found among the attributes. Otherwise, such regions will have "NA" as their ID.
*  `-forceUnique`    Forces all IDs to be unique by appending anumber to every occurence of a particular ID after the first. Ignores IDs equal to "NA".
*  `-h`              Print this help.

https://github.com/vidal-adrien/RNA-Rx-Pipeline-Pub/blob/d5bc9d242354fe4a5ff61abd4fafb116e30a1c07/perl/bedFromGff.pl#L1-L92
