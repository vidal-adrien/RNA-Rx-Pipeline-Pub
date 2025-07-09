#!/usr/bin/env perl

##################################################################################################
#~ tx2geneFromGff                                                                               ~#
#~ Author: Adrien Vidal                                                                         ~#
##################################################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
    print(
        "tx2GeneFromGff\n
        Creates a transcript to gene (tx2gene) two column table from a gff file.
        Usage: tx2GeneFromGff [options] -i <in.gff>\n
        Options:
            -i FILE\tPath of the input file (expects gff format.).
            -o FILE\tPath of the output file.
            -ida STR\tName of the attribute containing the identifier of the transcript. \"transcript_id\" by default.
            -gida STR\tName of the attribute containing the identifier of the transcript's gene. \"gene_id\" by default.
            -feat STR\tFeature type to select using the third column of the file. Use to select specific types of transcripts (e.g. mRNA). This argument may be repeated multipe times to select more than one feature. 
            -h\t\tPrint this help.\n"
    );
    exit 0;
}

my ($inFile) = "";
my ($tx2gene) = "";
my ($idAttr) = "transcript_id";
my ($geneIdAttr) = "gene_id";
my (@features) = ();
my ($help) = 0;

GetOptions( 'i=s'       =>  \$inFile,
            'o=s'       =>  \$tx2gene,
            'ida=s'     =>  \$idAttr,       # String
            'gida=s'    =>  \$geneIdAttr,   # String
            'feat=s'    =>  \@features,     # String array
            'h'         =>  \$help
) or usage();

usage() if $help;
usage() if ! scalar $inFile;

open (my $fd1, '<', "$inFile") or die("open: $!");
if( scalar $tx2gene ){
    # Open and connection to output file to test that it's possible to write there:
    open(my $test, '>', "$tx2gene") or die("open: $!");
    close($test);
    # Redirect STDOUT to output file:
    open(STDOUT, '>', $tx2gene);
}

while (<$fd1>){
    next if($_ eq "\n"); # Empty line, skip.
    next if($_ =~ m/^[#].*/); # Comment line, skip.    

    my @fields = split(/\t/, $_);

    # Filter by specified features:
    if( scalar @features ){ if(!grep( /^$fields[2]$/, @features )){
        next;
    } }

    my $id = "NA";
    if ($fields[8] =~ m/$idAttr[\s=]"?([\S]*?)"?[;,\n]/){
        $id = $1;
    }

    my $gid = "NA";
    if ($fields[8] =~ m/$geneIdAttr[\s=]"?([\S]*?)"?[;,\n]/){
        $gid = $1;
    }

    if( $id eq "NA" || $gid eq "NA" ){
        next
    }

    my $line = "$id\t$gid";

    print("$line\n")
}