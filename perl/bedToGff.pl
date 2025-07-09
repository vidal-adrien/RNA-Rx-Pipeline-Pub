#!/usr/bin/perl

##################################################################################################
#~ bedToGff                                                                                     ~#
#~ Author: Adrien Vidal                                                                         ~#
##################################################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
    print(
        "bedToGff\n
        Creates a gff table of the regions from a bed file.
        Usage: bedFromGff [options] -i <in.bed>\n
        Options:
            -i FILE\tPath of the input file (expects bed format.).
            -o FILE\tPath of the output file.
            -feat STRING\tEither a string like 'COL<x>' where x is the bed file's column to read as the feature type variable or any other string to write as the feature for all rows.
            -src STRING\tEither a string like 'COL<x>' where x is the bed file's column to read as the source variable or any other string to write as the feature for all rows.
            -h\t\tPrint this help.\n"
    );
    exit 0;
}

my ($inFile) = "";
my ($gff) = "";
my ($idAttr) = "ID";
my ($feature) = "locus";
my ($source) = "bed_file";
my ($phase) = ".";
my ($forceUnique) = 0;
my ($help) = 0;

GetOptions( 'i=s'       =>  \$inFile,
            'o=s'       =>  \$gff,
            'ida=s'     =>  \$idAttr,
            'feat=s'    =>  \$feature,
            'src=s'     =>  \$source,
            'pha=s'     =>  \$phase,
            'forceUnq'  =>  \$forceUnique,
            'h'         =>  \$help
) or usage();

usage() if $help;
usage() if ! scalar $inFile;

open (my $fd1, '<', "$inFile") or die("open: $!");
if( scalar $gff ){
    # Open and connection to output file to test that it's possible to write there:
    open(my $test, '>', "$gff") or die("open: $!");
    close($test);
    # Redirect STDOUT to output file:
    open(STDOUT, '>', $gff);
}

my %seen; # Hashmap to keep count of ids:
while (<$fd1>){
    next if($_ eq "\n"); # Empty line, skip.
    next if( $_ =~ m/^[#].*/ ); # Comment line, skip.
    my @fields = split(/\t/, $_);

    my $S = "";
    if ($source =~ m/^col([0-9]+)$/){
        $S = "$fields[$1]";
    } else {
        $S = $source;
    }

    my $F = "";
    if ($feature =~ m/^col([0-9]+)$/){
        $F = "$fields[$1]";
    } else {
        $F = $feature;
    }

    my $P = "";
    if ($phase =~ m/^col([0-9]+)$/){
        $P = "$fields[$1]";
    } else {
        $P = $phase;
    }

    my $id = $fields[3];

    if( $forceUnique ){
        # Append number to duplicate ids:
        if( $seen{$id}++ ){
            $id = "$id\_$seen{$id}";
        }
    }

    #                       Source
    #                       |   Feature                                             Phase
    #                       |   |                                                   |
    #           Sequence    V   V   Start       End         Score       Strand      V   Attributes:id
    my $line = "$fields[0]\t$S\t$F\t$fields[1]\t$fields[2]\t$fields[4]\t$fields[5]\t$P\t$idAttr=$id;";
    $line =~ s/\n//g;

    print("$line\n")
}

close($fd1)