#!/usr/bin/perl

##################################################################################################
#~ bedFromGff                                                                                   ~#
#~ Author: Adrien Vidal                                                                         ~#
##################################################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
    print(
        "bedFromGff\n
        Creates a bed table of the regions from a gff file.
        Usage: bedFromGff [options] -i <in.gff>\n
        Options:
          -i FILE\tPath of the input file (expects gff format.).
          -o FILE\tPath of the output file.
          -ida STR\tName of the attribute containing the identifier of a region. \"ID\" by default.
          -feat STR\tFeature type to select using the third column of the file. This argument may be repeated multipe times to select more than one feature.
          -nona\t\tRemove any region which ID could not be found among the attributes. Otherwise, such regions will have \"NA\" as their ID.
          -forceUnq\tForces all IDs to be unique by appending anumber to every occurence of a particular ID after the first. Ignores IDs equal to \"NA\".
          -h\t\tPrint this help.\n"
    );
    exit 0;
}

my ($inFile) = "";
my ($bed) = "";
my ($idAttr) = "ID";
my (@features) = ();
my ($forbidNa) = 0;
my ($forceUnique) = 0;
my ($help) = 0;

GetOptions( 'i=s'       =>  \$inFile,
            'o=s'       =>  \$bed,
            'ida=s'     =>  \$idAttr,   # String
            'feat=s'    =>  \@features, # String array
            'nona'      =>  \$forbidNa,
            'forceUnq'  =>  \$forceUnique,
            'h'         =>  \$help
) or usage();

usage() if $help;
usage() if ! scalar $inFile;

open (my $fd1, '<', "$inFile") or die("open: $!");
if( scalar $bed ){
    # Open and connection to output file to test that it's possible to write there:
    open(my $test, '>', "$bed") or die("open: $!");
    close($test);
    # Redirect STDOUT to output file:
    open(STDOUT, '>', $bed);
}

my %seen; # Hashmap to keep count of ids:
while (<$fd1>){
    next if($_ =~ m/^[#].*/); # Comment line, skip.
    my @fields = split(/\t/, $_);

    # Filter by specified features:
    if( scalar @features ){ if(!grep( /^$fields[2]$/, @features )){
        next;
    } }

    # Parse attributes field to find id attribute:
    my $id = "NA";
    if ($fields[8] =~ m/$idAttr[\s=]"?([\S]*?)"?[;,\n]/){
        $id = $1;
    }

    if( $id eq "NA" ){
        next if $forbidNa; # NA filtering (if enabled)
    } else {
        if( $forceUnique ){
            # Append number to duplicate ids:
            if( $seen{$id}++ ){
                $id = "$id\_$seen{$id}";
            }
        }
    }
    
    #           Sequence    Start       End         ID   Score       Strand
    my $line = "$fields[0]\t$fields[3]\t$fields[4]\t$id\t$fields[5]\t$fields[6]";
        
    print("$line\n")
}

close($fd1)