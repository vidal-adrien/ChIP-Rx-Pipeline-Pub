#!/usr/bin/perl

##################################################################################################
#~ bedFromFasta                                                                                 ~#
#~ Author: Adrien Vidal                                                                         ~#
##################################################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
    print(
        "bedfromFasta\n
        Creates a bed table of the full length of the sequences from a fasta file.
        Usage: bedFromFasta -i <in.fasta>\n
        Options:
          -i FILE\tPath of the input file (expects fasta format).
          -o FILE\tPath of the output file.
          -h\t\tPrint this help.\n"
    );
    exit 0;
}

my ($inFile) = "";
my ($bed) = "";
my ($help) = 0;

GetOptions( 'i=s'       =>  \$inFile,
            'o=s'       =>  \$bed,
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

my ($chr);
my ($count) = 0;
my ($firstChr) = 1;
while (<$fd1>){
    next if($_ =~ m/^[;].*/); # Comment line, skip.
    if ($_ =~ m/^>(.*)/) {  # Seq Id line, record the bp count of the previous sequence...
        if( ! $firstChr ){  # ...unless it's the first sequence of the file.
            print ("$chr\t0\t$count\n");
        } else {
            $firstChr = 0;
        }
        $chr = $1;  # New ID.
        $count = 0; # Reset count
    }
    else {
        my $line = $_ =~ s/\n//r;
        $count += length($line);
    }
}
print ("$chr\t0\t$count\n");

close($fd1);