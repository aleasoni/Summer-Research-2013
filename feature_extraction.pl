#!usr/bin/perl
#Alessandro Luca Asoni
#Machine Learning 600.475
#Final Project
#feature_creator.pl

use warnings;
use strict;

my $read   = $ARGV[0]; #enhancer sequences file
my $null   = $ARGV[1]; 
my $kmers  = $ARGV[2]; #reduced k-mer sequences
my $write  = $ARGV[3]; #svm_light format positive examples
my $cluster_size = $ARGV[4]; #base pair clustering length
my $clustering = $ARGV[5]; # if 0 don't cluster

open IN,  "< $read"  or die $!;
open KIN, "< $kmers" or die $!;
open OUT, "> $write" or die $!;

my %inv;

$inv{"A"}="T";
$inv{"C"}="G";
$inv{"G"}="C";
$inv{"T"}="A";

# subroutine to invert DNA sequence (e.g. ATTA is equivalent to TAAT)
sub invert {
    my ($a) = @_;
    my $seq = "";
    foreach my $i (1..length($a)) {
	$seq = $seq.$inv{substr($a,length($a)-$i,1)};
    }
    $seq;
}

my %features;
my $feature_id = 1;

while( <KIN> ) {
    chomp( $_ );
    $features{$_} = $feature_id;
    $features{&invert($_)} = $feature_id;
    $feature_id++;
}

my $line;
my $count = 0;
while($line = <IN>) {
    if( $count % 49 == 0 ) {
	my $progress = ($count/4906) * 100 * 2;
	printf STDERR "%.0f%% complete\n", $progress;
    }
    $count++;
    
    $line = <IN>;
    chomp( $line );
    print OUT $null . " ";
    my %kmer_frequencies;
    foreach my $i (0..(length($line) - 5 ))  {
	my $kmer = substr($line, $i, 5);
	if( not defined ($kmer_frequencies{$kmer}) ) { $kmer_frequencies{$kmer} = 1 }
	else {$kmer_frequencies{$kmer} = $kmer_frequencies{$kmer} + 1};
	$i += 1;
    }
    seek( KIN, 0, 0 );
    while( <KIN> ) {
	chomp( $_ );
	if( defined($kmer_frequencies{$_}) ) {
	    print OUT $features{$_} . ":" . $kmer_frequencies{$_} . " ";
	}
    }

    if( $clustering == 1 ) {

	my %cluster_hash;
	my $length;
	if( length( $line ) < $cluster_size ) {
	    $cluster_size = length( $line );
	}

	foreach my $j (0..(length($line) - $cluster_size )) {
	    my $chunk = substr($line, $j, $cluster_size); # cluster_sized chunck
	    my $reference = substr( $chunk, 0, 5);
	    my $ref_id = $features{$reference};

	    foreach my $k (5..($cluster_size - 5) ) {
		my $kmer = substr( $chunk, $k, 5);
		my $kmer_id = $features{$kmer};
		my $feature_number; 

		if( $ref_id > $kmer_id ) {
		    $feature_number = ( $kmer_id * 512 + $ref_id ) - ( $kmer_id )*($kmer_id - 1)/2;
		} else {
		    $feature_number = ( $ref_id * 512 + $kmer_id ) - ( $ref_id  )*($ref_id -  1)/2;
		}
			
		if( not defined( $cluster_hash{$feature_number}) ) {
		    $cluster_hash{$feature_number} = 1 
		} else {
		    $cluster_hash{$feature_number} = $cluster_hash{$feature_number} + 1;
		}
		    
	    }
	}
	my @keys;

	foreach my $key ( keys %cluster_hash )
	{
	    push( @keys, $key );
	}

	@keys = sort { $a <=> $b } (@keys);

	foreach( @keys ) {
	    print OUT $_ . ":" . $cluster_hash{$_} . " ";
	}
    }


    print OUT "\n";
}

close IN;
close KIN;
close OUT;
