#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);
use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


require "/home/sanchez/sda1/perl/MY/my.pl";
$input=$ARGV[0];
$output=$ARGV[1];


$n=0;
open(FH,"<$input");
open(OUT,">$output");
while ($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
#    print "$#data\n";
	if ($#data==2) {
	    $flux=$data[2];
	    $wave=$data[1];
	} else {
	    $flux=$data[1];
	$wave=$data[0];
	}
	$wave_out=$wave*10;
	$n++;	    
	print OUT "$n $wave_out $flux\n";
    }
}
close(FH);
close(OUT);
