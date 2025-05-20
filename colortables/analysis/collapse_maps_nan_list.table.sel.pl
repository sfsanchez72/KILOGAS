#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
use PDL::Graphics::LUT;
use Carp;
use PDL::Graphics::PGPLOT;

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";

if ($#ARGV<3) {
    print "USE: collapse_maps.pl SUBFIX TABLE PREFIX SELECTIONS(COLUMN,MIN,MAX)\n";
    exit;    
}

$subfix=$ARGV[0];
$table=$ARGV[1];
$prefix=$ARGV[2];
@columns=split(",",$ARGV[3]);
$ncol=int(($#columns+1)/3);
print "# (1) name\n";
print "# (2) survey\n";
print "# (3) Mass\n";
print "# (4) SFR_Ha\n";
print "# (5) EW_Ha_Re\n";
print "# (6) MORPH\n";
print "# (7) ION_TYPE\n";
print "# N.COL = $ncol\n";


open(FH,"<$table");
$n=0;
$n_DO=0;
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	my @data=split(",",$line);
	$name_now=$data[0];
	$survey{$name_now}=$data[1];	
	$DO{$name_now}=1;
	for ($i=0;$i<$ncol;$i++) {
	    my $nc=$columns[$i*3];
	    my $min=$columns[$i*3+1];
	    my $max=$columns[$i*3+2];
	    if (($data[$nc]<=$min)||($data[$nc]>$max)) {
		$DO{$name_now}=$DO{$name_now}*0;
	    }
	}
	$n++;
	if ($DO{$name_now}==1) {
	    $name_DO[$n_DO]=$name_now;
	    $n_DO++;
	}
    }
}
close(FH);

my @file;
$output=$prefix.$subfix;
$e_output=$prefix.".STDDEV".$subfix;
$nz=0;
for (my $i=0;$i<$n_DO;$i++) {
    $name_now=$name_DO[$i];    
    my $dir_pe="";
    if ($survey{$name_now} eq "CALIFA") {
	$dir_pe="/disk-a/sanchez/ppak/legacy/DATA/SSP/proc_elines_v2.2_NEW/";
    }
    if ($survey{$name_now} eq "SAMI") {
	$dir_pe="/disk-c/sanchez/SAMI/DR_v10/data/proc_elines/";
    }
    if ($survey{$name_now} eq "MUSE") {
	$dir_pe="/mnt/synology/bacterio/disk-c/sanchez/MUSE/AMUSING/analysis/proc_elines/";
    }
    if ($survey{$name_now} eq "MaNGA") {
	($manga,$chart,$bundle)=split(/\-/,$name_now);
	$dir_pe="/disk-a/manga/data/v2_4_3/".$chart."/".$bundle."/";
    }	
    my $file_now=$dir_pe."/".$name_now.$subfix;
    if (-e $file_now) {
	$file[$nz]=$file_now;
	$nz++;
    }
}
print "#$nz files in $table with @columns\n";
#exit;



#open(DIR,"ls *$subfix |");
#while($dir=<DIR>) {
#    chop($dir);
#    if ($dir !~ "all") {
#	$file[$nz]=$dir;
#	$nz++;
#    }
#}
#close(DIR);

print "$nz files\n";
for ($k=0;$k<$nz;$k++) {
    my $map=rfits($file[$k]);
    print "$k/$nz $file[$k] ";
    if ($k==0) {
	$h=$map->gethdr();
	($nx,$ny)=$map->dims();
	$pdl_cube=zeroes($nx,$ny,$nz);
	$pdl_cube_w=zeroes($nx,$ny);
    }
    #	$map->inplace->clip(0,1.01);
    $map->inplace->setnantobad;
    $ones=$map/$map;
    $ones->inplace->setnantobad;
    $ones->inplace->setbadtoval(0);
#    $map=$map*$ones;
    my $t = $pdl_cube_w->slice(":,:");
    $t .= $t+$ones;
    my $t = $pdl_cube->slice(":,:,($k)");
    $t .= $map;
    print "\n";
    
}
print "\n";
$old=0;
if ($old==1) {
    $pdl_cube=$pdl_cube/$pdl_cube_w;
    my $c=sumover($pdl_cube->xchg(0,2));
    my $d=$c->xchg(0,1);
    $d->sethdr($h);
    $d->wfits($output);
}

$pdl_cube->inplace->setvaltobad(0);
$pdl_cube->inplace->setnantobad;
$pdl_cube->inplace->clip(-1e12,1e12);
$pdl_cube->inplace->setvaltobad(-1e12);
$pdl_cube->inplace->setvaltobad(1e12);
($mean,$prms,$median,$min,$max,$adev,$rms) = statsover($pdl_cube->xchg(0,2));

my $d=$mean->xchg(0,1);
$d->sethdr($h);
$d->wfits($output);

my $d=$prms->xchg(0,1);
$d->sethdr($h);
$d->wfits($e_output);

#$pdl_cube_w->wfits("file.fits.gz");
exit;
