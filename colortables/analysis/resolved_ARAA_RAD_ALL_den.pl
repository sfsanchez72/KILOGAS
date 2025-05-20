#!/usr/bin/perl
use PDL::Core;
use PDL::Graphics::Gnuplot;
use PDL::Basic;
use PDL::Primitive;
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
use Statistics::OLS;
use PDL::Transform::Color;

use Statistics::R;
use Graphics::Gnuplot::Palettes qw(palette palette_names);
use PDL::IO::Pic;
use warnings;
no warnings "all";
use PDL::Graphics::PLplot;
use PDL::NiceSlice;

require("./my_ARAA_BPT_maps_den.pl");


print "Reading proc_elines files\n";


$n_CALIFA=0;
$n_MUSE=0;
$n_SAMI=0;
$n_MaNGA=0;
$n=0;
$SAMPLE_NOW="CALIFA";;
@name_CALIFA=read_get_proc($get_proc_CALIFA);
$n_CALIFA=$n;
print "CALIFA, $n_CALIFA\n";
$SAMPLE_NOW="MaNGA";;
@name_MaNGA=read_get_proc($get_proc_MaNGA);
$n_MaNGA=$n-$n_CALIFA;
print "MaNGA, $n_MaNGA\n";
$SAMPLE_NOW="SAMI";;
@name_SAMI=read_get_proc($get_proc_SAMI);
$n_SAMI=$n-$n_CALIFA-$n_MaNGA;
print "SAMI, $n_SAMI\n";
$SAMPLE_NOW="MUSE";
@name_MUSE=read_get_proc($get_proc_MUSE);
$n_MUSE=$n-$n_CALIFA-$n_MaNGA-$n_SAMI;
print "MUSE, $n_MUSE\n"; 

my @name=(@name_CALIFA,@name_MaNGA,@name_SAMI,@name_MUSE);
print "$n total cubes\n";


#
# READ PROC_ELINES
#
print "Reading mag_files files\n";


read_get_mag($get_mag_cubes_CALIFA);
read_get_mag($get_mag_cubes_MaNGA);
read_get_mag($get_mag_cubes_SAMI);
read_get_mag($get_mag_cubes_MUSE);
def_tables_get_mag();


#print "$name[0],$name[1000],$name[6000]";
#print "$name_CALIFA[0],$name_MaNGA[0],@name_CALIFA";

my @x;
my @x_2;
my @x_3;
my @e_x;
my @e_x_2;
my @e_x_3;
my @y;
my @e_y;
my @color;

##################################################################
# End reading tables
##################################################################



##################################################################
# Reading morphological tables!
##################################################################
#
# Morphology
#

#@morph_name=("E0","E1","E2","E3","E4","E5","E6","E7","S0","S0a","Sa","Sab","Sb","Sbc","Sc","Scd","Sd","Sdm","I","BCD");
@morph_name=("E","E","E","E","E","E","E","E","S0","S0a","Sa","Sab","Sb","Sbc","Sc","Scd","Sd","Sdm","I","BCD");
@bar_name=("A","AB","B");



#
# MaNGA
#
$N_PHOT=0;
$phot_file="get_mag_cubes_v2_2_0_cen.csv";
open(FH,"<$phot_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(",",$line);
	$name_now=$data[0];
	$PHOTOMETRY{$name_now}=$line;
	$MORPH{$name_now}="none";
	$morph{$name_now}=-1;
	$bar{$name_now}=-1;
	$N_PHOT++;
    }
}
print "## Cubes with photometry = $N_PHOT in MANGA\n";

#
# MaNGA, morphological classification
#
my $n_morph_MaNGA=0;
open(FH,"<clasificacion_manga_mpl4_mpl5.txt");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$n_morph_MaNGA++;
	my @data=split(" ",$line);
	my $name=$data[0];
	my $PHOT=$PHOTOMETRY{$name};
	my @phot_now=split(",",$PHOT);
	my $elip=int($phot_now[37]*10);
        my $line_new=$line;
	$line_new =~ s/$name//g;
#	$line_new =~ s///g;
	$line_new =~ s/\d//g;
	$line_new =~ s/\.//g;
	my @data_l=split(" ",$line_new);
	my $class1=$data_l[0];
#	print "$name|$class1\n";


#	my $class1=$data[4];
	$MORPH{$name}=$class1;
	$morph{$name}=-1;
#
# Presence of bars
#
#	print "#$name,$bar{$name}\n";

	if ($class1 =~ "S") {
	    if ($class1 =~ 'AB') { 	    
		$bar{$name} = 1; 
		$cut="AB";
		$class1 =~ s/$cut//g;
	    } else {
		if ($class1 =~ 'A') { 
		    $bar{$name} = 0; 
		    $cut="A";
		    $class1 =~ s/$cut//g;		
		} else {
		    if ($class1 =~ 'B') { 
			$bar{$name} = 2; 
			$cut="B";
			$class1 =~ s/$cut//g;		
		    } 
		}
	    }
	}

#
# Morphology
#

	if ($class1 =~ "BCD") {
	    $morph{$name}=19;
	}
#	if ($class1 =~ "edge") {
#	    $morph{$name}=16;
#	}
	if ($class1 =~ "S") {
	    $morph{$name}=13;
#	    if ($class1 =~ '0')  { $morph{$name}=8; }
	    if ($class1 =~ 'O')  { $morph{$name}=8; }
	    if ($class1 =~ 'a')  { $morph{$name}=10; }
	    if ($class1 =~ 'b')  { $morph{$name}=12; }
	    if ($class1 =~ 'c')  { $morph{$name}=14; }
	    if ($class1 =~ 'd')  { $morph{$name}=16; }
	    if ($class1 =~ 'I')  { $morph{$name}=18; }
	    if ($class1 =~ 'm')  { $morph{$name}=18; }
#	    if ($class1 =~ '0a') { $morph{$name}=9; }
	    if ($class1 =~ 'Oa') { $morph{$name}=9; }
	    if ($class1 =~ 'ab') { $morph{$name}=11; }
	    if ($class1 =~ 'bc') { $morph{$name}=13; }
	    if ($class1 =~ 'cd') { $morph{$name}=15; }
	    if ($class1 =~ 'dm') { $morph{$name}=17; }	    
	}
	if ($class1 =~ "I") {
	    $morph{$name}=18;
	}

	if ($class1 =~ "E") {
	    if ($elip>7) {
		$elip=7;
	    }
	    $morph{$name}=7;#$elip;
	}

	if ($class1 =~ "Ecompact") {
	    if ($elip>7) {
		$elip=7;
	    }
	    $morph{$name}=7;#$elip;
	}


#	if ($morph{$name}==0) {
#	    print "#$name,$MORPH{$name},$bar{$name},$morph{$name}\n"; 
#	}
    }
}

print "## Cubes with morphology = $n_morph_MaNGA in MANGA\n";

#######################################################
# CALIFA
######################################################
$N_PHOT=0;
#$phot_file="/home/sanchez/Documents/CALIFA/DR3/get_mag_cubes_v2.2.csv";
$phot_file="/home/sanchez/Documents/CALIFA/DR4/get_mag_cubes_v2.2.csv";
open(FH,"<$phot_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(",",$line);
	$name_now=$data[0];
	$C{$name_now}=$data[43];
	$e_C{$name_now}=$data[44];
	$Mabs_R{$name_now}=$data[29];
	$e_Mabs_R{$name_now}=$data[30];
	$BR{$name_now}=$data[54];
	$e_BR{$name_now}=$data[55];
	$PHOTOMETRY{$name_now}=$line;
	$MORPH{$name_now}="none";
	$morph{$name_now}=-1;
	$bar{$name_now}=-1;
	$N_PHOT++;
    }
}
print "## Cubes with photometry = $N_PHOT CALIFA\n";

$file="/home/sanchez/Documents/CALIFA/DR3/CALIFA_basic_joint.csv";
open(FH,"<$file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
        @data=split(/\,/,$line);
        my $id=int($data[0]*1.0);
        $name_now=$data[1];
        $NAME_ID{$id}=$name_now;
    }
}
close(FH);

#
# Need to be updated!
#
open(FH,"</home/sanchez/Documents/CALIFA/DR3/CALIFA_3_joint_classnum.csv");
$n_morph=0;
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$n_morph++;
	$line =~ s/ //g;
	my @data=split(/\,/,$line);
	my $name=$data[1];
	my $PHOT=$PHOTOMETRY{$name};
	my @phot_now=split(",",$PHOT);
	my $elip=int($phot_now[37]*10);
        my $line_new=$line;
	$line_new =~ s/$name//g;
	$line_new =~ s/\d//g;
	$line_new =~ s/\.//g;
	my @data_l=split(" ",$line_new);
	my $class1=$data_l[0];
	$morph{$name}=$data[5];
	$bar{$name}=$data[6];
	$MORPH{$name}=$morph_name[$data[5]];
    }
}
close(FH);
print "## $n_morph galaxies with morphology in CALIFA\n";


################################################################
# SAMI
################################################################


$N_PHOT_sami=0;
$phot_file="tables/get_mag_cubes_SAMI.csv";
open(FH,"<$phot_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(",",$line);
	$name_now=$data[0];
	$C{$name_now}=$data[43];
	$e_C{$name_now}=$data[44];
	$Mabs_R{$name_now}=$data[29];
	$e_Mabs_R{$name_now}=$data[30];
	$BR{$name_now}=$data[54];
	$e_BR{$name_now}=$data[55];
	$PHOTOMETRY{$name_now}=$line;
	$MORPH{$name_now}="none";
	$morph{$name_now}=-1;
	$bar{$name_now}=-1;
	$N_PHOT_sami++;
    }
}
print "## Cubes with photometry = $N_PHOT_sami SAMI\n";

$n_morph_SAMI=0;
$file="tables/SAMI_Morph_name.csv";
open(FH,"<$file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
        @data=split(/\,/,$line);
        $name_now=$data[0];
        $morph{$name_now}=$data[1];
	if ($data[1]>-1) {
	    $n_morph_SAMI++;
	}
    }
}
close(FH);

print "## $n_morph galaxies with morphology in SAMI\n";


###########################################################################
# MUSE
###########################################################################


$N_PHOT_MUSE=0;
#$phot_file="/home/sanchez/Documents/CALIFA/DR3/get_mag_cubes_v2.2.csv";
$phot_file="tables/get_mag_cubes_MUSE.csv";
open(FH,"<$phot_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(",",$line);
	$name_now=$data[0];
	$C{$name_now}=$data[43];
	$e_C{$name_now}=$data[44];
	$Mabs_R{$name_now}=$data[29];
	$e_Mabs_R{$name_now}=$data[30];
	$BR{$name_now}=$data[54];
	$e_BR{$name_now}=$data[55];
	$PHOTOMETRY{$name_now}=$line;
	$MORPH{$name_now}="none";
	$morph{$name_now}=-1;
	$bar{$name_now}=-1;
	$N_PHOT_MUSE++;
    }
}
print "## Cubes with photometry = $N_PHOT_MUSE MUSE\n";


#
# Need to be updated!
#
open(FH,"<tables/MUSE_morphology.csv");
$n_morph_MUSE=0;
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$line =~ s/ //g;
	my @data=split(/\,/,$line);
	my $name=$data[1];
	my $morph_now=$data[4];
	my $bar_now=$data[5];
	if ($morph_now ne "NO") {
	    $morph{$name}=-1;
	    $MORPH{$name}="NO";
#$morph_name[$data[5]];
	    
	    for (my $i=0;$i<$#morph_name+1;$i++) {
		if ($morph_now eq $morph_name[$i]) {
		    $morph{$name}=$i;
		    $MORPH{$name}=$morph_name[$i];
		}
	    }
	    if ($morph{$name}>-1) {
		if ($morph{$name}>7) {
		    if ($bar_now==0) {
			$bar{$name}="A";
		    } else {
			$bar{$name}="B";
		    }
		} else {
		    $bar{$name}="";
		}

		$n_morph_MUSE++;
#		print "$name,$morph{$name},$MORPH{$name},$bar{$name}\n";
	    }
	}
    }
}
close(FH);
print "## $n_morph_MUSE galaxies with morphology in MUSE\n";
#exit;






$n_morph=$n_morph+$n_morph_MaNGA+$n_morph_SAMI+$n_morph_MUSE;
print "## $n_morph galaxies with morphology in total\n";



#
# We glue all morphologies
#

for ($i=0;$i<$n;$i++) {
    $name_now=$name[$i];
    my @data=split(",",$LINE_PROC_ELINES{$name_now});
    if ($name_now =~ "manga") {
        def_tables_MaNGA();        
    } else {
        def_tables();    
    }
    $A_MORPH[$i]=$morph{$name_now};#$data[$n_Mass];
}
##################################################################
# ENDReading morphological tables!
##################################################################
# Output file
print "N = $n\n";

###################################################################
# Reading radial profiles


#******************************************************************

print "BPT by morphology\n";

#########################################################################
# BPT diagram PLPLOT MORPH
#########################################################################

# reading BPT diagrams!
my @diag_files;
my @diag_files_EW;

$nf=0;
my @diag_files;
my @diag_files_EW;
open(DIR,"ls fitsfiles/*.RAD_ALL.fits.gz |");
while($file=<DIR>) {
    chop($file);
    if ($file !~ "STDDEV") {
	$diag_files[$nf]=$file;
	$file =~ s/MAP/fMAP/g;
	$file =~ s/diag/diag_EW/g;
	$diag_files_EW[$nf]=$file;
	print "$nf,$diag_files[$nf]\n";
	$nf++;
    }
}
close(DIR);
print "$nf DIAG files\n";
exit;


my $pdl_cube;
my $pdl_fcube;
for ($k=0;$k<$nf;$k++) {
	my $map=rfits($diag_files[$k]);
	my $fmap=rfits($diag_files_EW[$k]);
	if ($k==0) {
	    $h=$map->gethdr();
	    ($NX,$NY)=$map->dims();
	    $pdl_cube=zeroes($NX,$NY,$nf);
	    $pdl_fcube=zeroes($NX,$NY,$nf);
	    $crval1=$$h{CRVAL1}; 
	    $crval2=$$h{CRVAL2};
	    $cdelt1=$$h{CDELT1}; 
	    $cdelt2=$$h{CDELT2};
	}
#	print "$NX,$ny,$k\n";
	my $t = $pdl_cube->slice(":,:,($k)");
	$t .= $map;
	my $t = $pdl_fcube->slice(":,:,($k)");
	$t .= $fmap;#->log10;
}
print "DONE\n"; 

my $pdl_fcube_EW=$pdl_fcube;
my $pdl_cube_ORG=$pdl_cube;


exit;

#******************************************************************

print "BPT by morphology\n";

#########################################################################
# BPT diagram PLPLOT MORPH
#########################################################################

# reading BPT diagrams!
my @diag_files;
my @diag_files_EW;

$nf=0;
my @diag_files;
my @diag_files_EW;
open(DIR,"ls fitsfiles/*.MAP_diag.fits.gz |");
while($file=<DIR>) {
    chop($file);
    if ($file !~ "STDDEV") {
	$diag_files[$nf]=$file;
	$file =~ s/MAP/fMAP/g;
	$file =~ s/diag/diag_EW/g;
	$diag_files_EW[$nf]=$file;
	print "$nf,$diag_files[$nf]\n";
	$nf++;
    }
}
close(DIR);
print "$nf DIAG files\n";

my $pdl_cube;
my $pdl_fcube;
for ($k=0;$k<$nf;$k++) {
	my $map=rfits($diag_files[$k]);
	my $fmap=rfits($diag_files_EW[$k]);
	if ($k==0) {
	    $h=$map->gethdr();
	    ($NX,$NY)=$map->dims();
	    $pdl_cube=zeroes($NX,$NY,$nf);
	    $pdl_fcube=zeroes($NX,$NY,$nf);
	    $crval1=$$h{CRVAL1}; 
	    $crval2=$$h{CRVAL2};
	    $cdelt1=$$h{CDELT1}; 
	    $cdelt2=$$h{CDELT2};
	}
#	print "$NX,$ny,$k\n";
	my $t = $pdl_cube->slice(":,:,($k)");
	$t .= $map;
	my $t = $pdl_fcube->slice(":,:,($k)");
	$t .= $fmap;#->log10;
}
print "DONE\n"; 

my $pdl_fcube_EW=$pdl_fcube;
my $pdl_cube_ORG=$pdl_cube;


$zmin=-0.49;
$zmax=2.65;
$x_min1=$crval1;
$x_max1=$crval1+$NX*$cdelt1;
$y_min1=$crval2;
$y_max1=$crval2+$NY*$cdelt2;



$pdl_cut_x=-3+0.02*pdl(0..300);
$pdl_cut_y=-0.7+0.2-3.67*$pdl_cut_x;
$pdl_cut_y2=-1.7+0.5-3.67*$pdl_cut_x;
$pdl_cut_y3=0.61/($pdl_cut_x-0.05)+1.3;
$pdl_cut_y4=0.61/($pdl_cut_x-0.47)+1.19;
$pdl_cut_y3->inplace->clip(-2,1.1);
$pdl_cut_y4->inplace->clip(-2,1.1);
$pdl_cut_y3->inplace->setvaltobad(1.1);
$pdl_cut_y4->inplace->setvaltobad(1.1);
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII_AGNs=1.89*($pdl_cut_x)+0.76;
$pdl_cut_y_OI=0.73/(($pdl_cut_x+0.59))+1.33;#+1.10;
$pdl_cut_y_OI_AGNs=1.18*($pdl_cut_x)+1.30;
$pdl_cut_y_OIII_AGNs=1.14*($pdl_cut_x)+0.36;

$pdl_cut_y_SII->inplace->clip(-2,1.1);
$pdl_cut_y_OI->inplace->clip(-2,1.1);
$pdl_cut_y_SII->inplace->setvaltobad(1.1);
$pdl_cut_y_OI->inplace->setvaltobad(1.1);

$label1_1='log([NII]/H#ga';
$label1_2='log([SII]/H#ga)';
$label1_3='log([OI]/H#ga)';
$label2='log([OIII]/H#gb)';
$label_cb='log|EW(H#ga)|';

$plot_file="plBPT_mass_morph_res".$pldev_suf; #$dev="pdfcairo";



#$alpha=0.5;
plsdev($pldev_name); plsfnam($plot_file); 
plscolbga(255,255,255,1);
my $xp1 = 100.;
my $yp1 = 100.;
my $xleng1 = 800;
my $yleng1 = 800;
my $xoff1 = 10;
my $yoff1 = 20;
plspage($xp1, $yp1, $xleng1, $yleng1, $xoff1, $yoff1);
plinit();
plsfci($fci[3]);
my $n_ang=6;
my $size=0.05;
READ_cmap_plplot0_r($pl_colormap0,$alpha0);
READ_cmap_plplot1_r($pl_colormap1,$alpha);
pladv(0);


#
 # All masses
#

#
# All galaxies
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.095, 0.295, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
#$symbol=0;
plot_BPT_res_masked($pdl_cube,$pdl_fcube,0,"ALL");

#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,5,"E/S0");
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,14,"Sa/Sb");

print "DONE\n"; 
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,19,"Sc/Sd");



@Mass_limit=(7,9.5,10.5,11,13.5);
#@order=(7,2,11,16, 6,1,10,16, 9,4,13,18, 8,3,12,17);
$nord=0;
@order=(8,3,12,17,  9,4,13,18, 6,1,10,15, 7,2,11,16);

for (my $I=0;$I<$#Mass_limit;$I++) {
    my $min_M=$Mass_limit[$I];
    my $max_M=$Mass_limit[$I+1];
    my $y1_view=0.10+0.17*($I+1);
    my $y0_view=0.10+0.17*($I);
    print "$min_M<M<$max_M #=@dims\n";
    plscol0a(0,0,0,0,1);
    plcol0(0);
    plvpor (0.095, 0.295, $y0_view, $y1_view);
    plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcnstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
    }
    my $label_now="".$min_M."-".$max_M."";
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],$label_now); $nord++;


#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;
}

######################

#READ_cmap_plplot1_r("/home/sanchez/sda2/code/colortables/magma.csv",0.7);
READ_cmap_plplot1_r($pl_colormap1,0.7);
#my @stats=stats($pdl_color);
#my $zmin=$stats[3];
#my $zmax=$stats[4];
plcol0(0);
plcolorbar(0.90, 0.93, 0.10, 0.95,0,$zmin,$zmax,$label_cb,2.65,0,-2);
plvpor (0.13, 0.895, 0.10, 0.95);
pllab($label1_1,$label2,"");

#plmtex(3.2,0,-1.8,"l",$label2);
plend();

#exit;
print "BPT by distance\n";
#########################################################################
# BPT diagram PLPLOT Distance
#########################################################################

# reading BPT diagrams!
my @diag_files;
my @diag_files_EW;

$nf=0;
my @diag_files;
my @diag_files_EW;
open(DIR,"ls fitsfiles/*.MAP_OIII_NII.fits.gz |");
while($file=<DIR>) {
    chop($file);
    if ($file !~ "STDDEV") {
	$diag_files[$nf]=$file;
	$file =~ s/MAP/fMAP/g;
	$file =~ s/NII/NII_r/g;
	$diag_files_EW[$nf]=$file;
	print "$nf,$diag_files[$nf]\n";
	$nf++;
    }
}
close(DIR);
print "$nf DIAG files\n";

my $pdl_cube;
my $pdl_fcube;
for ($k=0;$k<$nf;$k++) {
	my $map=rfits($diag_files[$k]);
	my $fmap=rfits($diag_files_EW[$k]);
	if ($k==0) {
	    $h=$map->gethdr();
	    ($NX,$NY)=$map->dims();
	    $pdl_cube=zeroes($NX,$NY,$nf);
	    $pdl_fcube=zeroes($NX,$NY,$nf);
	    $crval1=$$h{CRVAL1}; 
	    $crval2=$$h{CRVAL2};
	    $cdelt1=$$h{CDELT1}; 
	    $cdelt2=$$h{CDELT2};
	}
#	print "$NX,$ny,$k\n";
	my $t = $pdl_cube->slice(":,:,($k)");
	$t .= $map;
	my $t = $pdl_fcube->slice(":,:,($k)");
	$t .= $fmap;#->log10;
}
print "DONE\n"; 

$zmin=0;
$zmax=3;
$pdl_fcube->inplace->clip(-0.1,$zmax);
$x_min1=$crval1;
$x_max1=$crval1+$NX*$cdelt1;
$y_min1=$crval2;
$y_max1=$crval2+$NY*$cdelt2;



$pdl_cut_x=-3+0.02*pdl(0..300);
$pdl_cut_y=-0.7+0.2-3.67*$pdl_cut_x;
$pdl_cut_y2=-1.7+0.5-3.67*$pdl_cut_x;
$pdl_cut_y3=0.61/($pdl_cut_x-0.05)+1.3;
$pdl_cut_y4=0.61/($pdl_cut_x-0.47)+1.19;
$pdl_cut_y3->inplace->clip(-2,1.1);
$pdl_cut_y4->inplace->clip(-2,1.1);
$pdl_cut_y3->inplace->setvaltobad(1.1);
$pdl_cut_y4->inplace->setvaltobad(1.1);
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII_AGNs=1.89*($pdl_cut_x)+0.76;
$pdl_cut_y_OI=0.73/(($pdl_cut_x+0.59))+1.33;#+1.10;
$pdl_cut_y_OI_AGNs=1.18*($pdl_cut_x)+1.30;
$pdl_cut_y_OIII_AGNs=1.14*($pdl_cut_x)+0.36;

$pdl_cut_y_SII->inplace->clip(-2,1.1);
$pdl_cut_y_OI->inplace->clip(-2,1.1);
$pdl_cut_y_SII->inplace->setvaltobad(1.1);
$pdl_cut_y_OI->inplace->setvaltobad(1.1);

$label1_1='log([NII]/H#ga';
$label1_2='log([SII]/H#ga)';
$label1_3='log([OI]/H#ga)';
$label2='log([OIII]/H#gb)';
$label_cb='log|EW(H#ga)|';

$plot_file="plBPT_mass_dist_res".$pldev_suf; #$dev="pdfcairo";


#$alpha=0.5;
plsdev($pldev_name); plsfnam($plot_file); 
plscolbga(255,255,255,1);
my $xp1 = 100.;
my $yp1 = 100.;
my $xleng1 = 800;
my $yleng1 = 210;
my $xoff1 = 10;
my $yoff1 = 20;
plspage($xp1, $yp1, $xleng1, $yleng1, $xoff1, $yoff1);
plinit();
plsfci($fci[3]);
my $n_ang=6;
my $size=0.05;
READ_cmap_plplot0_r($pl_colormap0,$alpha0);
READ_cmap_plplot1_r($pl_colormap1,$alpha);
pladv(0);


#
 # All masses
#

#
# All galaxies
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.095, 0.295, 0.17, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bncst", "bcnstv");
#$symbol=0;
plot_BPT_res_masked($pdl_cube,$pdl_fcube,0,"ALL");
#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, 0.17, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bncst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,5,"E/S0");
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, 0.17, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bncst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,14,"Sa/Sb");

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, 0.17, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bncst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,19,"Sc/Sd");



#READ_cmap_plplot1_r("/home/sanchez/sda2/code/colortables/magma.csv",0.7);
READ_cmap_plplot1_r($pl_colormap1,0.7);
#my @stats=stats($pdl_color);
#my $zmin=$stats[3];
#my $zmax=$stats[4];
$label_cb="R/Re";
plcol0(0);
plcolorbar(0.90, 0.93, 0.17, 0.95,0,$zmin,$zmax,$label_cb,2.65,0,-2);
plvpor (0.13, 0.895, 0.18, 0.95);
pllab($label1_1,$label2,"");


plend();

#################################################
# Mass Morph Dist Res
#################################################

print "#Mass_Morph_Dist_Res I\n";



$zmin=0;
$zmax=3.5;
$x_min1=$crval1;
$x_max1=$crval1+$NX*$cdelt1;
$y_min1=$crval2;
$y_max1=$crval2+$NY*$cdelt2;



$pdl_cut_x=-3+0.02*pdl(0..300);
$pdl_cut_y=-0.7+0.2-3.67*$pdl_cut_x;
$pdl_cut_y2=-1.7+0.5-3.67*$pdl_cut_x;
$pdl_cut_y3=0.61/($pdl_cut_x-0.05)+1.3;
$pdl_cut_y4=0.61/($pdl_cut_x-0.47)+1.19;
$pdl_cut_y3->inplace->clip(-2,1.1);
$pdl_cut_y4->inplace->clip(-2,1.1);
$pdl_cut_y3->inplace->setvaltobad(1.1);
$pdl_cut_y4->inplace->setvaltobad(1.1);
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII_AGNs=1.89*($pdl_cut_x)+0.76;
$pdl_cut_y_OI=0.73/(($pdl_cut_x+0.59))+1.33;#+1.10;
$pdl_cut_y_OI_AGNs=1.18*($pdl_cut_x)+1.30;
$pdl_cut_y_OIII_AGNs=1.14*($pdl_cut_x)+0.36;

$pdl_cut_y_SII->inplace->clip(-2,1.1);
$pdl_cut_y_OI->inplace->clip(-2,1.1);
$pdl_cut_y_SII->inplace->setvaltobad(1.1);
$pdl_cut_y_OI->inplace->setvaltobad(1.1);

$label1_1='log([NII]/H#ga';
$label1_2='log([SII]/H#ga)';
$label1_3='log([OI]/H#ga)';
$label2='log([OIII]/H#gb)';
$label_cb='log|EW(H#ga)|';

$plot_file="plBPT_mass_morph_dist_res_OLD".$pldev_suf; #$dev="pdfcairo";


#$alpha=0.5;
plsdev($pldev_name); plsfnam($plot_file); 
plscolbga(255,255,255,1);
my $xp1 = 100.;
my $yp1 = 100.;
my $xleng1 = 800;
my $yleng1 = 800;
my $xoff1 = 10;
my $yoff1 = 20;
plspage($xp1, $yp1, $xleng1, $yleng1, $xoff1, $yoff1);
plinit();
plsfci($fci[3]);
my $n_ang=6;
my $size=0.05;
READ_cmap_plplot0_r($pl_colormap0,$alpha0);
READ_cmap_plplot1_r($pl_colormap1,$alpha);
pladv(0);


#
 # All masses
#

#
# All galaxies
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.095, 0.295, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
#$symbol=0;
plot_BPT_res_masked($pdl_cube,$pdl_fcube,0,"ALL");
#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,5,"E/S0");
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,14,"Sa/Sb");

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked($pdl_cube,$pdl_fcube,19,"Sc/Sd");



@Mass_limit=(7,9.5,10.5,11,13.5);
#@order=(7,2,11,16, 6,1,10,16, 9,4,13,18, 8,3,12,17);
$nord=0;
@order=(8,3,12,17,  9,4,13,18, 6,1,10,15, 7,2,11,16);
for (my $I=0;$I<$#Mass_limit;$I++) {
    my $min_M=$Mass_limit[$I];
    my $max_M=$Mass_limit[$I+1];
    my $y1_view=0.10+0.17*($I+1);
    my $y0_view=0.10+0.17*($I);
    print "$min_M<M<$max_M #=@dims\n";
    plscol0a(0,0,0,0,1);
    plcol0(0);
    plvpor (0.095, 0.295, $y0_view, $y1_view);
    plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcnstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
    }
    my $label_now="".$min_M."-".$max_M."";
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],$label_now); $nord++;


#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked($pdl_cube,$pdl_fcube,$order[$nord],""); $nord++;
}

######################

#READ_cmap_plplot1_r("/home/sanchez/sda2/code/colortables/magma.csv",0.7);
READ_cmap_plplot1_r($pl_colormap1,0.7);
#my @stats=stats($pdl_color);
#my $zmin=$stats[3];
#my $zmax=$stats[4];
$label_cb="R/Re";
plcol0(0);
plcolorbar(0.90, 0.93, 0.10, 0.95,0,$zmin,$zmax,$label_cb,2.65,0,-2);
plvpor (0.13, 0.895, 0.10, 0.95);
pllab($label1_1,$label2,"");

#plmtex(3.2,0,-1.8,"l",$label2);
plend();


#################################################
# Mass Morph Dist Res
#################################################

print "#Mass_Morph_Dist_Res II\n";



#$zmin=0;
#$zmax=3.5;
$zmin=-0.49;
$zmax=2.65;
$x_min1=$crval1;
$x_max1=$crval1+$NX*$cdelt1;
$y_min1=$crval2;
$y_max1=$crval2+$NY*$cdelt2;



$pdl_cut_x=-3+0.02*pdl(0..300);
$pdl_cut_y=-0.7+0.2-3.67*$pdl_cut_x;
$pdl_cut_y2=-1.7+0.5-3.67*$pdl_cut_x;
$pdl_cut_y3=0.61/($pdl_cut_x-0.05)+1.3;
$pdl_cut_y4=0.61/($pdl_cut_x-0.47)+1.19;
$pdl_cut_y3->inplace->clip(-2,1.1);
$pdl_cut_y4->inplace->clip(-2,1.1);
$pdl_cut_y3->inplace->setvaltobad(1.1);
$pdl_cut_y4->inplace->setvaltobad(1.1);
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII=0.61/($pdl_cut_x-0.3)+1.3;
$pdl_cut_y_SII_AGNs=1.89*($pdl_cut_x)+0.76;
$pdl_cut_y_OI=0.73/(($pdl_cut_x+0.59))+1.33;#+1.10;
$pdl_cut_y_OI_AGNs=1.18*($pdl_cut_x)+1.30;
$pdl_cut_y_OIII_AGNs=1.14*($pdl_cut_x)+0.36;

$pdl_cut_y_SII->inplace->clip(-2,1.1);
$pdl_cut_y_OI->inplace->clip(-2,1.1);
$pdl_cut_y_SII->inplace->setvaltobad(1.1);
$pdl_cut_y_OI->inplace->setvaltobad(1.1);

$label1_1='log([NII]/H#ga';
$label1_2='log([SII]/H#ga)';
$label1_3='log([OI]/H#ga)';
$label2='log([OIII]/H#gb)';
$label_cb='log|EW(H#ga)|';

$plot_file="plBPT_mass_morph_dist_res".$pldev_suf; #$dev="pdfcairo";


#$alpha=0.5;
plsdev($pldev_name); plsfnam($plot_file); 
plscolbga(255,255,255,1);
my $xp1 = 100.;
my $yp1 = 100.;
my $xleng1 = 800;
my $yleng1 = 800;
my $xoff1 = 10;
my $yoff1 = 20;
plspage($xp1, $yp1, $xleng1, $yleng1, $xoff1, $yoff1);
plinit();
plsfci($fci[3]);
my $n_ang=6;
my $size=0.05;
READ_cmap_plplot0_r($pl_colormap0,$alpha0);
READ_cmap_plplot1_r($pl_colormap1,$alpha);
pladv(0);
#$color=255;
#$color1=205;
#$color2=152;

#$color=255;
#$color1=230;
#$color2=205;


$color=255;
$color1=220;
$color2=185;



#
 # All masses
#


#
# All galaxies
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.095, 0.295, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
#$symbol=0;
#my $pdl_fcube_EW=$pdl_fcube;
#my $pdl_fcube_ORG=$pdl_fcube;
plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,0,"ALL");

#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,5,"E/S0");
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,14,"Sa/Sb");

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, 0.78, 0.95);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,19,"Sc/Sd");



@Mass_limit=(7,9.5,10.5,11,13.5);
#@order=(7,2,11,16, 6,1,10,16, 9,4,13,18, 8,3,12,17);
$nord=0;
@order=(8,3,12,17,  9,4,13,18, 6,1,10,15, 7,2,11,16);
for (my $I=0;$I<$#Mass_limit;$I++) {
    my $min_M=$Mass_limit[$I];
    my $max_M=$Mass_limit[$I+1];
    my $y1_view=0.10+0.17*($I+1);
    my $y0_view=0.10+0.17*($I);
    print "$min_M<M<$max_M #=@dims\n";
    plscol0a(0,0,0,0,1);
    plcol0(0);
    plvpor (0.095, 0.295, $y0_view, $y1_view);
    plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcnstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcnstv");
    }
    my $label_now="".$min_M."-".$max_M."";
    plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,$order[$nord],$label_now); $nord++;

   
    
#
# E/S0
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.295, 0.495, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,$order[$nord],""); $nord++;

#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.495, 0.695, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,$order[$nord],""); $nord++;
#
# Sa/Sab
#
plscol0a(0,0,0,0,1);
plcol0(0);
plvpor (0.695, 0.895, $y0_view, $y1_view);
plwind ($x_min1,$x_max1,$y_min1,$y_max1);
    if ($I==0) {
	plbox (1.0, 0.25, 1.0, 0.25, "bcnst", "bcstv");
    } else {
	plbox (1.0, 0.25, 1.0, 0.25, "bcst", "bcstv");
    }
    plot_BPT_res_masked_EW($pdl_cube_ORG,$pdl_fcube,$pdl_fcube_EW,$order[$nord],""); $nord++;

    if ($I==0) {
	plcol0(1);
	my $pdl_x_box=pdl([$x_max1-0.37*($x_max1-$x_min1),$x_max1-0.03*($x_max1-$x_min1),$x_max1-0.03*($x_max1-$x_min1),$x_max1-0.37*($x_max1-$x_min1),$x_max1-0.37*($x_max1-$x_min1)]);
	my $pdl_y_box=pdl([$y_min1+0.05*($y_max1-$y_min1),$y_min1+0.05*($y_max1-$y_min1),$y_min1+0.68*($y_max1-$y_min1),$y_min1+0.68*($y_max1-$y_min1),$y_min1+0.05*($y_max1-$y_min1)]);
	plfill($pdl_x_box,$pdl_y_box);
	plcol0(0);
	plline($pdl_x_box,$pdl_y_box);
	plptex($x_max1-0.35*($x_max1-$x_min1),$y_min1+0.60*($y_max1-$y_min1),0,0,0,"R/Re");
	plcol0($color);
	plptex($x_max1-0.35*($x_max1-$x_min1),$y_min1+0.45*($y_max1-$y_min1),0,0,0,"<1");
	plcol0($color1);
	plptex($x_max1-0.35*($x_max1-$x_min1),$y_min1+0.3*($y_max1-$y_min1),0,0,0,"1-2");
	plcol0($color2);
	plptex($x_max1-0.35*($x_max1-$x_min1),$y_min1+0.15*($y_max1-$y_min1),0,0,0,">2");
	
    }


}





######################

#READ_cmap_plplot1_r("/home/sanchez/sda2/code/colortables/magma.csv",0.7);
READ_cmap_plplot1_r($pl_colormap1,0.7);
#my @stats=stats($pdl_color);
#my $zmin=$stats[3];
#my $zmax=$stats[4];

plcol0(0);
plcolorbar(0.90, 0.93, 0.10, 0.95,0,$zmin,$zmax,$label_cb,2.65,0,-2);
plvpor (0.13, 0.895, 0.10, 0.95);
pllab($label1_1,$label2,"");

#plmtex(3.2,0,-1.8,"l",$label2);
plend();

exit;

sub plot_BPT_res {
    my $pdl_cube=$_[0];
    my $pdl_fcube=$_[1];
    my $INDEX=$_[2];
    my $LABEL=$_[3];
#    print "$INDEX,$LABEL\n";
    $pdl_fden_BPT_NII=$pdl_fcube->slice(":,:,($INDEX)");
    $pdl_den_BPT_NII=$pdl_cube->slice(":,:,($INDEX)");
    
    $pdl_fden_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_fden_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_den_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->setbadtoval(0);
    plimagefr($pdl_fden_BPT_NII,$x_min1,$x_max1,$y_min1,$y_max1,$zmin,$zmax,$zmin,$zmax,0,0);
    my $dx=($x_max1-$x_min1)/$NX;
    my $dy=($y_max1-$y_min1)/$NY;
    @tr=($dx,0,$x_min1+$dx,0,$dy,$y_min1+$dy);
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_BPT_NII);
    my @stats_MAP=stats($pdl_den_BPT_NII);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_BPT_NII->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    plcol0(0);
    plwidth(2.0);
    plcont ($pdl_den_BPT_NII, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);
    plcol0(0);
    plwidth(2);
    pllsty(1);
    plline($pdl_cut_x,$pdl_cut_y4);
    pllsty(2);
    plline($pdl_cut_x,$pdl_cut_y3);
    pllsty(1);
    $pdl_cut_x_sec=$pdl_cut_x->where($pdl_cut_x>-0.1);
    $pdl_cut_y_OIII_AGNs_sec=$pdl_cut_y_OIII_AGNs->where($pdl_cut_x>-0.1);
    pllsty(3);
    plline($pdl_cut_x_sec,$pdl_cut_y_OIII_AGNs_sec);
    pllsty(1);
    plwidth(1);
    plcol0(0); plptex($x_min1+0.05*($x_max1-$x_min1),$y_min1+0.1*($y_max1-$y_min1),0,0,0,$LABEL);
    
    return 1;
}


sub plot_BPT_res_masked {
    my $pdl_cube=$_[0];
    my $pdl_fcube=$_[1];
    my $INDEX=$_[2];
    my $LABEL=$_[3];
#    print "$INDEX,$LABEL\n";
    $pdl_fden_BPT_NII=$pdl_fcube->slice(":,:,($INDEX)");
    $pdl_den_BPT_NII=$pdl_cube->slice(":,:,($INDEX)");
    
    $pdl_fden_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_fden_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_den_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->setbadtoval(0);
    my $kernel=ones(3,3)/9;
    $pdl_den_BPT_NII=conv2d($pdl_den_BPT_NII, $kernel, {Boundary => Reflect});

    my $dx=($x_max1-$x_min1)/$NX;
    my $dy=($y_max1-$y_min1)/$NY;
    @tr=($dx,0,$x_min1+$dx,0,$dy,$y_min1+$dy);
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_BPT_NII);
    my @stats_MAP=stats($pdl_den_BPT_NII);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_BPT_NII->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    $level_cut=$pdl_levels->at(90);
    my $pdl_my_levels=$pdl_den_BPT_NII->lclip($level_cut);
    $pdl_my_levels->inplace->setvaltobad($level_cut);
    $pdl_my_levels=$pdl_my_levels/$pdl_my_levels;
    $pdl_my_levels->inplace->setnantobad;
    $pdl_my_levels->inplace->setbadtoval(0);
    $pdl_fMAP=$pdl_fden_BPT_NII*$pdl_my_levels;
    $pdl_fMAP->inplace->setvaltobad(0);
    
#    (my $pdl_fMAP,my $pdl_my_levels)=where($pdl_fden_BPT_NII,$pdl_den_BPT_NII,$pdl_den_BPT_NII>$level_cut);
    
    #    plimagefr($pdl_fden_BPT_NII,$x_min1,$x_max1,$y_min1,$y_max1,$zmin,$zmax,$zmin,$zmax,0,0);
    plimagefr($pdl_fMAP,$x_min1,$x_max1,$y_min1,$y_max1,$zmin,$zmax,$zmin,$zmax,0,0);



    plcol0(0);
    plwidth(2.0);
    plcont ($pdl_den_BPT_NII, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);
    plcol0(0);
    plwidth(2);
    pllsty(1);
    plline($pdl_cut_x,$pdl_cut_y4);
    pllsty(2);
    plline($pdl_cut_x,$pdl_cut_y3);
    pllsty(1);
    $pdl_cut_x_sec=$pdl_cut_x->where($pdl_cut_x>-0.1);
    $pdl_cut_y_OIII_AGNs_sec=$pdl_cut_y_OIII_AGNs->where($pdl_cut_x>-0.1);
    pllsty(3);
    plline($pdl_cut_x_sec,$pdl_cut_y_OIII_AGNs_sec);
    pllsty(1);
    plwidth(1);
    plcol0(0); plptex($x_min1+0.05*($x_max1-$x_min1),$y_min1+0.1*($y_max1-$y_min1),0,0,0,$LABEL);
    
    return 1;
}



sub plot_BPT_res_masked_EW {
    my $pdl_cube=$_[0];
    my $pdl_fcube=$_[1];
    my $pdl_fcube_EW=$_[2];
    my $INDEX=$_[3];
    my $LABEL=$_[4];    

    $pdl_fden_BPT_NII=$pdl_fcube_EW->slice(":,:,($INDEX)");
    $pdl_fdist=$pdl_fcube->slice(":,:,($INDEX)");
    $pdl_den_BPT_NII=$pdl_cube->slice(":,:,($INDEX)");
    my @junk=$pdl_fdist->dims;
    $pdl_fden_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_fden_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_fden_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->clip(-1e12,1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(-1e12);
    $pdl_den_BPT_NII->inplace->setvaltobad(1e12);
    $pdl_den_BPT_NII->inplace->setnantobad;
    $pdl_den_BPT_NII->inplace->setbadtoval(0);
    my $kernel=ones(3,3)/9;
    $pdl_den_BPT_NII=conv2d($pdl_den_BPT_NII, $kernel, {Boundary => Reflect});

    my $dx=($x_max1-$x_min1)/$NX;
    my $dy=($y_max1-$y_min1)/$NY;
    @tr=($dx,0,$x_min1+$dx,0,$dy,$y_min1+$dy);
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_BPT_NII);
    my @stats_MAP=stats($pdl_den_BPT_NII);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_BPT_NII->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    my @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    $level_cut=$pdl_levels->at(90);
    my $pdl_my_levels=$pdl_den_BPT_NII->lclip($level_cut);
    $pdl_my_levels->inplace->setvaltobad($level_cut);
    $pdl_my_levels=$pdl_my_levels/$pdl_my_levels;
    $pdl_my_levels->inplace->setnantobad;
    $pdl_my_levels->inplace->setbadtoval(0);
    $pdl_fMAP=$pdl_fden_BPT_NII*$pdl_my_levels;
    $pdl_fMAP->inplace->setvaltobad(0);
    plimagefr($pdl_fMAP,$x_min1,$x_max1,$y_min1,$y_max1,$zmin,$zmax,$zmin,$zmax,0,0);

    plcol0(0);
    plwidth(2.0);

#    plcont ($pdl_den_BPT_NII, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);

    
    ################################################################################
    # R<Re;
    #   ($pdl_den_now,$pdl_diag_now,$pdl_fdist_now)=where($pdl_den_BPT_NII,$pdl_fden_BPT_NOW,$pdl_fdist,$pdl_fdist<0.5);



    my $cut=1.3;
    my $pdl_fdist_now=$pdl_fdist->hclip($cut);
    $pdl_fdist_now=$pdl_fdist_now*$pdl_my_levels;
    $pdl_fdist_now->inplace->setvaltobad($cut);
    $pdl_fdist_now=$pdl_fdist_now/$pdl_fdist_now;
    $pdl_fdist_now->inplace->setnantobad;
    $pdl_fdist_now->inplace->setbadtoval(0);
#    $pdl_fdist_now=$pdl_fdist_now*1.0;
    my $pdl_den_now=$pdl_den_BPT_NII/$pdl_fdist_now;
    $pdl_den_now->inplace->setnantobad;
    $pdl_den_now->inplace->setbadtoval(0);
    my $kernel=ones(5,5)/25;
    $pdl_den_now=conv2d($pdl_den_now, $kernel, {Boundary => Reflect});

    
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_now);
    my @stats_MAP=stats($pdl_den_now);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_now->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    #my @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    my @levels=($pdl_levels->at(65));
    plcol0($color);
    plwidth(2.0);
    plcont ($pdl_den_now, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);

    ################################################################################
    # 1<R<2Re;
    #   ($pdl_den_now,$pdl_diag_now,$pdl_fdist_now)=where($pdl_den_BPT_NII,$pdl_fden_BPT_NOW,$pdl_fdist,$pdl_fdist<0.5);
    my $cut0=1.2;
    my $cut1=2.1;
    my $pdl_fdist_now=$pdl_fdist->clip($cut0,$cut1);
    $pdl_fdist_now=$pdl_fdist_now*$pdl_my_levels;
    $pdl_fdist_now->inplace->setvaltobad($cut0);
    $pdl_fdist_now->inplace->setvaltobad($cut1);
    $pdl_fdist_now=$pdl_fdist_now/$pdl_fdist_now;
    $pdl_fdist_now->inplace->setnantobad;
    $pdl_fdist_now->inplace->setbadtoval(0);
#    $pdl_fdist_now=$pdl_fdist_now*1.0;
    my $pdl_den_now=$pdl_den_BPT_NII/$pdl_fdist_now;
    $pdl_den_now->inplace->setnantobad;
    $pdl_den_now->inplace->setbadtoval(0);
    my $kernel=ones(5,5)/25;
    $pdl_den_now=conv2d($pdl_den_now, $kernel, {Boundary => Reflect});
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_now);
    my @stats_MAP=stats($pdl_den_now);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_now->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    #    my @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    my @levels=($pdl_levels->at(65));
    plcol0($color1);
    plwidth(2.0);
    plcont ($pdl_den_now, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);

    ################################################################################
    # R>2Re;
    #   ($pdl_den_now,$pdl_diag_now,$pdl_fdist_now)=where($pdl_den_BPT_NII,$pdl_fden_BPT_NOW,$pdl_fdist,$pdl_fdist<0.5);
    my $cut0=1.8;
    my $cut1=3.5;
    my $pdl_fdist_now=$pdl_fdist->clip($cut0,$cut1);
    $pdl_fdist_now=$pdl_fdist_now*$pdl_my_levels;
    $pdl_fdist_now->inplace->setvaltobad($cut0);
    $pdl_fdist_now->inplace->setvaltobad($cut1);
    $pdl_fdist_now=$pdl_fdist_now/$pdl_fdist_now;
    $pdl_fdist_now->inplace->setnantobad;
    $pdl_fdist_now->inplace->setbadtoval(0);
#    $pdl_fdist_now=$pdl_fdist_now*1.0;
    my $pdl_den_now=$pdl_den_BPT_NII/$pdl_fdist_now;
    $pdl_den_now->inplace->setnantobad;
    $pdl_den_now->inplace->setbadtoval(0);
    my $kernel=ones(5,5)/25;
    $pdl_den_now=conv2d($pdl_den_now, $kernel, {Boundary => Reflect});
    my $pdl_levels=zeroes(100);
    my $ii=0;
    my $sum_all=sum($pdl_den_now);
    my @stats_MAP=stats($pdl_den_now);
    for (my $ii=99;$ii>-1;$ii--) {
	my $f1=$ii/100;
	for (my $jj=99;$jj>-1;$jj--) {
	    my $f2=$jj/100;
	    my $cut=$f2*$stats_MAP[4]->at(0);
	    my $pdl_MAP_cut=$pdl_den_now->lclip($cut);
	    $pdl_MAP_cut->inplace->setvaltobad($cut);
	    $pdl_MAP_cut->inplace->setbadtoval(0);
	    my $sum_now=sum($pdl_MAP_cut);
	    if ($sum_all!=0) {
		my $f_now=$sum_now/$sum_all;
		if (($f_now>$f1)&&($pdl_levels->at($ii)==0)) {
		    set($pdl_levels,$ii,$cut);
		}
	    } else {
		set($pdl_levels,$ii,0);
	    }
	}
    }
    #    my @levels=($pdl_levels->at(10),$pdl_levels->at(50),$pdl_levels->at(85));
    my @levels=($pdl_levels->at(65));
    plcol0($color2);
    plwidth(2.0);
    plcont ($pdl_den_now, 1, $NX, 1, $NX, pdl(@levels), \&mypltr, 0);



    plcol0(0);
    plwidth(2);
    pllsty(1);
    plline($pdl_cut_x,$pdl_cut_y4);
    pllsty(2);
    plline($pdl_cut_x,$pdl_cut_y3);
    pllsty(1);
    $pdl_cut_x_sec=$pdl_cut_x->where($pdl_cut_x>-0.1);
    $pdl_cut_y_OIII_AGNs_sec=$pdl_cut_y_OIII_AGNs->where($pdl_cut_x>-0.1);
    pllsty(3);
    plline($pdl_cut_x_sec,$pdl_cut_y_OIII_AGNs_sec);
    pllsty(1);
    plwidth(1);
    plcol0(0); plptex($x_min1+0.05*($x_max1-$x_min1),$y_min1+0.1*($y_max1-$y_min1),0,0,0,$LABEL);
    
    return 1;
}
