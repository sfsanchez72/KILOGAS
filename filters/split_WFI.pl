#!/usr/bin/perl

$n=0;
open(FH,"<WFI_filters.txt");
while($line=<FH>) {
    chop($line);
    $lines[$n]=$line;
    $n++;
}
close(FH);

$nf=0;
$count=0;
for ($i=0;$i<$n;$i++) {
    @data=split(" ",$lines[$i]);
#    print "$lines[$i] $count\n";
    if ($count==1) {
	if ($data[0]!=1) {
#	    print "$ni[$nf]\n";
	    $NI++;
	} else {
	    @DAT=split(" ",$lines[$i-2]);
	    $ni[$nf-1]=$DAT[0];
	#    $i=$i-1;
	    $count=0;
	    print "$start[$nf] $ni[$nf]\n";
#	    <stdin>;
	}
    }


    if (($data[0]==1)&&($count==0)) {
	$namef[$nf]=$lines[$i-1].".txt";
	$namef[$nf] =~ s/\+//;
	$namef[$nf] =~ s/ /_/g;
	$start[$nf]=$i;
	$count=1;
	$NI=0;
	print "$nf $namef[$nf] ";
	$nf++;
    }

}
print "PASO $ni[5]\n";
print "NF=$nf\n";
for ($i=0;$i<$nf;$i++) {
    print "$i $namef[$i] $ni[$i] $start[$i]\n";
    open(OUT,">$namef[$i]");
    for ($j=0;$j<$ni[$i];$j++) {
	$k=$start[$i]+$j;
	@data=split(" ",$lines[$k]);
#	print "$lines[$k] @data\n";
#	print "$data[0] $data[1] $data[2]\n";
	print OUT "$data[0] $data[1] $data[2]\n";
    }
    close(OUT);
}



exit;
