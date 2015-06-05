#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($codonfile,$fastafile,$c1,$c2,$outfile,$cutoff,$tag,$listfile,$protein);
my ($aminoacid,$ignorefile,@expressions);
my $verbose = 1 ;
GetOptions(
            "tag=s"=>\$tag ,
            "codonfile=s"=>\$codonfile ,
            "fastafile=s"=>\$fastafile ,
            "c1=s"=>\$c1 ,
            "c2=s"=>\$c2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -c1 ") if(!defined $c1);
usage( "Need to give a input file name => option -c2 ") if(!defined $c2);
usage( "Need to give a input file name => option -tag ") if(!defined $tag);

my $codontable = util_getCodonTable();



my ($info1,$mapValues1) = util_ReadCodonBiasFile($c1);
my ($info2,$mapValues2) = util_ReadCodonBiasFile($c2);

print "Diffing\n";

my $ofhabsdiff = util_write("DIFFS/absdiff.$tag");
my $ofhatgcdiff = util_write("DIFFS/atgcdiff.$tag");
my $FINALABSSUM = 0 ;
foreach my $k (sort keys %{$mapValues1}){
my $atsumthird = 0 ;
my $gcsumthird = 0 ;
	my $t1 = $mapValues1->{$k};
	my $t2 = $mapValues2->{$k};
	my $sum = 0;

	my $A = 0;
	my $B = 0;
	foreach my $c (sort keys %{$t1}){
		my $v1 = $t1->{$c};
		my $v2 = $t2->{$c};
		my $absdiff = util_format_float(abs($v1 - $v2),1);
		$sum  = $sum + $absdiff ;

		my $diff = util_format_float(($v1 - $v2),1);
	    if($c =~ /..T/ || $c =~ /..A/){
			$atsumthird = util_format_float($atsumthird + $diff,1) ;
		    $A = $A+$v1 ;
		    $B = $B+$v2 ;
		}
		else{
			$gcsumthird = util_format_float($gcsumthird+$diff,1) ;
		}
	}
	my $D = util_format_float($A - $B,1) ;
	if(abs($D) < 0.1){
		$D = 0 ;
	}
	print $ofhabsdiff "$k $sum \n";
	$FINALABSSUM = $FINALABSSUM + $sum ;

	#print "$k $atsumthird $gcsumthird \n";
	print $ofhatgcdiff "($k,$D) \n" if($k ne "M" && $k ne "W");
}
print $ofhabsdiff "%(SUMOFDIFFS, $FINALABSSUM ) \n";




sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
