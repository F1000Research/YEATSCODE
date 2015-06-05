#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($forR,$reverse,$infile,$cutoff,$outfile,$which_tech,$seperator,$listfile);
my (@idx);
$seperator = ",";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "cutoff=s"=>\$cutoff ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "forR=s"=>\$forR ,
            "seperator=s"=>\$seperator ,
            "idx=i"=>\@idx,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
die if(!@idx);

my $idx1 = $idx[0];
my $idx2 = $idx[1];
my $ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();

my @X = util_R_extractIndex($infile,"X",$idx1);
my @Y = util_R_extractIndex($infile,"Y",$idx2);
my $NX = @X ;
my $NY = @X ;
print "$NX $NY\n";
die if($NX ne $NY);

my $XSTR = "X = c ( " ;
my $YSTR = "Y = c ( " ;
my $seperator = "";
my $first = 1 ;
my $CNTignored = 0 ;
while(@X){
	my $x = shift @X ;
	my $y = shift @Y ;
	if($x < $cutoff && $y < $cutoff){
		$CNTignored++;
		next ;
	}
	if($first){
		$seperator = "";
		$first = 0;
	}
	else{
		$seperator = ",";
	}
	$XSTR = $XSTR .  "$seperator $x";
	$YSTR = $YSTR .  "$seperator $y";
}

print "ignored $CNTignored\n";
print $ofh "$XSTR ) \n";
print $ofh "$YSTR ) \n";




print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
