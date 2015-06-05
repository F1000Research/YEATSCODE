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
my ($infile,$outfile,$which_tech,$listfile);
my (@expressions);
my $idx = 0;
$outfile = "ttt";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "idx=i"=>\$idx ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my $ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();


my @numbers  ;
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ;
	 push @numbers ,  $l[$idx] ;
}
close($ifh);


my @sorted= sort { $a <=> $b } (@numbers);
my $N = @sorted ; 
my $MIN = $sorted[0];
my $MAX = $sorted[$N - 1];

my ($m, $sd) = util_GetMeanSD(\@numbers);
my $upper = $m + $sd ; 
my $lower = $m - $sd ; 
print "mean sd upper lower \n";
my $total = $m* $N ; 
print "$total $m $sd upper=$upper lower=$lower. There are $N numbers, max = $MAX and min = $MIN. Total = $total \n";


chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
