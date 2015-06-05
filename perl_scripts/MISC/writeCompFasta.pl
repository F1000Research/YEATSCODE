#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
  use Time::HiRes qw( usleep ualarm gettimeofday tv_interval
   clock_gettime clock_getres  clock
   );


use PDB;
use Atom;
use Residue;

use POSIX qw(floor);
use Math::Combinatorics;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($pdb1,$pdb2,$infile,$outfile,$atomidx,$dontrunpymol);
my ($interactive,$annotate,$keepname,$dist,$mutatefile,$maxresults,$inconf,$outconf,$resultfile,$checkself);
my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my $MINDIST = 2 ;
$, = "  ";
GetOptions(
            "pdb1=s"=>\$pdb1 ,
            "infile=s"=>\$infile ,
            "mutatefile=s"=>\$mutatefile ,
            "keepname"=>\$keepname ,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give infile ") if(!defined $infile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();

my ($fe,$x) =  util_readfasta($infile);
chomp $x ;

$outfile = $infile . ".comp.fasta" if(!defined $outfile);
my $ofh = util_write($outfile);


my $rev = util_getComplimentaryString($fe) ;

my $nm = defined $keepname ? "$x rev" : "${x}_rev" ;

print $ofh "$nm\n";
print $ofh "$rev\n";


print "Wrote complimentary DNA seq in $outfile\n";



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
