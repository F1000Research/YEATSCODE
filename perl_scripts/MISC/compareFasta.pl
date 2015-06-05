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
my ($p1,$p2,$pdb1,$pdb2,$infile,$outfile,$atomidx,$dontrunpymol);
my ($interactive,$annotate,$dist,$mutatefile,$maxresults,$inconf,$outconf,$resultfile,$checkself);
my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my $MINDIST = 2 ;
$, = "  ";
GetOptions(
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "mutatefile=s"=>\$mutatefile ,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();


usage( "Need to pdb name => option -p1 ") if(!defined $p1);
usage( "Need to pdb name => option -p2 ") if(!defined $p2);


system ("touch $outfile") if(! -e $outfile);
my $ofh = util_append($outfile);

my $seq1 = util_GetFasta($p1);
my $seq2 = util_GetFasta($p2);
if($seq1 ne $seq2){
	#print "Seq1= $seq1 \n";
	#print "Seq2= $seq2 \n";
	print $ofh "diff $p1 $p2\n";
}
else{
	print $ofh "same $p1 $p2\n";
}





sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
