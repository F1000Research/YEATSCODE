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
use ConfigPDB;


use POSIX qw(floor);
use Math::Combinatorics;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($pdb1,$pdb2,$infile,$outfile,$atomidx,$dontrunpymol);
my ($anndir,$interactive,$annotate,$aalist,$dist,$findresidues,$maxresults,$inconf,$outconf,$resultfile,$checkself);
my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my $MINDIST = 2 ;
$, = "  ";
my $LEN = 10 ; 
GetOptions(
            "pdb1=s"=>\$pdb1 ,
            "anndir=s"=>\$anndir ,
            "interactive"=>\$interactive ,
            "checkself"=>\$checkself ,
            "dontrunpymol"=>\$dontrunpymol ,
            "findresidues"=>\$findresidues ,
            "atomidx=s"=>\$atomidx ,
            "resultfile=s"=>\$resultfile ,
            "maxresults=i"=>\$maxresults ,
            "dist=f"=>\$dist ,
            "len=i"=>\$LEN ,
            "infile=s"=>\$infile ,
            "aalist=s"=>\$aalist ,
            "grpconfig=s"=>\$grpconfig ,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a aalist => option -aalist ") if(!defined $aalist);
#usage( "Need to give a residue number") if(!defined $atomidx);
my $ofh = util_write($outfile);
my $ifh = util_read($infile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

ConfigPDB_Init($grpconfig);

my ($aainfo,$grplist,$map3to1,$map1to3) = util_ParseAAGroups($aalist);

my $info = {};
while(<$ifh>){
	 my ($nm,$junk) = split ; 
	 my $NM = $map3to1->{$nm} ; 
	 my $mappednm = $aainfo->{$NM} or die ; 
	 if(!defined $info->{$mappednm}){
		$info->{$mappednm} = 0 ; 
	}
	$info->{$mappednm} = $info->{$mappednm} + 1 ;
}

my $cnt = 1 ;
foreach my $k (sort keys %{$info}){
     $k = uc($k);
	 my $val = $info->{$k} ; 
	 print $ofh "$cnt  $val \n";
	 $cnt++;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;

}
