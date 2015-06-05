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
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $postfix = "ALL.1.fasta";
my $what = "nr";
GetOptions(
            "postfix=s"=>\$postfix ,
            "protein=s"=>\$protein ,
            "what=s"=>\$what ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $PWD = cwd;

my $out = "$BLASTOUT/$protein.blast.$what";
if(-e $out){
    my $v = utils_parseBlastOut($out);
    print "$out exists, so not running blast\n" ;
    print "$v\n";
    exit ;
}


my $EXEC = $what eq "nr" || $what eq "pdbaa" ? "$SRC/BLAST/blastp":"$SRC/BLAST/blastn";
my $exec = "$EXEC -db $what -query $FASTADIR/$protein.$postfix -out $out";
print "Running $exec\n";
system($exec);



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
