#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($start,$ending,$name,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "name=s"=>\$name ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofhend = util_open_or_append("END.fasta");
my $ofhbegin = util_open_or_append("BEGIN.fasta");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my ($str,$firstline) = util_readfasta($infile);
my $len = length($str);


my @l = split " ", $firstline ;
$name = $l[0];

my $end = util_findLength($str,100,"end");
my $begin = util_findLength($str,100,"begin");

print $ofhend "${name}end\n";
print $ofhend "$end\n";

print $ofhbegin "${name}begin\n";
print $ofhbegin "$begin\n";

chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
