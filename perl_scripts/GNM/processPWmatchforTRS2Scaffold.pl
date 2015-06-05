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
my ($infile,$f1,$f2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($outdir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "outdir=s"=>\$outdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a listfile -option -f1  ") if(!defined $f1);
usage( "Need to give a listfile -option -f2  ") if(!defined $f2);
usage( "Need to give a outdir pdb id -option -outdir  ") if(!defined $outdir);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $ifh = util_read($infile);

while(<$ifh>){
	my ($i,$j) = split ;
	print $ofh "myblastcompare2Fastafiles.csh $f1/$i.ALL.1.fasta $f2/$j.ALL.1.fasta $outdir/$i.pw N \n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
