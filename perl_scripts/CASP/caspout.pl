#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($model,$target,$infile,$p1,$p2,$model1,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "model=s"=>\$model ,
            "target=s"=>\$target ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a p1 pdb id -option -p1  ") if(!defined $p1);
usage( "Need to give a p2 pdb id -option -p2  ") if(!defined $p2);
usage( "Need to give a target pdb id -option -target  ") if(!defined $target);
usage( "Need to give a model pdb id -option -model  ") if(!defined $model);

print "PFRMAT QA\n";
print "TARGET $target\n";
print "AUTHOR 3663-5056-2522\n";
print "METHOD Escapist and Proqud - http://f1000research.com/articles/2-243 and  http://f1000research.com/articles/2-211/v3\n";
print "MODEL $model\n";
print "QMODE 1\n";
print "$p1 1 \n";
print "$p2 0.9 \n";
print "END\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
