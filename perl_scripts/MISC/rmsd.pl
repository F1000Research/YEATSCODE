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
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
my @l1 = qw (2.017 3.109 6.715 4.889 7.599 5.970);
my @l2 = qw (1.994 4.513 4.967 5.319 6.534    5.854);
my $N = @l1 ;


my $sum = 0 ;
while(@l1){
	my $a = shift @l1 ;
	my $b = shift @l2 ;

	my $diff = $a-$b;
	my $diffsq = $diff*$diff;
	$sum = $sum + $diffsq;

}

my $tmp = $sum/$N;
my $rmsd = sqrt($tmp);

print "RMDS = $rmsd\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
