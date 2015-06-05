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
$outfile = "$infile.resultstyle";
my $ofh = util_write($outfile);
print "writing to $outfile\n";

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $info = {};
while(<$ifh>){
     next if(/^\s*$/);
	 s/,/   /g;
	 my (@l) = split ; 
	 foreach my $i (@l){
		my ($name,$number,$type) = ($i =~ /([a-zA-Z]+)([0-9]+)(.*)/);
		print $ofh "$name/$number/$type ";
	 }
}
close($ifh);
print $ofh "\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
