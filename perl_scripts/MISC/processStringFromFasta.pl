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
my $cutoff ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
my $ifh = util_read($infile);
my $info ={};
my $pdbnm = "";

my $CNT = 0;
while(<$ifh>){
     
     chop ;
     next if(/^\s*$/);
     if(/^>\s*/){
	 	$pdbnm = $_ ;
		next ;
	 }
	 my ($nm,$junk) = split ; 
	 $CNT++;
	$info->{$nm} = $pdbnm ;
}


my $N = (keys %{$info});
print "There were $N unique in $CNT ones\n";
foreach my $k (keys %{$info}){
		my $pdbnm = $info->{$k};
		my $len = length($k);
		next if($len < $cutoff);
		print $ofh "$pdbnm \n";
		print $ofh "$k\n\n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
