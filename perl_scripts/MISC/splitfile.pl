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
my ($ofh,$infile,$outfile,$which_tech,$number,$listfile,$protein);
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
            "number=i"=>\$number ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -number ") if(!defined $number);
my $ifh = util_read($infile);

my $info = {};
my $CNT = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 $CNT++ ;
}
close($ifh);

my $smallcount = $CNT/$number;
$ifh = util_read($infile);
$CNT = 0 ;
my $index = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 if($CNT eq 0){
         $ofh = util_write($infile . ".$index");
		 $index++;
	 }
	 print $ofh $_ ; 
	 $CNT++ ;
	 if($CNT > $smallcount){
	 	$CNT = 0 ; 
	 }
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
