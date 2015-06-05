#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;
use Primer;
  use Bio::Tools::CodonTable;
 my $myCodonTable   = Bio::Tools::CodonTable->new();
 

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my ($to,$idx) ;
my $howmany = 100000 ;
my $verbose = 0 ;
my $TEMP = 78 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "to=s"=>\$to ,
            "idx=i"=>\$idx ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
util_CmdLine("infile",$infile);
util_CmdLine("outfile",$outfile);
my $ofh = util_write($outfile);
my $ofhwarning = util_append("warning.stop");
my $ifh = util_read($infile);
my $info = {};
my $line = "";

my $preheader = "";
while(<$ifh>){
	chop ; 
     next if(/^\s*$/);
	 if(/\s*>/){
	 	die "found two preheader" if($preheader ne "");
		$preheader = $_ ; 
		next ;
	 }
	 s/\s*//g; 
	$line = $line . $_;
}

die "Expected nucleotide seq " if(IsAA($line));


my @l = ($line =~ /(...)(...)(...)/g);

my $warned = 0 ; 
while(@l){
	my $l = shift @l ;
   my $len = length($l);
   die "wrong len $len" if($len ne 3);
   my  $aaName = $myCodonTable->translate($l);
   #print "$l";
   if ($aaName eq "*" && @l){
   	   #print "$aaName \n";
       print "Warning: $infile found stop codon\n" if(!$warned);
       print $ofhwarning "Warning: $infile found stop codon\n" if(!$warned);
	   $warned = 1; 
   }
   print $ofh "$aaName";
}
print $ofh "\n";

