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
my ($f1,$f2,$het2pdb,$pdbseqres,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $size ;
my $verbose = 0 ;
GetOptions(
            "het2pdb=s"=>\$het2pdb ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
print "Writing to $outfile - the commands to convert premon.in's\n";
my $ofh = util_append($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);

my (@l) = (keys %{$info});
if(@l ne 2){
	print $ofh "something not right $infile \n";
	die ;
}
my $A = shift @l ;
my $B = shift @l ;

my $a = $info->{$A};
my $b = $info->{$B};

my $allgood = 1 ;
my $subset = 0 ;
if($a ne $b){
	if($a =~ /$b/){
       $subset = 1 ;
	   print $ofh "$a $b\n";
	}
	elsif($b =~ /$a/){
       $subset = 2 ;
	   print $ofh "$a $b\n";
	}
	else{
		$allgood = 0 ;
	}
}

if(!$allgood){
	print $ofh "$infile is not the same \n";
}
else{
	if($subset){
	    print $ofh "$infile subset $subset \n";
	}
}
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
