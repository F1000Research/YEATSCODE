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
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


my $info = {};
my $prev ;
my $total = 0;
while(<$ifh>){
	chomp;
     next if(/^\s*$/);
	 if(/(R|K)...(R|K)/){
	 	my (@l) = /((R|K)...(R|K))/g;
		my $N = @l ;
        my $CNT = 0 ; 

		my @ll ;
		foreach my $l (@l){
			next if(!($l =~ /(E|D).(E|D)/));
			my $len = length($l);
			next if ($len < 5);
			push @ll, $l ;
			$CNT++;
		}

		next if(!$CNT);
		$total++;
		print $ofh "$prev ";
		foreach my $l (@ll){
			print $ofh "$l "
		}
		print $ofh "\n";
	 }
	 $prev = $_ ;
}

print "found $total\n";
close($ifh);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
