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
my $ofh = util_append($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $ifh = util_read($infile);

my $info = {};
my @l ;
while(<$ifh>){
	 if(/^\s*POINTS/){
	 	@l = split ;
		shift @l ;
	 }
	 if(/^\s*DIST/){
	 	my ($junk,$a,$b,$j,$dist) = split ;
		$info->{$a} = {} if(! defined $info->{$a});
		$info->{$b} = {} if(! defined $info->{$b});

		$info->{$a}->{$b} = $dist ;
		$info->{$b}->{$a} = $dist ;
	}
}

my @vals ;
my $TAB = {};
foreach my $k (sort keys %{$info}){
	my $v = $info->{$k} ; 
	my @V = (sort  { $a <=> $b } values %{$v});
	my $N = @V - 1;
	my $lowval = $V[$N];
	#print "$k $lowval \n";
	$TAB->{$lowval} = $k ;
	push @vals, $lowval ; 
}

my @sortedvals = sort  { $a <=> $b } @vals ;
my $MIN = $sortedvals[0] ;

my $OOO = {};

$OOO->{a} = 0 ; 
$OOO->{b} = 1 ; 
$OOO->{c} = 2 ; 
$OOO->{d} = 3 ; 

my $IDX = $OOO->{$TAB->{$MIN}} ; 
my $AA = $l[$IDX] ; 

print "MMMMM =$MIN $TAB->{$MIN} $AA\n";
print $ofh "$protein $MIN $AA \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
