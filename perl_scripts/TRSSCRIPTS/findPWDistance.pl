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
my ($infile,$p1,$mapfile,$cutoff,$p2,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=i"=>\$cutoff ,
            "mapfile=s"=>\$mapfile ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -ignorefile ") if(!defined $ignorefile);


$outfile = "$infile.PW";
print "Writing to $outfile\n";
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $meaninfo = {};
my $sdindfo = {};

my @ignorelist= util_read_list_sentences($ignorefile);
my $ignorelist = {};
map { s/\s*//g ; $ignorelist->{lc($_)} = 1 ; } @ignorelist ;
#my @NAMES = qw(CE CI CK EM FL HC HL HP HU IF LE LM LY PK PL PT RT SE TZ VB);
my @NAMES = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);
		

		my $ofhpw = util_write("pairwisenames");
while(@NAMES){
	my $a = shift @NAMES ;
	foreach my $i (@NAMES){
		print $ofhpw "$a-$i\n";
	}
}

my @NUMS ;
my $firstone = 1 ;
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;
	 if(exists $ignorelist->{$name}){
	 	#print "Ignoring not existent $name\n";
		next ;
	 }
	 my $N = @l - 1;
	 foreach my $idx (0..$N){

		my $CNT = 0 ;
		while(@l){
		    my $a = shift @l ;
			foreach my $i (@l){
				my $diff = abs ($a - $i);
				if($firstone){
					push @NUMS, $diff ;
				}
				else{
					$NUMS[$CNT] = $NUMS[$CNT] + $diff ;
				}
			    $CNT++;
			}
		}
	 }
	 $firstone = 0 ;
}


my $NN = @NUMS ;
print "There were $NN items\n";
$, = "\n";
print $ofh @NUMS ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
