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
my ($infile,$p1,$ignoreband,$cutoff,$p2,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
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
            "ignoreband=i"=>\$ignoreband ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -ignorefile ") if(!defined $ignorefile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my @ignorelist= util_read_list_sentences($ignorefile);
my $ignorelist = {};
map { s/\s*//g ; $ignorelist->{lc($_)} = 1 ; } @ignorelist ;
#my @NAMES = qw(CE CI CK EM FL HC HL HP HU IF LE LM LY PK PL PT RT SE TZ VB);
my @NAMES = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);


my @FFF ; 
my @FFFNEG ; 
foreach my $n (@NAMES){
	my $ooo = "STAT/SPLIT.$n";
	my $ofhmm = util_write($ooo);
	push @FFF, $ofhmm;

	$ooo = "STAT/NEGSPLIT.$n";
	$ofhmm = util_write($ooo);
	push @FFFNEG, $ofhmm;
}

my $POS = {};
my $NEG = {};
my $MERGED = {};
my $IGN = 0 ;
my $NAME2TOTALCNT = {};
my $cnt = 0 ;
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;
	 if(exists $ignorelist->{$name}){
	 	#print "Ignoring not existent $name\n";
		$IGN++;
		next ;
	 }
	 my $N = @NAMES -1 ;

	 my $FINALN = 0 ; 
	 my $TOT = 0 ;
	 my @codesused ; 
	 my @VV ; 
	 foreach my $idx (0..$N){
	 	my $code = $NAMES[$idx];
		my $a = shift @l ;
		#push @codesused, $code ;
		push @codesused, $a ;
		$TOT = $TOT + $a ;
		if($a){
			$FINALN++;
			push @VV, $a ;
			my $FH = $FFF[$idx];
			print $FH "$name $a\n";
		}

	 }
	 if(!@VV){
		print "Nothing here for $name, so nexting\n" if($verbose);
	 	next ;
	 }

	 my ($mean,$sd) = util_GetMeanSD(\@VV);

	 $mean = int($mean);
	 $sd = int($sd);
	 $name = uc($name);

	 $sd = 1 if($FINALN eq 1);
	 $sd = 1 if($sd eq 0);

	 if(!$sd){
	    print  "$name $TOT $FINALN $mean $sd\t@codesused \n";
	 	die ;
	 }
	 my ($mean2sdratio) = int(100*($mean/$sd));
	 
	 next if(defined $cutoff && $mean < $cutoff);
	 $cnt++;
	 print $ofh "$name $TOT $FINALN $mean $sd $mean2sdratio\t@codesused \n";


	 die "@l" if(@l);
}
print "Wrote $cnt kkkk\n";



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
