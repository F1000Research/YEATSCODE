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
my ($infile,$outfile,$which_tech,$fpocket,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "fpocket=s"=>\$fpocket ,
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
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

my $list = {};
my @list ;
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

foreach my $i (@list){
	print "$i\n";
}

#my $Match = "GLU/205/OE1 GLU/206/OE1 SER/630/OG TYR/662/OH";
my $Match = "GLU/205/OE1 GLU/206/OE1 TYR/662/OH";
my $info = {};


my ($info,$infoperatom) ;
if(defined $fpocket){
    ($info,$infoperatom) = parseFPocketFile($fpocket);
}


my $cnt = 0 ;
while(<$ifh>){
   my $orig = $_ ;
   next if(/^\s*$/);
   chomp ;
	 $cnt++ ;
	 $outfile = "pymol.$protein.in.$cnt";
     next if(/^\s*$/);
	 s/@.*// ;
	 s/#/ /g;
	 s/-/ /g;
	 my @l = split ; 

	 my $cntpocket = 0 ; 
	 if(defined $fpocket){
		foreach my $i (@l){
			if(exists $infoperatom->{$i}){
				$cntpocket++;
			}
		}
	 }


	 next if(defined $fpocket && $cntpocket < 2);
	 #if($cntpocket eq 4){
	 	print "$cnt $orig";
	 #}

     print "cntpocket = $cntpocket \n";

	 #my $NEW = "$l[0] $l[1] $l[3]";

     my $ofhTMP = util_write($outfile);
     my $ofhTMPpot = util_write("$outfile.pot.in");

	 print $ofhTMP "$_ \n" ;
	 print $ofhTMPpot "$_ \n" ;

	 print $ofhTMP "$Match\n";
	 print $ofh "alignProteins.pl -out $outfile.p1m -in $outfile -p1 1PTGA -p2 3W2TA \n";
	 print $ofh "printPairwise.pl -out $outfile.pot.out -c \$CONFIGGRP -ra 222 -pr $protein -in $outfile.pot.in\n";

	 close($ofhTMP);
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
