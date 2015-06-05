#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($p1,$p2,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $maxdist = 3 ;
my $verbose = 1 ;
my ($moveZ,$verify,$decaaf,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "maxdist=i"=>\$maxdist ,
            "moveZ=i"=>\$moveZ ,
            "decaaf=i"=>\$decaaf ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofhclose = util_write("log.close");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a protein 2 id -option -p2  ") if(!defined $p2);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my ($D) = util_AlignAndMatchRemaining($p1,$p2,$infile,$maxdist);  


sub parseSingleLine{
	my ($line) = @_ ; 
	my ($num,$restype,$resnum,$atom,$x,$y,$z) = split " " , $line ; 
	return ($num,$restype,$resnum,$atom,$x,$y,$z);
}
print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
sub getResultsLine{
    my ($in,$p1,$p2) = @_ ; 
    my $ifh = util_read($in);
    my @l ; 
	my ($a1,$a2);
    while(<$ifh>){
	    next if(/RESULT/);
		print $_ ;
		if(!defined $a1){
		    $a1 = $p1->ParseResultLine($_) ;
		    next ;
		}
		if(!defined $a2){
		    $a2 = $p2->ParseResultLine($_) ;
		     next;
		}
    }
    return ($a1,$a2);
}

