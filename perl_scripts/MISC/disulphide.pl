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
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
use Scalar::Util qw(looks_like_number);


use Math::Trig;




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($cutoff,$exception,$ann,$config,$p1,$p2,$infile,$score,$outfile,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;
my $DBCA ;
my $DBCB ;
my $DBCN ;

$cutoff = 0.012 ;
my $MINSAMPLE = 1 ;
my $NUMBEROFCTERMINALIGNORE = 0 ; 
my $NUMBEROFNTERMINALIGNORE = 0 ; 
my $ignorepro = 0 ;

my ($verify,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "protein=s"=>\$protein ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "cutoff=f"=>\$cutoff ,
            "config=s"=>\$config,
            "score=s"=>\$score,
            "ignorepro=i"=>\$ignorepro,
            "radii=i"=>\$radii ,
            "exception=s"=>\$exception ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a protein -option -protein  ") if(!defined $protein);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();

    my $i = $protein ;
    my @proteins ; 
    push @proteins, $i ; 

    my $PWD = cwd;

    my $pdb = "$PDBDIR/$protein.pdb";
    my $pdb1 = new PDB();
    $pdb1->ReadPDB($pdb);
    
	my $cnt = util_GetDisulphide ($protein,$pdb1,1);
	print "There were $cnt disulphide bonds\n";
    

    

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}


