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
my ($justseqeunce,$aalist,$cutoff,$outfile,$exception,$ann,$config,$p1,$p2,$infile,$score,$ignorepro,,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;




GetOptions(
            "which_tech=s"=>\$which_tech ,
            "p1=s"=>\$p1 ,
            "aalist=s"=>\$aalist ,
            "p2=s"=>\$p2 ,
            "outfile=s"=>\$outfile ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "protein=s"=>\$protein ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "cutoff=f"=>\$cutoff ,
            "config=s"=>\$config,
            "justseqeunce=s"=>\$justseqeunce,
            "score=s"=>\$score,
            "exception=s"=>\$exception ,
           );


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;

usage( "Need to give protein => option -protein ") if(!defined $protein);

die "Dont recognize command line arg @ARGV " if(@ARGV);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);


system ("touch $protein.helixvalues") if(! -e "$protein.helixvalues");
my $donefh = util_append("$protein.helixvalues");
my $ofhcircular = util_write("$protein.circularfasta");
my $writetex = 1 ;
util_helixWheel($protein,$pdb1,$donefh,$ofhcircular,$writetex,$justseqeunce);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}


