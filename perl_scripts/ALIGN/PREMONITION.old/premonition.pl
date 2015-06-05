#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Algorithm::Combinatorics qw(combinations) ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($mlength,$p1,$p2,$infile,$outfile,$dist,$which_tech,$listfile,$protein);
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
            "mlength=i"=>\$mlength ,
            "moveZ=i"=>\$moveZ ,
            "dist=f"=>\$dist ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $ofhclose = util_write("log.close");
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a dist -option -dist  ") if(!defined $dist);
usage( "Need to give a mlength -option -mlength  ") if(!defined $mlength);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my $file1 = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);

my @res = $pdb1->GetResidues();
my $N = @res;
my $cnt = 0 ;
my $sum = 0 ;
my $cntmatch = 0 ; 

my $DB = {};
my $done = {};
foreach my $res (@res){
	my $resnum = $res->GetResNum();
	system ("$SRC/ALIGN/premonitionsingleresidues.pl -outf $p1.out -p1 $p1 -dis $dist -ml $mlength -resnum $resnum ");
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

