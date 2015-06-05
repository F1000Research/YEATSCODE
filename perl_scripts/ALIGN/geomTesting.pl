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
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $ofhclose = util_write("log.close");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my $file1 = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);

if(defined $verify){
	$before1 = $pdb1->DistanceMatrix();
}


my $info = {};





#print "Parsing results file\n";
my ($atoms1,$atoms2) = pymolin_getResultsLine($infile,$pdb1,$pdb1);
my @allatoms = $pdb1->GetAtoms();

my @atoms1 = @{$atoms1};
my $a0 = $atoms1[0];
my $a1 = $atoms1[1];
my $a2 = $atoms1[2];
#my ($newX,$newY,$newZ) = $pdb1->AlignXto2Atoms($a0,$a1);
my ($newX,$newY,$newZ) = $pdb1->AlignXto2AtomsInSet($a0,$a1,\@allatoms);

$pdb1->WritePDB("AAA",1);

my @l = geom_GetSetOfVectorsForFullCircle(180,10);
while(@l){
	my $y = shift @l ;
	my $z = shift @l ;
    my $p2 = vector( 0, $y, $z);
    my ($XXX,$newY1,$ZZZ,$vec,$R) = geom_RotateWithXConstantAndZwillbeZero($p2);
    $pdb1->ApplyRotationMatrix($R,\@allatoms);
	last ;
}

$pdb1->WritePDB("BBB",1);

foreach my $a (@atoms1){
	$a->Print();
}

if(defined $verify){
	my $after1 = $pdb1->DistanceMatrix();
	$pdb1->VerifyDistanceMatices($after1,$before1);
}

$pdb1->WritePDB($outfile,1);

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

