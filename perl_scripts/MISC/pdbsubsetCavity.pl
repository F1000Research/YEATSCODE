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
my ($infile,$outfile,$start,$end,$which_tech,$listfile,$protein);
my ($castp,@expressions);
my $subtract = 0 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "castp=s"=>\$castp ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "subtract=i"=>\$subtract ,
            "start=i"=>\$start ,
            "end=i"=>\$end ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a castp pdb id -option -castp  ") if(!defined $castp);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);


# in case you use fpocket
#my ($infofpocket,$infoperatom,$volume) = parseFPocketFile($fpocket);
my ($infocastp,$infobiggest) = parseCastp($castp,$subtract);

my $tableofres = {};
foreach my $i (keys %{$infobiggest}){
    $tableofres->{$i} = 1 ;
}

$pdb1->ReadPDBAndWriteSubset($pdb,$tableofres,$ofh);

chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
