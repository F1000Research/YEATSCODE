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


$, = "\t";
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($tag,$close2activesite,$rundir,$anndir,$infile,$outfile,$premon,$errlog,$readpotential,$listfile,$protein);
my ($config,@expressions);
my $size ;
my $verbose = 0 ;
my $MAXMATCHES = 10000 ;
my $MAXALLOWEDDISTDEV = 3 ;
my $DISTFORCLOSEATOMS = 5 ;
GetOptions(
            "premon=s"=>\$premon ,
            "rundir=s"=>\$rundir ,
            "protein=s"=>\$protein ,
            "config=s"=>\$config,
            "anndir=s"=>\$anndir,
            "infile=s"=>\$infile ,
            "tag=s"=>\$tag ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "errlog=s"=>\$errlog ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
            "close2activesite=i"=>\$close2activesite ,
            "readpotential=i"=>\$readpotential ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a protein -option -protein  ") if(!defined $protein);
usage( "Need to give a protein pdb id -option -outfile  ") if(!defined $outfile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();


my $ofh = util_append($outfile);
my $PWD = cwd;

    my $pdb = "$PDBDIR/$protein.pdb";
    my $pdb1 = new PDB();
    $pdb1->ReadPDB($pdb);
	print "reading\n";
	my $numberofchains = $pdb1->ReadPDBAndSplit("$PDBDIR/$protein.pdb",$protein,"tmpfile");
	print $ofh "$protein\n" if($numberofchains ne 1);
	print "$protein\n" if($numberofchains ne 1);
	 


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
