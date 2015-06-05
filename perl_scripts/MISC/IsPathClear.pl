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
use Algorithm::Combinatorics qw(combinations) ;




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ann,$config,$p1,$p2,$infile,$outfile,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;

my ($verify,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "config=s"=>\$config,
            "radii=f"=>\$radii ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a radii file name => option -radii ") if(!defined $radii);
usage( "Need to give a protein file name => option -protein ") if(!defined $protein);


my $ofh = util_write($outfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);
ConfigPDB_Init($config);


my $info = {};
my @resultlines ;
my $ifh = util_read($infile);
my @proteins ;
push @proteins, $protein; 


$, = " ";
while(<$ifh>){
     next if(/^\s*$/);
     next if(/RESULT/);
     chop ;

	 #my ($dist,$pots,$name) = util_ProcessSingleLine($pdb1,$pqr1,$pots1,$_);
	 my $atomlist = $pdb1->ParseResultLine($_);
	 foreach my $atom (@{$atomlist}){
	 	print "$atom \n";
	 }
     my $iter = combinations($atomlist, 2);
     while (my $c = $iter->next) {
        my @combo = @{$c} ; 
		my ($a,$b) = @combo ; 
		$a->Print();
		$b->Print();
		$pdb1->IsPathClearBetweenTwoAtoms($a,$b,$radii);
	 }

	 ### just do one: 
	 die ;
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
