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
my ($infile,$radii,$outfile,$which_tech,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "radii=f"=>\$radii ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

$outfile = "$protein.maxdist.out" if(!defined $outfile);
my $ofh = util_write($outfile);

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

 my @res = $pdb1->GetResidues();
 my @allatoms ;
foreach my $res (@res){
		my @atoms1 = $res->GetAtoms();
		push @allatoms, @atoms1;
}
print "Ignoring consecutive residues\n";
while(@allatoms){
		    my $a2 = shift @allatoms ;
		    foreach my $a1 (@allatoms){
		           my $d = $pdb1->DistanceAtoms($a1,$a2); 
				   my $nm1 = $a1->GetNameSlashSep();
				   my $nm2 = $a2->GetNameSlashSep();

				   my $N1 = $a1->GetResNum();
				   my $N2 = $a2->GetResNum();
				   my $DIFF = abs($N1-$N2);
				   #next if($N1 eq $N2 || $DIFF eq 1 );
				   next if($N1 eq $N2 );

				   my $nm = "$nm1.$nm2";
				   print $ofh "$protein $nm $d\n";
	        }

}

close($ofh);
system("sort.pl -in $outfile -out $outfile.sort -idx 2");



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
