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
my ($infile,$radii,$p1,$p2,$outfile,$which_tech,$hetatm,$listfile,$protein);
my (@expressions);
my $cutoff = 4 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "radii=f"=>\$radii ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a protein 2 id -option -p2  ") if(!defined $p2);
my $CNT = 0 ; 
$outfile = "$p1.$p2.maxdist.out" if(!defined $outfile);
my $ofh = util_write($outfile);

my $file1 = "$PDBDIR/$p1.pdb";
my $file2 = "$PDBDIR/$p2.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);
my $pdb2 = new PDB();
$pdb2->ReadPDB($file2);

my @CA1 = GetCA($pdb1);
my @CA2 = GetCA($pdb2);

my @allatoms1 ;
my @allatoms2 ;
my $done1 = {};
my $done2 = {};
while(@CA1){
		    my $a1 = shift @CA1 ;
		    foreach my $a2 (@CA2){
		           my $d = $pdb1->DistanceAtoms($a1,$a2); 
				   if($d < 20){
				   #if(1){
				   	   my $N1 = $a1->GetResNum();
				   	   my $N2 = $a2->GetResNum();
					   my $r1 = $pdb1->GetResidueIdx($N1);
					   my $r2 = $pdb2->GetResidueIdx($N2);

					   if(! exists $done1->{$N1}){
		                   my @atoms1 = $r1->GetAtoms();
		                   push @allatoms1, @atoms1;
						   $done1->{$N1} = 1 ;
					   }
					   if(! exists $done2->{$N2}){
		                   my @atoms2 = $r2->GetAtoms();
		                   push @allatoms2, @atoms2;
						   $done2->{$N2} = 1 ;
					  }

				   }
	        }

}

my $A = {};
my $B = {};
my $AB = {};
while(@allatoms2){
		    my $a2 = shift @allatoms2 ;
		    foreach my $a1 (@allatoms1){
		           my $d = $pdb1->DistanceAtoms($a1,$a2); 
				   next if($d > $cutoff);
				   my $nm1 = $a1->GetNameSlashSep();
				   my $nm2 = $a2->GetNameSlashSep();

				   my $N1 = $a1->GetResNum();
				   my $N2 = $a2->GetResNum();

				   my $X = $a1->GetResName();
				   my $Y = $a2->GetResName();

				   $A->{$N1} = $X ;
				   $B->{$N2} = $Y ;
				   $AB->{$X.$N1.$Y.$N2} = 1 ;
				   #next if($N1 eq $N2);

				   my $nm = "$nm1.$nm2";
				   print $ofh "$p1 $p2 $nm $d\n";

	        }

}

print "For first: ";
foreach my $k (sort {$a <=> $b} keys %{$A}){
	my $v = $A->{$k};
	print "$v$k ";
}
print "\n";

print "For second: ";
foreach my $k (sort {$a <=> $b} keys %{$B}){
	my $v = $B->{$k};
	print "$v$k ";
}

print "\n";
print "For BOTH: ";
foreach my $k (keys %{$AB}){
	my $v = $AB->{$k};
	print "$k ";
}
print "\n";

close($ofh);
system("wc -l $outfile");
system("sort.pl -in $outfile -out $outfile.sort -idx 3");


sub GetAllAtoms{
	my ($P) = @_ ;
    my @res = $P->GetResidues();
    my @allatoms ;
   foreach my $res (@res){
		   my @atoms1 = $res->GetAtoms();
		   push @allatoms, @atoms1;
   }
   return @allatoms ;
}
sub GetCA{
	my ($P) = @_ ;
    my @res = $P->GetResidues();
    my @allatoms ;
   foreach my $res (@res){
		   #my @atoms1 = $res->GetAtoms();
		   my $a = $res->GetAtomType("CA");
		   push @allatoms, $a if(defined $a);
   }
   return @allatoms ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
