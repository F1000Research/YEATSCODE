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
use Algorithm::Combinatorics qw(combinations permutations variations) ;




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ann,$config,$p1,$p2,$in1,$in2,$infile,$outfile,$which_tech,$listfile,$protein);
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
            "in1=s"=>\$in1 ,
            "in2=s"=>\$in2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a protein 2 id -option -p2  ") if(!defined $p2);
usage( "Need to give a input 1 id -option -in1  ") if(!defined $in1);
usage( "Need to give a input 2 id -option -in2  ") if(!defined $in2);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my $pdb = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

$pdb = "$PDBDIR/$p2.pdb";
my $pdb2 = new PDB();
$pdb2->ReadPDB($pdb);

my $ofh = util_write($outfile);

ConfigPDB_Init($config);
$, = " , ";

my @atoms1 = GetAtoms($pdb1,$in1);
my @atoms2 = GetAtoms($pdb2,$in2);
my $N1 = @atoms1 ;
my $N2 = @atoms2 ;
die if($N1 ne $N2);
my @distlist2 = @{$pdb2->DistanceInGivenSetOfAtoms(\@atoms2)};

my $iter1 = variations(\@atoms1,$N1);
my $bestmatches = {};
while (my $c1 = $iter1->next) {
        my @combo1 = @{$c1} ; 
        my @distlist1 = @{$pdb1->DistanceInGivenSetOfAtoms(\@combo1)};
        my $diff = util_PairwiseDiff(\@distlist1,\@distlist2);
        print @distlist1, "\n"; ;
        print " J $N1 $N2 $diff\n";
		$bestmatches->{$diff} = \@combo1;
}

foreach my $k (sort {$a <=> $b} keys %{$bestmatches}){
	print "Best = $k\n";
	my $v = $bestmatches->{$k};
    my @distlist1 = @{$pdb1->DistanceInGivenSetOfAtoms(\@{$v})};
    my @distlist2 = @{$pdb2->DistanceInGivenSetOfAtoms(\@atoms2)};
	foreach my $a (@{$v}){
		$a->Print();
	}
	print "Other\n";
	foreach my $a (@atoms2){
		$a->Print();
	}
	print @distlist1, "\n"; ;
	print @distlist2, "\n"; ;

	my ($d1,$pot1) = GetPairWisePotDiff($p1,\@{$v});
	my ($d2,$pot2) = GetPairWisePotDiff($p2,\@atoms2);
	my @pot1 = @{$pot1};
	my @pot2 = @{$pot2};
	print @pot1, "\n"; ;
	print @pot2, "\n"; ;
	# just one
	last ;
}

#######################################


sub GetAtoms{
	my ($PDB,$infile) = @_ ; 
    my @list= util_read_list_sentences($infile);
    my $list = {};
    my @reallist ; 
    map { s/,/ /g ; my @j = split ; push @reallist, @j;  } @list ;
    map { s/\s*//g ; $list->{$_} = 1 ; } @reallist ;

    my @atoms ; 
    foreach my $i (@reallist){
	    $i =~ s/,//g;
	    my ($name,$number) = ($i =~ /([a-zA-Z]+)([0-9]+)/);
	    $name = uc($name);
	    my $len = length($name); 
	    die "Wrong length $len" if($len != 1 && $len != 3);
	    if($len == 1){
	          $name = $PDB->GetThreeLetter($name);
	    }
	    my ($res) = $PDB->GetResidueIdx($number);
	    my $type = ConfigPDB_GetAtom($res->GetName()) or die;
	    my ($atom1) = $PDB->GetAtomFromResidueAndType($number,$type);
	    push @atoms, $atom1 ;
    }
	return @atoms ;
}


#my ($dist1,$pot1) = GetPairWise($p1,$in1);
#my ($dist2,$pot2) = GetPairWise($p2,$in2);
#
sub GetPairWisePotDiff{
   my ($protein,$atomlist) = @_ ;
   my $info = {};
   my @resultlines ;
   my @proteins ;
   push @proteins, $protein; 
   
      
   my @info = util_ReadPdbs($PDBDIR,$APBSDIR,1,@proteins) ; 
   my $info = shift @info ;
   my $pdb1 = $info->{PDBOBJ};
   my $pqr1 = $info->{PQR};
   my $pots1 = $info->{POTS};
      
   $, = " ";
   my ($dist,$pots,$name) = util_ProcessSingleLineAtomlist($pdb1,$pqr1,$pots1,$atomlist);
   my @dist = @{$dist}; 
   my @pots = @{$pots}; 
   
   return (\@dist,\@pots);
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
