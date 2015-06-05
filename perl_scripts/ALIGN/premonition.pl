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


my $EXTENDRADIUSBY = 2 ; 
my $DEFAULTRADIUS = 10 ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($mlength,$p1,$p2,$infile,$outfile,$dist,$which_tech,$listfile,$protein);
my (@expressions);
my $maxdist = 3 ;
my $verbose = 1 ;
my ($moveZ,$verify,$decaaf,$before1,$before2);
GetOptions(
            "protein=s"=>\$protein ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "mlength=i"=>\$mlength ,
            "dist=f"=>\$dist ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#my $ofh = util_write($outfile);
system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a mlength -option -mlength  ") if(!defined $mlength);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
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

my $info = {};
my $ifh = util_read($listfile);
while(<$ifh>){ 
     next if(/^\s*$/);
     chop ;
     my ($nm,$v) = split ;
     $info->{$nm} = $v + $EXTENDRADIUSBY; 
}
$info->{$p1} = $DEFAULTRADIUS if(! exists $info->{$p1});


my @files = ();
#print $ofh "$SRC/ALIGN/premonitionsingleresidues.pl -outf <out> -p1 $p1 -ml $mlength -resnum <all> \n";
foreach my $res (@res){
	next  if($res->GetAtomStr() ne "ATOM");
	my $resnum = $res->GetResNum();

	if(0){

    my $CAatom = $pdb1->GetAtomFromResidueAndType($resnum,"CA");

	#### Find residues in a distance $dist from the CA atom
	    my @atomlist ;
	    push @atomlist, $CAatom ;
	    my ($results,$combined) = $pdb1->GetNeighbourHoodAtom(\@atomlist,11);
	    my $sort ;
		my @goodatoms ; 
		my $residues = {};
		my @Rlist ;
        foreach my $j (@{$combined}){
		    my $atomstr = $j->GetAtomStr();
			next if($atomstr eq "HETATM");
		    my $NNN = $j->GetResNum(); 
			 $residues->{$NNN} = 1 ;
		}

		my $DDD = {};
        foreach my $j (keys %{$residues}){
		   my ($RRR) = $pdb1->GetResidueIdx($j);
		   if(! exists $DDD->{$RRR}){
			    $DDD->{$RRR} = 1 ;
			    push @Rlist, $RRR ; 
			}
			else{
				die "WJHHHHHHHY?";
			}
		}

		#push @Rlist, $res ;
		my $N = @Rlist ; 

		print $ofhlog "$p1 $resnum has $N residues, getting combinations of $mlength \n";
		next ;
	}

	my $file = "$p1.$resnum.singleout";
	push @files, $file ;
	if( ! -e $file){
	    $dist = $info->{$p1} ; 
	    print "$SRC/ALIGN/premonitionsingleresidues.pl -outf $file -p1 $p1 -dist $dist -ml $mlength -resnum $resnum\n";
	    system ("perl $SRC/ALIGN/premonitionsingleresidues.pl -outf $file -p1 $p1 -dist $dist -ml $mlength -resnum $resnum ");

	}
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

