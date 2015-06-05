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

my $MYCNT = 2 ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($mlength,$p1,$p2,$infile,$outfile,$dist,$which_tech,$listfile,$protein);
my (@expressions);
my $maxdist = 3 ;
my $verbose = 1 ;
my ($moveZ,$verify,$decaaf,$before1,$resnum,$before2);
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
            "resnum=i"=>\$resnum ,
            "mlength=i"=>\$mlength ,
            "moveZ=i"=>\$moveZ ,
            "dist=f"=>\$dist ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $ofhlog = util_append("logclose");
usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
usage( "Need to give a mlength -option -mlength  ") if(!defined $mlength);
usage( "Need to give a resnum -option -resnum  ") if(!defined $resnum);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my $alphatable = {};
my @ALPHA = qw( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z );
my $acnt = 0 ; 
foreach my $a (@ALPHA){
	$alphatable->{$a} = $acnt++;
}

my $file1 = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);

my @res = $pdb1->GetResidues();
my $N = @res;
my $sum = 0 ;

my $DB = {};
my $done = {};
my ($res) = $pdb1->GetResidueIdx($resnum);
{
    next  if($res->GetAtomStr() ne "ATOM");

    my $CAatom = $pdb1->GetAtomFromResidueAndType($resnum,"CA");

	#### Find residues in a distance $dist from the CA atom
	    my @atomlist ;
	    push @atomlist, $CAatom ;
	    my ($results,$combined) = $pdb1->GetNeighbourHoodAtom(\@atomlist,$dist);
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
			if(RepeatedNumbers(@Rlist)){
				$, = "\n";
				print @Rlist , "\n"; ;
				die "repeat" ;
			}

		print $ofhlog "$p1 $resnum has $N residues, getting combinations of $mlength \n";

		my $iter = combinations(\@Rlist, $mlength);
		my $CNT = 0 ;
		while (my $c = $iter->next) {

		    my @combo = @{$c} ;

			## sort based on single letter
		    my @resultssorted = sort { $a->SingleLetterValue($pdb1,$alphatable) <=> $b->SingleLetterValue($pdb1,$alphatable) } @{$c};
			if(0){
			   foreach my $r ( @{$c}){
				   my $nm = $r->GetName();
				   my $num = $r->GetResNum();
				   print "$nm$num, ";
			   }
				   print " \n";
			   foreach my $r ( @resultssorted){
				   my $nm = $r->GetName();
				   my $num = $r->GetResNum();
				   print "$nm$num, ";
			   }
				   print " \n";
   
			   die if($MYCNT > 3);
			   $MYCNT++;
			}
			
			

			## create the index string for the amino acid string
		    my $nm = "";
		    my $index = "";
			
			my $DONE = {};
			my $nextkaro = 0 ;
			foreach my $r (@resultssorted){
				my $n = $r->GetResNum();
				if(exists $DONE->{$n}){
					die "$n exist" ;
					$nextkaro = 1 ;
				}

				$DONE->{$n} = 1 ;

				my $type = $r->PrintSingleLetter($pdb1);
				$nm = $nm . $type ;
				$index = $index . "." . $n ; 
			}
			#print "NM = $nm\n";
			next if($nextkaro);

			## insert in the DB for that particular amino acid string
			next if(exists $done->{$index});
			$CNT++;
			if(!exists $DB->{$nm}){
				    $DB->{$nm} = [];	
		    }
			$done->{$index} = 1 ;
		    push @{$DB->{$nm}}, $index ; 
		}
		print "$p1 $resnum has $N residues, getting combinations of $mlength - total = $CNT \n";
		print $ofhlog "$p1 $resnum has $N residues, getting combinations of $mlength - total = $CNT \n";
}

print $ofh "# mlength=$mlength dist=$dist \n";
foreach my $k (keys %{$DB}){
	my @l = @{$DB->{$k}};
    print $ofh "$k ";
	foreach my $i (@l){
		print $ofh "$i ";
	}
    print $ofh "\n";
}

sub RepeatedNumbers{
	my (@n) = @_ ;
	my $done  = {};
	foreach my $i (@n){
		if(exists $done->{$i}){
			print "RepeatedNumbers $i\n";
		return 1 ;
         } 
		$done->{$i} = 1 ;
	}
	return 0 ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

