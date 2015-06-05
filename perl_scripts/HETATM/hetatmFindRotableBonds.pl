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
my $cutoff = 2 ;
my $MAXDISTDB_O_C = 1.25;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "radii=f"=>\$radii ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=f"=>\$cutoff ,
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

my $N = @allatoms ;

my @ROT ; 
my @REACTIVE_O_S ; 
my @DB_O_C ; 

my $TOPINFO = {};
my $ROTATABLE = {};

foreach my $a2 (@allatoms){
			my $nm2 = $a2->GetNameSlashSep();
			my @l = split "/", $nm2 ;
			$nm2 = $l[2];
			#next if(!($nm2 =~ /N/));

			my $info = {};
		    foreach my $a1 (@allatoms){
				next if($a1 eq $a2);
		           my $d = $pdb1->DistanceAtoms($a1,$a2); 
				   my $nm1 = $a1->GetNameSlashSep();
			       my @l = split "/", $nm1 ;
			       $nm1 = $l[2];
				   $info->{$nm1} = $d ; 


				   my $nm = "$nm1.$nm2";
	        }
			my @sorted = sort { $info->{$a} <=> $info->{$b} } (keys %{$info});

			my $ekhogaya = 0 ;
			my $cnt = 0 ;

			my $connections = "";
			my $lastconnection ;
			my $lastdist ;
			my @connections ;
			foreach my $i (@sorted){
				my $d = $info->{$i};
				last if($ekhogaya && $d > $cutoff);
				push @connections, $i;
				$lastconnection = $i ;
				$lastdist = $d ;
				$ekhogaya = 1 ;
				$connections = $connections . " $i $d ";
				$cnt++;
			}
			$TOPINFO->{$nm2} = \@connections ;
			die "$nm2 has no connections" if($cnt eq 0);
			print "$nm2 $connections \n ";



			## rotable bonds
			if($cnt eq 2){
				## ignore Carbons that are part of a ring
			    if(!($nm2 =~ /C/)){
				    push @ROT, "ROT $nm2 $connections";
					$ROTATABLE->{$nm2} = 1 ;
				}
				
			}

			## O connected to sulphur
			## O double bond to carbon
			if($cnt eq 1){
			    if(($nm2 =~ /O/)){
			        if(($lastconnection =~ /S/)){
				        push @REACTIVE_O_S, "REACTIVE_O_S $nm2 $lastconnection $lastdist"
					}
			        elsif(($lastconnection =~ /C/)){
					    warn "$lastdist is larger than $MAXDISTDB_O_C for $nm2 and $lastconnection" if($lastdist > $MAXDISTDB_O_C);
				        push @DB_O_C, "DB_O_C $nm2 $lastconnection $lastdist";
					}
				}
				
			}

}

PrintArray(@ROT);
PrintArray(@DB_O_C);
PrintArray(@REACTIVE_O_S);

foreach my $rot (keys %{$ROTATABLE}){
	my @l = @{$TOPINFO->{$rot}} ;
	print "$rot @l  \n";
	my $one = $l[0];
	my $other = $l[1];

	my $oneside = GetAtomsFromRotable($rot,$one);
	my $otherside = GetAtomsFromRotable($rot,$other);
	print "One side ";
	foreach my $ATOM (keys %{$oneside}){
		print "$ATOM ";
	}
	print "\n";
	print "Other side ";
	foreach my $ATOM (keys %{$otherside}){
		print "$ATOM ";
	}
	print "\n";

}


sub GetAtomsFromRotable{
	my ($rot,$one) = @_ ;
	my $ONESIDE ;
	$ONESIDE->{$rot} = 1 ;
	$ONESIDE->{$one} = 1 ;

	my $done = {};
	$done->{$rot} = 1 ;

	while(1){
		my $added = 0 ;
		foreach my $ATOM (keys %{$ONESIDE}){
			next if(exists $done->{$ATOM});
			next if(exists $ROTATABLE->{$ATOM});

	        $done->{$ATOM} = 1 ;
	        my @conns = @{$TOPINFO->{$ATOM}} ;
		    foreach my $conn (@conns){
			     next if(exists $ROTATABLE->{$conn});
			     if(! exists $ONESIDE->{$conn}){
			         $ONESIDE->{$conn} = 1 ; 
					 $added = 1 ;
				 }
		     }
		}
		last if(!$added);
	}

	return $ONESIDE ;
}



sub CreateSet{
	my ($atom) = @_ ;

}


sub PrintArray{
	my (@l) = @_ ;
	foreach my $i (@l){
		print $ofh "$i\n";
	}
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
