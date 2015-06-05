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

use AAConfig;

my $aaconfig = new AAConfig("/home/sandeepc/aa.config");



my $ADDONLYREACTIVE = 0 ; 
my $MATCHREVERSE = 0 ;
my $POTENTIALMATCH = 0 ;
my $DISTANCEMATCH = 1  ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ignoremissingresidues,$ignoreGLY,$debugres,$doenergy,$ann,$config,$p1,$p2,$infile,$threshPD,$threshsign,$threshDist,$outfile,$readpotential,$which_tech,$listfile,$protein);
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 0 ;

$threshPD = 150 ;
#$threshsign = 150 ;
$threshsign = 100 ;
$threshDist = 3 ;
$threshDist = 2 ;
my ($tag,$onlypolar,$radii,$before1,$before2);
$readpotential = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "onlypolar=i"=>\$onlypolar ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "ignoregly"=>\$ignoreGLY ,
            "listfile=s"=>\$listfile ,
            "tag=s"=>\$tag ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
            "threshPD=i"=>\$threshPD ,
            "threshsign=i"=>\$threshsign ,
            "readpotential=i"=>\$readpotential ,
            "debugres=s"=>\$debugres ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a radii file name => option -radii ") if(!defined $radii);
usage( "Need to give a tag file name => option -tag ") if(!defined $tag);
usage( "Need to give a tag file name => option -tag ") if(!defined $tag);
usage( "Need to give a onlypolar file name => option -onlypolar ") if(!defined $onlypolar);


my $debugreslist = {};
my $debugresfh ;
if(defined $debugres){
    my @list= util_read_list_sentences($debugres);
    map { s/\s*//g ; $debugreslist->{$_} = 1 ; } @list ;
	$debugresfh = util_write($debugres . ".out");
}

my $ofh = util_write($outfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

ConfigPDB_Init($config,$ofh);

my @proteins ;
push @proteins, $p1 ;
push @proteins, $p2 ;

my @info = util_ReadPdbs($PDBDIR,$APBSDIR,$readpotential,@proteins);
my $pdb1 = $info[0]->{PDBOBJ};
my $pdb2 = $info[1]->{PDBOBJ};
my $pqr1 = $info[0]->{PQR};
my $pqr2 = $info[1]->{PQR};
my $pots1 = $info[0]->{POTS};
my $pots2 = $info[1]->{POTS};


my $idxtable1 = ProcessOnePDB($pdb1,$pqr1,$pots1,$p1);
my $idxtable2 = ProcessOnePDB($pdb2,$pqr2,$pots2,$p2);


my $fhmeanPD = util_write("meanPD.$p1.$p2.$radii");
my $fhmeanPDtag = util_write("meanPD.$p1.$p2.$radii.$tag");
my $fhmeanPDabs = util_write("meanabsPD.$p1.$p2");
my $fhmeanabsPDtag = util_write("meanabsPD.$p1.$p2.$radii.$tag");
my $fhsdPD = util_write("sdPD");
my $fhsdPDabs = util_write("sdabsPD");

my $fhmeanDD = util_write("meanDD.$p1.$p2.$radii");
my $fhmeanDDtag = util_write("meanDD.$p1.$p2.$radii.$tag");
my $fhmeanDDabs = util_write("meanabsDD.$p1.$p2.$radii");
my $fhmeanabsDDtag = util_write("meanabsDD.$p1.$p2.$radii.$tag");
my $fhsdDD = util_write("sdDD");
my $fhsdDDabs = util_write("sdabsdD");

my $fhthreshPD = util_write("threshPD.$radii");
my $fhthreshDD = util_write("threshDD.$radii");

my $fhchangesign = util_write("changesign.$p1.$p2");
my $fhchangesigngraph = util_write("changesign.$p1.$p2.$radii");
my $fhvector = util_write("vector.$radii");
my $fhenergy = util_write("energeydiff.$radii");
my $CNT = 0 ;
my $doneRes = {};


# We have the info for both proteins in idxtable - each residue 
# has a table that gives the distance/potential  for every other atom
my @greaterthan;
my @lesserthan;
my @postivePD;
my @negativePD;
foreach my $i (sort {$a <=> $b}  keys %{$idxtable1}){

	my $p1 = $idxtable1->{$i}->{PDIFF} ;
	my $d1 = $idxtable1->{$i}->{DDIFF} ;
	my $resname1 = $idxtable1->{$i}->{RESNAME} ;
	my $totE1 =	$idxtable1->{$i}->{TOTALENERGY} ;

	my $p2 = $idxtable2->{$i}->{PDIFF} ;
	my $d2 = $idxtable2->{$i}->{DDIFF} ;
	my $totE2 =	$idxtable2->{$i}->{TOTALENERGY} ;
	

	next if(!defined $p2 || !defined $d2);
	
	my $energeydiff = $totE1 - $totE2 ;
	print $fhenergy "$i $energeydiff\n";


	#my $fnamep = "$i.PD";
    #my $fhp = util_write($fnamep);
    my @allPD;
    my @allPDabs ;
    my @allDD;
    my @allDDabs ;

	my @processedAtoms ; 
    foreach my $j (sort {$a <=> $b}  keys %{$idxtable2}){
		if($i eq $j){
		    $CNT = $CNT + CheckSanity($pdb1,$pdb2,$i);
			die "Too many mismatches $CNT" if($CNT > 10);
			next ;
		}

		my $str = "$i.$j";


		$ignoremissingresidues = 1 ;
		if(defined $ignoremissingresidues){
		if(!exists $p1->{$str} || !exists $p2->{$str} || !exists $d1->{$str} || !exists $d2->{$str}){
			 print "Not there for $str \n" if($verbose);
			 next ;
		}
		}
		else{
		die "$str" if(!exists $p1->{$str});
		die "$str not in p2" if(!exists $p2->{$str});
		die "$str" if(!exists $d1->{$str});
		die "$str" if(!exists $d2->{$str});
		}

		my $DD = util_format_float($d1->{$str} - $d2->{$str}, 1);
		my $DDstr = "$d1->{$str} - $d2->{$str}";
		

		## ignore atoms beyond a certain distance
		if(abs($d1->{$str}) > $radii){
			next ;
		}

		my $PD = util_format_float($p1->{$str} - $p2->{$str}, 1);
		my $PDstr = "$p1->{$str} - $p2->{$str}";


		my $r = $pdb1->GetResidueIdx($j);
		my $resname2 = $r->GetName();
		my $aIsPolar = $aaconfig->IsPolar($resname1) ;
		my $bIsPolar = $aaconfig->IsPolar($resname2) ;
		my $IsPolar = $aIsPolar && $bIsPolar ;


		next if($onlypolar eq 1 && !$IsPolar);

		push @processedAtoms, $j ;

		my $abs= abs($PD);
		push @allPD, $PD;
		push @allPDabs , $abs;

		my $absDD = abs($DD);
		push @allDD, $DD;
		push @allDDabs, $absDD;

		if(defined $debugres && exists $debugreslist->{$i}){
		    my $A = $pdb1->GetResidueIdx($i);
		    my $B = $pdb1->GetResidueIdx($j);
			my $name1 = $A->GetName();
			my $name2 = $B->GetName();
	        my $typeA = ConfigPDB_GetAtom($A->GetName());
	        my $typeB = ConfigPDB_GetAtom($B->GetName());
			print $debugresfh "KKKKK$name1/$i/$typeA $name2/$j/$typeB $DDstr diff= $DD $PDstr diff= $PD \n";
		}


		my $changesign = 0 ; 
		$changesign = 1 if($p1->{$str} > 0 && $p2->{$str} < 0);
		$changesign = 1 if($p2->{$str} > 0 && $p1->{$str} < 0);

		# for a single residue, map the potential deviations
		#print $fhp "$j $PD";


		## dont do abs here - else we will repeat for other direction
		if($changesign && $PD > $threshsign){
			 if($aIsPolar && $bIsPolar && !exists $doneRes->{$i}){
		         my $A = $pdb1->GetResidueIdx($i);
		         my $B = $pdb1->GetResidueIdx($j);
	             my $typeA = ConfigPDB_GetAtom($A->GetName());
	             my $typeB = ConfigPDB_GetAtom($B->GetName());
	             #$typeA = "CA";
	             #$typeB = "CA";
	             my $atom1 = $pdb1->GetAtomFromResidueAndType($i,$typeA);
	             my $atom2 = $pdb1->GetAtomFromResidueAndType($j,$typeB);
				 my $DIST1 = $pdb1->DistanceAtoms($atom1,$atom2);
	             $atom1 = $pdb2->GetAtomFromResidueAndType($i,$typeA);
	             $atom2 = $pdb2->GetAtomFromResidueAndType($j,$typeB);
				 my $DIST2 = $pdb2->DistanceAtoms($atom1,$atom2);
		         print $fhchangesign "$resname1/$i/$typeA $resname2/$j/$typeB $p1->{$str} $p2->{$str} $PD $DIST1 $DIST2  \n";
		         print $fhchangesigngraph "$i $j \n";
				 $doneRes->{$i} = 1 ;
			 }
		}
	}


	#close($fhp);

	if(@allPD){

	my ($final,$l1,$l2) = GetVectorForSingleAtom($i,\@processedAtoms,$pdb1,$pqr1,$pots1);
	my $length = $final->length ;
	print $fhvector "$i $length\n";

	    my ($meanPD,$PDabs) = PrintData($i,\@allPD,\@allPDabs,$fhmeanPD,$fhmeanPDtag,$fhmeanPDabs,$fhmeanabsPDtag,$fhsdPD,$fhsdPDabs);
	    my ($meanDD,$DDabs) = PrintData($i,\@allDD,\@allDDabs,$fhmeanDD,$fhmeanDDtag,$fhmeanDDabs,$fhmeanabsDDtag,$fhsdDD,$fhsdDDabs);
    
		my $r = $pdb1->GetResidueIdx($i);
		my $nm = $r->GetName();
        if(abs($meanPD) > $threshPD){
		    print $fhthreshPD "$nm$i $meanPD \n";
			push @greaterthan, $i ; 
			if($meanPD > 0){
			     push @postivePD, $i ; 
			}
			else{
			     push @negativePD, $i ; 
			}
    
        }
		else{
			push @lesserthan, $i ; 
		}
        if(abs($DDabs) > $threshDist){
		    print $fhthreshDD "$nm$i $DDabs \n";
    
        }
	}
}

my $gtstr = join " ", @greaterthan ;
my $posstr = join " ", @postivePD ;
my $negstr = join " ", @negativePD ;
print "red 0 $gtstr \n";
print "red 0 $posstr \n";
print "blue 0 $negstr \n";


sub PrintData{

	my ($i,$all,$allabs,$fhmean,$fhmeantag,$fhmeanabs,$fhmeanabstag,$fhsd,$fhsdabs) = @_ ;
	my $mean = Math::NumberCruncher::Mean($all) or warn "Mean not found" ;
	my $sd = Math::NumberCruncher::StandardDeviation($all) or warn "sd not found" ;
	my $meanabs = Math::NumberCruncher::Mean($allabs) or warn "Mean not found" ;
	my $sdabs = Math::NumberCruncher::StandardDeviation($allabs) or warn "sd not found" ;

	$mean = util_format_float($mean,1);
	$sd = util_format_float($sd,1);
	$meanabs = util_format_float($meanabs,1);
	$sdabs = util_format_float($sdabs,1);

	print $fhmean "$i $mean\n";
	print $fhmeantag "$i $mean\n";
	print $fhmeanabs "$i $meanabs\n";
	print $fhmeanabstag "$i $meanabs\n";
	print $fhsd "$i $sd\n";
	print $fhsdabs "$i $sdabs\n";
	return ($mean,$meanabs) ;
}

sub ProcessOnePDB{
    my ($pdb,$pqr,$pots,$protein) = @_ ;
    my @reslist1 = $pdb->GetResidues();
    my $idxtable = {};
    foreach my $r (@reslist1){
	    my $nm = $r->GetName();
		next if($nm eq "GLY" && defined $ignoreGLY);
	    my $idx = $r->GetIdx();
	    my $atomstr = $r->GetAtomStr();
	    next if($atomstr eq "HETATM");
    
	    my $type = ConfigPDB_GetAtom($r->GetName());
	    die "undefined type for $nm" if(! defined $type);
	    my $reactiveatom = $pqr->GetAtomFromResidueAndType($idx,$type);
	    my $potatom  = util_GetPotForAtom($reactiveatom,$pqr,$pots) ;
	    #print "$nm $idx $potatom \n";
	    $idxtable->{$idx} = {};
	    #$idxtable->{$idx}->{ATOM} = $reactiveatom ;
	    $idxtable->{$idx}->{POT} = $potatom ;
	    $idxtable->{$idx}->{RESNAME} = $nm ;

        my $list = util_make_list($reactiveatom);
	    my ($junk,$neigh)  = $pqr->GetNeighbourHoodAtom($list,$radii);

	    my $done ;
	    my $sort = {} ;
		my $totE = 0  ; 
		my @l = @{$neigh};
		push @l , $reactiveatom ;
		if(defined $doenergy){
	    foreach my $a (@l){
			my $nm = $a->GetName();
			next if(exists $done->{$nm});
			$done->{$nm} = 1 ;
			if($a->GetCharge()){
	             my $pot  = util_GetPotForAtom($a,$pqr,$pots) ;
				 my $E = $pot*$a->GetCharge();
				 $totE = $E + $totE ;
			}
		}
		}
		$idxtable->{$idx}->{TOTALENERGY} = $totE ; 
    }




    my $warned = {};
	### find pairwise properties
    foreach my $i (sort {$a <=> $b}  keys %{$idxtable}){
	     #my $a1 = $idxtable->{$i}->{ATOM} ;
	     my $p1 = $idxtable->{$i}->{POT} ;
		 if(!defined $p1){
			 	print STDERR "Potential not defined for $i\n" if(! exists $warned->{$i});
				$warned->{$i} = 1 ;
				next ;
		 }



	     my ($p,$d,$NORM) ; 
        foreach my $j (sort {$a <=> $b}  keys %{$idxtable}){
	  	     next if($i eq $j);
	         #my $a2 = $idxtable->{$j}->{ATOM} ;

		         my $A = $pdb->GetResidueIdx($i);
		         my $B = $pdb->GetResidueIdx($j);
	             my $typeA = ConfigPDB_GetAtom($A->GetName());
	             my $typeB = ConfigPDB_GetAtom($B->GetName());
	             #$typeA = "CA";
	             #$typeB = "CA";
	             my $a1 = $pdb->GetAtomFromResidueAndType($i,$typeA);
	             my $a2 = $pdb->GetAtomFromResidueAndType($j,$typeB);
	         my $p2 = $idxtable->{$j}->{POT} ;
			 if(!defined $p2 ){
			 	print STDERR "Potential not defined for $j\n" if(! exists $warned->{$j});
				$warned->{$j} = 1 ;
				next ;
			 }
		     my $pdiff = $p1 - $p2 ; 
		     my $dist = $a1->Distance($a2);

		if(defined $debugres && $debugres eq $i){
		}


		     my $str = "$i.$j";
		     $p->{$str} =   util_format_float($pdiff,1) ;
		     $d->{$str} =   util_format_float($dist,1);

			 #print "AAAAA $str $dist $protein \n";
	    }
	    $idxtable->{$i}->{PDIFF} = $p ;
	    $idxtable->{$i}->{DDIFF} = $d ;
    }
    return $idxtable ;
}

sub CheckSanity{
	my ($P1,$P2,$i) = @_ ;
	my $r1 = $P1->GetResidueIdx($i);
	my $r2 = $P2->GetResidueIdx($i);
	my $n1 = $r1->GetName();
	my $n2 = $r2->GetName();
	print STDERR "Mutation:$n1 $n2 $i\n" if( $n1 ne $n2);
	return 0 if($r1->GetName() eq $r2->GetName());
	return 1 ; 
}


sub usage{
my ($msg) = @_ ;
print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
die ;
}

sub GetVectorForSingleAtom{
	my ($i,$ly,$pdb,$pqr,$pots) = @_ ;
	my $r = $pdb->GetResidueIdx($i);
	my $type = ConfigPDB_GetAtom($r->GetName());
	die "undefined type for " if(! defined $type);
	my $x = $pqr->GetAtomFromResidueAndType($i,$type);
	
	
	my $XV = $pdb->MakeVector($x);
	my $potx  = util_GetPotForAtom($x,$pqr,$pots) ;
	return if(! defined $potx) ;

	my @y = @{$ly};
	my @YVs ;
	my @Ypots ;
	foreach my $ii (@y){
	     my $r = $pdb->GetResidueIdx($ii);
	     my $type = ConfigPDB_GetAtom($r->GetName());
	     die "undefined type for " if(! defined $type);

	     my $y = $pqr->GetAtomFromResidueAndType($ii,$type);
		my $YV = $pdb->MakeVector($y);
	    my $poty  = util_GetPotForAtom($y,$pqr,$pots) ;
		next if(! defined $poty);

		push @YVs, $YV ;
		push @Ypots, $poty ;
	}
	my ($final,$l1,$l2) = GetVectorForSinglePoint($XV,\@YVs,$potx,\@Ypots);
	return ($final,$l1,$l2) ;
}


sub GetVectorForSinglePoint{
	my ($XV,$YVs,$potx,$lpoty) = @_ ; 
	my @YVs = @{$YVs} ;
	my @Ypots = @{$lpoty} ;
	
	my $N = @YVs - 1 ; 
	my @V ;
	my @Vabs ;
	my $final ; 
	foreach my $idx (0..$N){
		my $poty = $Ypots[$idx];
		my $YV = $YVs[$idx];
		my $V ; 
		my $diff = $potx - $poty ;
		if($diff > 0 ){
		   $V = $YV - $XV ;	
		}
		else{
		   $V = $XV - $YV ;	
		}
	    my $norm = $V->norm;
	    my $length = $V->length;

		my $Vabs = $norm * abs($diff) ;
		push @V, $norm ;
		push @Vabs, $Vabs ;

		if(!defined $final){
			$final = $Vabs ;
		}
		else{
			$final = $final + $Vabs ;
		}
	}



	return ($final,\@V,\@Vabs);
}
