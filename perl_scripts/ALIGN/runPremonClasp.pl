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
my $errfh = util_append($errlog);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a anndir -option -anndir  ") if(!defined $anndir);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a protein pdb id -option -premon  ") if(!defined $premon);
usage( "Need to give a protein pdb id -option -size  ") if(!defined $size);
usage( "Need to give a protein pdb id -option -tag  ") if(!defined $tag);
usage( "Need to give a protein pdb id -option -rundir  ") if(!defined $rundir);
usage( "Need to give a protein pdb id -option -readpotential  ") if(!defined $readpotential);
usage( "Need to give a protein pdb id -option -close2activesite  ") if(!defined $close2activesite);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;
    my @proteins ; 
    push @proteins, $protein ; 
	my $i = $protein;
    #next if( ! -e "$APBSDIR/$i/$i.pqr");
    #next  if( ! -e "$APBSDIR/$i/pot1.dx.atompot" && ! -e "$APBSDIR/$i/pot1.dx.atompot");
    my @info = util_ReadPdbs($PDBDIR,$APBSDIR,$readpotential,@proteins) ; 
    my $info = shift @info ;
    my $pdb1 = $info->{PDBOBJ};
    my $pqr1 = $info->{PQR};
    my $pots1 = $info->{POTS};
	my $numberofchains = $pdb1->ReadPDBAndSplit("$PDBDIR/$protein.pdb",$protein,"tmpfile");
	die "Please ensure each PDB has one chain" if($numberofchains ne 1);
	 



my $closeresidues  = {};
if($close2activesite){
	my $file = "$ANNDIR/$protein.$size.$tag.premon.in";
	if(-e $file){
         my ($status,$premoninInfo)  = util_ParsePremonIn($file,$errfh);
		 my @l = @{$premoninInfo->{REALSTRING}};
		 my @atoms = $pdb1->CreateAtomListFromString(@l);
		 foreach my $a (@atoms){
		 	 my ($results,$combined) = $pdb1->GetNeighbourHoodAtom(\@atoms,$DISTFORCLOSEATOMS);
             foreach my $j (@{$combined}){
		         my $resnum = $j->GetResNum(); 
				 $closeresidues->{$resnum} =1 ;
			  }
		  }
		 my $NNN = keys %{$closeresidues} ;
		 print "Found $NNN close atoms : @l \n";
	}
	else{
		print "Warning: Did not find close2 file $file\n";
		$close2activesite = 0 ; 
	}
}


if(-z $premon){
   print $errfh "$premon size is 0 - please check\n"; 
   die ;
}
$info = util_ReadPremonOut($premon,$size);



my $mapfileforPBDsHolo2Apo = util_read("$rundir/$anndir/mapfileHolo2Apo");
my $mapfileforPBDsHolo2ApoInfo = {};
while(<$mapfileforPBDsHolo2Apo>){
     next if(/^\s*$/);
     chop ;
	 my ($nm,$junk) = split ; 
	 $mapfileforPBDsHolo2ApoInfo->{$nm} = $junk ;
	 #print "$nm $junk\n";
}


my $fhfinallist = util_write("finallist");

my @list= util_read_list_sentences($listfile);
#### Process each reference 
my $cumulativeScores  = {};
foreach my $i (@list){
	my $claspoutdist = "$protein.$i.dist.out";
	my $claspoutfinal = "$protein.$i.pdb.out";
	my $file = "$rundir/$anndir/$i.$size.$tag.premon.in";
	if(! -e $file || -z $file){
		print $errfh "No file $file\n";
		next ;
	}

    my $pdb2 = new PDB();

	my $file2 ;
	if(exists $mapfileforPBDsHolo2ApoInfo->{$i}){
		my $X = $mapfileforPBDsHolo2ApoInfo->{$i} ;
		die "how come $i is not defiend" if(! defined $X);
        $file2 = "$PDBDIR/$X.pdb";
		print "Using $X instead of $i\n";
	}
	else{
        $file2 = "$PDBDIR/$i.pdb";
	}
    $pdb2->ReadPDB($file2);

	my @configs ; 
	my $justonce = 1 ;
	foreach my $IDX (1..4){
	     my $configtmp = "$rundir/$anndir/$i.$size.$tag.premon.in.config$IDX";
		 if(! -e $configtmp){
		 	 $configtmp = $config ;
		 }
		 else{
		 	 if($justonce){
		 	     print "Info:Using personalized configs $configtmp. Just reporting once.\n";
			     $justonce = 0 ;
			 }
		 }
		 push @configs, $configtmp;
	}

	print $fhfinallist "$i\n";
	next if(-e $claspoutfinal);
	#print "Processing $claspoutfinal \n" ;
	my $pymolin = "$protein.$i.pymol.in";
	my $logout = "$protein.$i.log";
	my $claspfhdist = util_write($claspoutdist); 
	my $claspfhfinal = util_write($claspoutfinal); 
	my $logfh = util_write($logout); 
	my $fhpymolin = util_write($pymolin); 

    my ($status,$premoninInfo)  = util_ParsePremonIn($file,$errfh);
	my @RESULTSTYLE = @{$premoninInfo->{RESULTSTYLE}};

	die "No resultstyle. Please rerun config generator" if(!@RESULTSTYLE);

	my @ALIST  ; 
	{
	   foreach my $i (0..3){
	       my $II = $RESULTSTYLE[$i];
	       my $a = $pdb2->ParseAtomLine($II);
	       push @ALIST, $a;
	   }
	}

	if(!$status){
		print $errfh "$file syntax not correct\n";
	}
	else{
		my $STRLIST = $premoninInfo->{STR} ;

		my @matches ;
		my $FINALSTR = "";
		foreach my $STR (@{$STRLIST}){
			$FINALSTR = $FINALSTR . "/" . $STR ;
		     if(!exists $info->{$STR}){
			     next ;
		     }
		     my $matchestmp = $info->{$STR};
		     push @matches,  @{$matchestmp} ;
		}

		if(!@matches){
			     #print $claspfhdist "# $FINALSTR does not exist in $premon\n";
			     #print $claspfhfinal "# $FINALSTR does not exist in $premon\n";
			     close($claspfhdist);
			     close($claspfhfinal);
				 next ;
		}




		my $D = $premoninInfo->{D} ;
		my $PD = $premoninInfo->{PD} ;

		my $sortedDistInfo = {};
		my $sortedPDInfo = {};
		my $sortedFinalScore = {};
		my $matchNames = {};
		my $matchPD = {};
		my $matchD = {};
		my $matchMaxDist = {};
		foreach my $match (@matches){

			$match =~ s/\./ /g;
			my (@numbers) = split " ", $match ;
			die if(@numbers ne $size);

			next if(RepeatedNumbers(@numbers));
			if($close2activesite){
				my $existsinclose = 0 ; 
				foreach my $n (@numbers){
				    if(exists $closeresidues->{$n}){
						$existsinclose = 1 ;
						last ;
					}
				}
				if(!$existsinclose){
					#print "Not in close so next\n";
					next ;
				}
			}
			print "Numbers are: @numbers \n" if($verbose);


			my @atoms ;
			my $nm = "";

			my $missingatom = 0 ;
			my $INDEX = 0 ;
			foreach my $number (@numbers){
	             my ($res) = $pdb1->GetResidueIdx($number);
				 my $configtmp = $configs[$INDEX];
				 $INDEX++;
                 ConfigPDB_Init($configtmp);
	             my $type = ConfigPDB_GetAtom($res->GetName()) or die;
				 #print "$type $number\n";
	             my ($atom1) = $pdb1->GetAtomFromResidueAndType($number,$type) ;
				 if(!defined $atom1){
				 	print "atom missing for $number and $type\n" if($verbose);
				 	$missingatom = 1 ;
				 }
				 $nm  = $nm . $res->GetName(). "/".  $number . "/$type" .  " ";
	             push @atoms, $atom1 ;
			}
			next if($missingatom);

			$matchNames->{$match} = $nm;
            my @distlist = @{$pdb1->DistanceInGivenSetOfAtoms(\@atoms)};
            my @pdlist  ; 
			if($readpotential){
			    @pdlist = @{$pdb1->PDInGivenSetOfAtoms(\@atoms,$pqr1,$pots1)};
            }

			my ($distscore,$maxdist) = util_ScoreDistance($D,\@distlist,$MAXALLOWEDDISTDEV);
			my ($pdscore) = 0 ;
			if($readpotential){
			  $pdscore = util_ScorePD($PD,\@pdlist);
			}

			$matchD->{$match} = \@distlist;
			$matchPD->{$match} = \@pdlist;
			$matchMaxDist->{$match} = $maxdist ;
			$sortedDistInfo->{$match} = $distscore ;
			$sortedPDInfo->{$match} = $pdscore ;
			$sortedFinalScore->{$match} = $distscore + $pdscore ;
		}


		## this can happen only if close2 is specified
		my $NNN = keys %{$sortedDistInfo} ;
		my ($junk,$head) = ($file =~ /(.*\/)(.*)/);
		if($verbose){
		    print "$head with $NNN for $FINALSTR\n";
		}
		else{
			print "...\n";
		}
		if(!$NNN){
		    die "Dont expect zero if close2activesite is $close2activesite" if(!$close2activesite);
			print $claspfhdist "# nothing close to active site\n";
			print $claspfhfinal "# nothing close to active site\n";
			close($claspfhdist);
			close($claspfhfinal);
			next ;
		}


		### sort first based on distance 
		my @sortMatchesDist =  sort { $sortedDistInfo->{$a} <=> $sortedDistInfo->{$b} } (keys %{$sortedDistInfo});
		my $cnt = 0  ;
		#print $logfh " ================= DISTANCE =====================\n";
		foreach my $match (@sortMatchesDist){
			my $distscore= $sortedDistInfo->{$match};
			my $nm = $matchNames->{$match};
			my $pdlist = $matchPD->{$match};
			my $dlist = $matchD->{$match};
			my $maxdist = $matchMaxDist->{$match};

			
			my $pdscore = 0 ; 
			my $total = $pdscore + $distscore ;
			print $claspfhdist "$nm @ $distscore $pdscore\t\t$total \n";


			if(0){
			print $logfh "\n\n$nm @ $distscore $pdscore\t\t$total \n";
			print $logfh  "MaxDist = $maxdist\n";
			print $logfh  "@{$D} \n";
			print $logfh  "@{$dlist} \n";
			print $logfh  "@{$PD} \n";
			print $logfh  "@{$pdlist} \n";
			}
			$cnt++;
			last if($cnt > $MAXMATCHES);
		}

		print $logfh " ================= DISTANCE and PD =====================\n";
		my @sortMatchesFinal =  sort { $sortedFinalScore->{$a} <=> $sortedFinalScore->{$b} } (keys %{$sortedFinalScore});
		$cnt = 0  ;
		foreach my $match (@sortMatchesFinal){
			my $distscore= $sortedDistInfo->{$match};
			my $nm = $matchNames->{$match};
			my $pdlist = $matchPD->{$match};
			my $dlist = $matchD->{$match};
			my $maxdist = $matchMaxDist->{$match};

			

			my $DDDD = 1000 ;
			if($cnt < 10 ){
			    my @l = split " ", $nm ; 
			    my @BLIST  ; 
			    foreach my $i (0..3){
				    my $II = $l[$i];
			        my $a = $pdb1->ParseAtomLine($II);
				    push @BLIST, $a;
			    }
                $DDDD = util_AlignAndMatchRemainingAtoms($pdb1,$pdb2,\@BLIST,\@ALIST,10);
			}

			my $pdscore = $sortedPDInfo->{$match} ; 
			my $total = $sortedFinalScore->{$match};
			print $claspfhfinal "$nm @ $distscore $pdscore\t\t$total\t\t$DDDD \n";

			if($cnt eq 0){
				print $fhpymolin "$nm \n";
				print $fhpymolin "@RESULTSTYLE \n";
			}

			print $logfh "\n\n$nm  @ $distscore $pdscore\t\t$total\n";
			print $logfh  "MaxDist = $maxdist\n";
			print $logfh  "@{$D} \n";
			print $logfh  "@{$dlist} \n";
			print $logfh  "@{$PD} \n";
			print "@{$PD} \n";
			print $logfh  "@{$pdlist} \n";
			$cnt++;
			last if($cnt > $MAXMATCHES);
		}


		#exit ; # if just one
	}
	close($claspfhdist);
	close($claspfhfinal);
	close($logfh);
	close($fhpymolin);
}
print "\n";



sub RepeatedNumbers{
	my (@n) = @_ ;
	my $done  = {};
	foreach my $i (@n){
		return 1 if(exists $done->{$i});
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
