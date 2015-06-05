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


my $distwater = 1.4 ; 
my $deltacube = 3 ;
my $cutoffview = 3 ;

my $RECALIBRATE = 0 ;

my $padding = 3 ; 
my $DISTSURFACE = 1 ; 
my $DOSURFACE = 0 ; 
my $DISTATOMS   = 1 ; 
my $BOUNDARYDELTA   = 1000 ; 
my $RECALIBRATECOUNT= 100 ; 
my $CNT1000 = 5 ;
my $addatomsclose2 = 0 ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($debugresidue,$ann,$config,$p1,$p2,$infile,$threshPD,$threshsign,$threshDist,$outfile,$readpotential,$which_tech,$listfile,$protein);
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;

my ($onlypolar,$radii,$before1,$before2);
$readpotential = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "debugresidue=s"=>\$debugresidue ,
            "onlypolar=i"=>\$onlypolar ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
            "distwater=f"=>\$distwater,
            "padding=i"=>\$padding ,
            "deltacube=i"=>\$deltacube ,
            "cutoffview=i"=>\$cutoffview ,
            "addatomsclose2=f"=>\$addatomsclose2 ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a p1 file name => option -p1 ") if(!defined $p1);


print "padding = $padding, distwater = $distwater, deltacube = $deltacube cutoffview = $cutoffview \n";

my $ofh = util_write($outfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

ConfigPDB_Init($config,$ofh);

my $pdb = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my @atoms1 = $pdb1->GetAtomsNoHet();
my ($minx,$miny,$minz);
my ($maxx,$maxy,$maxz);
$minx = $miny = $minz = 100000 ;
$maxx = $maxy = $maxz = -10000 ;
foreach my $a1 (@atoms1){
	my ($x,$y,$z) = $a1->Coords();
	next if(!defined $x);
	$maxx = $x if($x > $maxx); $maxy = $y if($y > $maxy); $maxz = $z if($z > $maxz);
	$minx = $x if($x < $minx); $miny = $y if($y < $miny); $minz = $z if($z < $minz);
}

### round off 
$maxx = int($maxx); $maxy = int($maxy); $maxz = int($maxz);
$minx = int($minx); $miny = int($miny); $minz = int($minz);

### pad 
$maxx = $maxx + $padding ; $maxy = $maxy + $padding ; $maxz = $maxz + $padding ;
$minx = $minx - $padding ; $miny = $miny - $padding ; $minz = $minz - $padding ;

print "maxx maxy maxz\n";
print "$maxx $maxy $maxz\n";
print "minx miny minz\n";
print "$minx $miny $minz\n";
my $diffx = $maxx - $minx ; my $diffy = $maxy - $miny  ; my $diffz = $maxz - $minz  ;
my $CNT = 0 ; 
my $pseudoatom = new Atom();
$pseudoatom->SetIdx(10000);


my $gridresilution = 1 ; 

my $allPointTable = {};
my @allPointList  = ();
my $boundary ={};
my $tablefilled = {};


##############################################################
### Accumulate all the grid points and mark the boundaries ###
##############################################################
my @Xs ;
my @Ys ;
my @Zs ;
my $CNTBoundary = 0 ;
foreach my $p (0..$diffx){
	my $X = $gridresilution * $p + $minx ; 
    foreach my $q (0..$diffy){
	    my $Y = $gridresilution * $q + $miny ; 
        foreach my $r (0..$diffz){
	        my $Z = $gridresilution * $r + $minz ; 

			push @Xs, $X ;
			push @Ys, $Y ;
			push @Zs, $Z ;
	        my $str = geom_MakeKeyFromCoord($X,$Y,$Z);
			$allPointTable->{$str} = {};
			push @allPointList,$str;
			$CNT++ ;

			## Mark the boundaries
			if ($p eq 0 || $q eq 0 || $r eq 0 || $p eq $diffx || $q eq $diffy || $r eq $diffz){
			    $boundary->{$str} = 1 ;
				$CNTBoundary++;
			}
        }
    }
}

print "There were $CNT points, and $CNTBoundary boundary points \n";


my $natoms = @atoms1 ;
my $ngridpoints = @allPointList ;
### choose the closest grid point as the atom - hence the distance
my $cntfilled = 0 ; 

print "Assigning $ngridpoints grid points to $natoms atoms - deltacube is $deltacube\n";
foreach my $a1 (@atoms1){
	my ($x,$y,$z) = $a1->Coords();
	my $NM = $a1->GetName();

	die if(!defined $x); ## why would this ever happen?

    my $intx = int($x); 
	my $inty = int($y); 
	my $intz = int($z);

	my $dx = abs($x - $intx) ; 
	my $dy = abs($y - $inty); 
	my $dz = abs($z - $intz) ; 

	$x = $dx > 0.5 ? $intx+1: $intx;
	$y = $dy > 0.5 ? $inty+1: $inty;
	$z = $dz > 0.5 ? $intz+1: $intz;

	my $Sx  = $x- $deltacube ; my $Ex  = $x + $deltacube ;
	my $Sy  = $y- $deltacube ; my $Ey  = $y + $deltacube ;
	my $Sz  = $z- $deltacube ; my $Ez  = $z + $deltacube ;
	my $found = 0 ;

	### there may be multiple atoms in each grid point
	my $str = geom_MakeKeyFromCoord($x,$y,$z);
	$tablefilled->{$str} = [] if(!defined $tablefilled->{$str} || !$tablefilled->{$str} );
	push @{$tablefilled->{$str}}, $a1 ; 

	my $distable = {};

	my @dl = ();
	my $dtable = {};
    foreach my $p ($Sx..$Ex){
        foreach my $q ($Sy..$Ey){
            foreach my $r ($Sz..$Ez){
                  $pseudoatom->SetCoords($p,$q,$r);
				   my $d = $pseudoatom->Distance($a1) ;
				   
	               my $str = geom_MakeKeyFromCoord($p,$q,$r);
				   #die if(!defined $allPointTable->{$str});
				   next if(!defined $allPointTable->{$str});
				   next if(defined $boundary->{$str});

				   if($d <= $distwater){
				      push @dl, $str ; 
				   }
            }
        }
    }



	foreach my $str (@dl){
		if(! exists $tablefilled->{$str}){
	       $tablefilled->{$str} = [] if(!defined $tablefilled->{$str} || !$tablefilled->{$str} );
	       push @{$tablefilled->{$str}}, $a1 ; 
		}
	}
}

my $tablesurface  = {};
my $tablesurfacevalues  = {};
my  $JJJ = (keys %{$tablefilled});

print "Scanning for each filled point to find surface, cutoffview = $cutoffview\n";
CreateTableSurfaceByScanningOnEachSide();

my $KKK = (keys %{$tablesurface});

my $cavity = {};
my $exposedsurface = {};

foreach my $k (keys %{$tablesurface}){
	my $v = $tablesurface->{$k};
	my $value = $tablesurfacevalues->{$k};
	next if(!$v);
	my @l = @{$v};
	if($value < $cutoffview){
		$cavity->{$k} = 1 ;
	}
	else{
		$exposedsurface->{$k} = 1 ;
	}
	foreach my $a (@l){
		my $nm = $a->GetNameSlashSep();
		#print "$nm VVV = $value\n";
	}
}

my $NCAV = (keys %{$cavity});
print " number of initial cavity points = $NCAV from $KKK surface points\n";

my $staticvalues = {};
my $recalibrated = {};

my $ncavity = (keys %{$cavity});
print "There were $ncavity cavity points \n";



### this looks close grid points to the cavity - and tries to increase it.
## dont remember why
if(1){
my @Addtocavity = Add2Cavity($cavity);
my $naddtocavity = @Addtocavity ; 
print "Adding $naddtocavity points to cavity\n";
foreach my $k (@Addtocavity){
	$cavity->{$k} = 1 ; 
}
}

if(defined $debugresidue){
    my @points1 ; 
	my $points = {};
    my $IFH = util_read($debugresidue);
    my $OFH = util_write("$debugresidue.out");
    while(<$IFH>){
         next if(/^\s*$/);
	     my ($nm,$key) = split ; 
	     push @points1, $key ;
		 $points->{$key} = 1 ;
    }
	foreach my $p (@points1){
	    my @l = geom_GetPointsAroundOnePoint($p,1);
		foreach my $k (@l){
		     if(exists $recalibrated->{$k}){
	             my @values = @{$recalibrated->{$k}};
	             my @origvalues = @{$staticvalues->{$k}};
			     my $iscavity = exists $cavity->{$k} ? 1 : 0 ; 
			     print $OFH "$k @values $iscavity\n";
			     print $OFH "$k @origvalues $iscavity\n";
		     }
		}
	}
	die "debugresidue $debugresidue done \n";
}


print "Assigning cavities\n";
my $numcavities = 0 ; 
my $CAVITIES = {};
my $CAVITIESATOMS = {};
my $donecavities = {};
foreach my $k (keys %{$cavity}){
	next if(exists $donecavities->{$k});
	$donecavities->{$k} = 1 ;

    $numcavities++;
	$CAVITIES->{$k} = $numcavities; 
	$CAVITIESATOMS->{$numcavities} = [] ;

	my @l ;
	push @l, $k ; 
	while(@l){
		my $newkey = shift @l ;
	    my @ll = geom_GetPointsTopDownNorthSouthEastWest($newkey);
		foreach my $l (@ll){
			next if(! exists $cavity->{$l});
	        next if(exists $donecavities->{$l});

	        $CAVITIES->{$l} = $numcavities; 
	        $donecavities->{$l} = 1 ;
	        push @l, $l ; 
			if(exists $tablesurface->{$l}){
		 	      push @{$CAVITIESATOMS->{$numcavities}} ,@{$tablesurface->{$l}}; 
			}
		}
	}
}

print "numcavities = $numcavities \n";

foreach my $k (keys %{$CAVITIESATOMS}){
	my @l = @{$CAVITIESATOMS->{$k}};
    my $tableres = {};
	foreach my $a (@l){
		my $resnum = $a->GetResNum();
		$tableres->{$resnum} = 1; 
	}
	my $N = (keys %{$tableres});
	next if($N eq 0);
	next if($N eq 1);
	print "$k has $N residues \n";
}
die ;



foreach my $k (keys %{$cavity}){
	if($cavity->{$k} eq 1 ){
		my @atoms ;
		my @ret ; 
  
        RecurseCreateCavity($k,\@ret,\@atoms);

		$numcavities++ ; 

		my $N  = @ret ; 
		$CAVITIES->{$numcavities} = {};
		$CAVITIES->{$numcavities}->{POINTS} = \@ret ; 
		$CAVITIES->{$numcavities}->{ATOMS} = \@atoms ; 
		$CAVITIES->{$numcavities}->{NUM} = $N ; 
	}
}



my @sorted = sort { $CAVITIES->{$a}->{NUM} <=> $CAVITIES->{$b}->{NUM} } (keys %{$CAVITIES});
foreach my $numcavities (@sorted){    
	my $N = $CAVITIES->{$numcavities}->{NUM} ;
	my $atoms = $CAVITIES->{$numcavities}->{ATOMS} ;
	my $points = $CAVITIES->{$numcavities}->{POINTS} ;
	my $done = {};
	my $doneatom = {};
	my $str = "";
	my $numres = 0 ; 
	my $natomsincavity = 0 ; 


	foreach my $a (@{$atoms}){
	    my $NM = $a->GetName();
		if(! exists $doneatom->{$NM}){
			$natomsincavity++ ; 
		}
	    my $nm = $a->GetResName();
	    my $number = $a->GetResNum();
		my $num = "$nm.$number";
		 if(! exists $done->{$num}){
	          $str = $str . " $num "  ;
			  $numres++;
		 }
	     $done->{$num} = 1 ;
        #$pseudoatom->SetCoords($x,$y,$z);
		#my $list = util_make_list($pseudoatom);
		#my ($junk,$neigh)  = $pdb1->GetNeighbourHoodAtom($list,2);
		#foreach my $a (@{$neigh}){
		    #my $num = $a->GetResNum();
		#}
	}
	print $ofh "$numcavities has size $N, and $numres residues, and $natomsincavity atoms - $str \n";
}

exit ; 
###############  END ################################


### this is quite a complicated function - so pay attention ### 
### As we go up, we dont want to go down - just east, west , north and south
sub RecalibrateValues{
	my ($k,$list,$idx,$dir,$otherdir) = @_ ; 
	my @l = @{$list} ; 
	my @origvalues = @{$staticvalues->{$k}};
	my $cnt = 0 ; 
	
	## this might as well have been while 1
	while($cnt < $RECALIBRATECOUNT){
		$cnt++ ; 
		if($dir eq 0){
	          $l[$idx] = $l[$idx] -1 ; 
		 }
		 else{
	          $l[$idx] = $l[$idx] + 1 ; 
		 }
	     my ($x,$y,$z) = @l ;
	     my $newkey = geom_MakeKeyFromCoord($x,$y,$z);
		 if(exists $boundary->{$newkey} || exists $tablefilled->{$newkey}){
		    $recalibrated->{$k} = \@origvalues;
		 	return ; 
		 }
		 my $values = $staticvalues->{$newkey};
		 my @values = @{$values};
		 foreach my $IDX (0..5){
		 	next if($IDX eq $idx || $IDX eq $otherdir);
			my $newval = $values[$IDX];
			$origvalues[$IDX] = $origvalues[$IDX] + int($newval/($cnt)) ; 
		 }
	}
}

sub MoveInOnedirection{
	my ($k,$list,$idx,$dir,$ret) = @_ ; 
	my @l = @{$list} ; 
	my $cnt = 0 ; 
	while(1){
		if($dir eq 0){
	          $l[$idx] = $l[$idx] -1 ; 
		 }
		 else{
	          $l[$idx] = $l[$idx] + 1 ; 
		 }

	     my ($x,$y,$z) = @l ;
		 #print "$x $y $z \n";
	     my $newkey = geom_MakeKeyFromCoord($x,$y,$z);
		 if(exists $boundary->{$newkey}){
		 	$cnt = $cnt+ $BOUNDARYDELTA ; 
			push @{$ret}, $cnt ; 
			return $cnt; 
		 }
		 elsif(exists $tablefilled->{$newkey}){
			push @{$ret}, $cnt ; 
		    return $cnt ; 
		 }
		 ## could make a move 
		 $cnt++ ; 
	}
}







sub usage{
my ($msg) = @_ ;
print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
die ;
}


sub GetSurface{
	my ($tf,$distsurface) = @_ ; 
	my $ts = {};
    foreach my $pointfilled (keys %{$tf}){
		my @l = geom_GetPointsAroundOnePoint($pointfilled,$distsurface);
		my $atom = $tf->{$pointfilled} ; 
		foreach my $k (@l){
	        $ts->{$k} = $atom if(! exists $tf->{$k});
		}
    }
    
    #foreach my $k (keys %{$ts}){
	    #print "surface $k \n";
    #}
	return $ts ; 
}



sub RecurseCreateCavity{
	my ($k,$ret,$atoms) = @_ ; 
    die "Did not expect to get to the boundary $k" if(exists $boundary->{$k});
	if(defined $cavity->{$k}){
	   if($cavity->{$k} eq 1 ){
	      $cavity->{$k} = 0 ; 
		  push @{$ret}, $k ;
	      my ($x,$y,$z) = geom_MakeCoordFromKey($k);
		  my $k1 = geom_MakeKeyFromCoord($x+1,$y,$z);

		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }
		  $k1 = geom_MakeKeyFromCoord($x-1,$y,$z);
		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }
		  $k1 = geom_MakeKeyFromCoord($x,$y+1,$z);
		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }
		  $k1 = geom_MakeKeyFromCoord($x,$y-1,$z);
		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }
		  $k1 = geom_MakeKeyFromCoord($x,$y,$z+1);
		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }
		  $k1 = geom_MakeKeyFromCoord($x,$y,$z-1);
		 if(!exists $boundary->{$k1}){
		  RecurseCreateCavity($k1,$ret,$atoms) if(!defined $cavity->{$k1} || $cavity->{$k1} eq 1);
		  }


	    }
    }		


	{
		  if($addatomsclose2){
	           my ($x,$y,$z) = geom_MakeCoordFromKey($k);
               $pseudoatom->SetCoords($x,$y,$z);
		       my $list = util_make_list($pseudoatom);
		       my ($junk,$neigh)  = $pdb1->GetNeighbourHoodAtom($list,$addatomsclose2);
		       foreach my $a (@{$neigh}){
		           push @{$atoms} ,$a ;
		       }
		  }

		 if( exists $tablesurface->{$k}){
		 	push @{$atoms} ,@{$tablesurface->{$k}}; 
		 }
		 elsif (exists $tablefilled->{$k}){
		 	push @{$atoms} ,@{$tablefilled->{$k}} if($tablefilled->{$k});
		 }
		return ;
	}
}

sub Add2Cavity{
	my ($CAV) = @_ ; 
    my @addtocavity ; 
    foreach my $k (keys %{$CAV}){
		# also includes the main point - which will be ignored in the next if loop
	    my @l = geom_GetPointsTopDownNorthSouthEastWest($k);
        foreach my $newkey (@l){
		       if(!exists $exposedsurface->{$newkey} && !exists $boundary->{$newkey} 
			             && ! exists $tablefilled->{$newkey} && !exists $CAV->{$newkey}){
			   		   my $ret = GetVisibilityValue($newkey);
					   if($ret < 1){
				            push @addtocavity, $newkey ;
					   }
					   else{
					   	$exposedsurface->{$newkey} = 1;
					   }
			   }
         }
    }
	return @addtocavity ; 
}

sub CreateTableSurfaceByScanningOnEachSide{

	my $CNTTTT = 0 ; 
	foreach my $k (keys %{$tablefilled}){
		
		my $ret = GetVisibilityValue($k);

	    if($ret){
		    $tablesurface->{$k} = $tablefilled->{$k} ;
		    $tablesurfacevalues->{$k} = $ret ;
	    }
	}

}


sub GetVisibilityValue{
	my ($k) = @_ ;
	    my ($x,$y,$z) = geom_MakeCoordFromKey($k);

	    my $ret = MoveOneSide($x,$y,$z,1,"x");
	    $ret = $ret + MoveOneSide($x,$y,$z,0,"x") ;
	    $ret = $ret + MoveOneSide($x,$y,$z,1,"y") ;
	    $ret = $ret + MoveOneSide($x,$y,$z,0,"y") ;
	    $ret = $ret + MoveOneSide($x,$y,$z,1,"z") ;
	    $ret = $ret + MoveOneSide($x,$y,$z,0,"z") ;

		return $ret ;

}

sub MoveOneSide{
    my ($x,$y,$z,$dir,$COORD) = @_ ;
	while(1){
		 if($COORD eq "x")  {
		    $x = $dir ? $x+1 : $x-1;
		 }
		 elsif ($COORD eq "y"){
		    $y = $dir ? $y+1 : $y-1;
		 }
		 elsif ($COORD eq "z"){
		    $z = $dir ? $z+1 : $z-1;
		 }
		my $newkey = geom_MakeKeyFromCoord($x,$y,$z);
		if(exists $boundary->{$newkey}){
			#print "found boundary $dir $COORD $newkey\n";
			return 1 ;
		}
		if(exists $tablefilled->{$newkey}){
			return 0 ;
		}
	}
}

if(0){

print "Computing static distances to boundary/atoms for each point\n";
my $IGNOREALREADY = 0 ;
foreach my $k (@allPointList){
	my $cnt1000 = 0 ;
	if(! exists $tablefilled->{$k} && ! exists $boundary->{$k}){
	    my ($x,$y,$z) = geom_MakeCoordFromKey($k);
	    my $list = util_make_list($x,$y,$z);
		my @ret ; 
	    #params for MoveInOnedirection  --> ($k,$list,$idx,$dir,$ret) = @_ ; 
		my $north = MoveInOnedirection($k,$list,0,0,\@ret);
	   	    $cnt1000++ if($north < 1000);
		my $south = MoveInOnedirection($k,$list,0,1,\@ret);
	   	    $cnt1000++ if($south < 1000);
		my $east = MoveInOnedirection($k,$list,1,1,\@ret);
	   	    $cnt1000++ if($east < 1000);
		my $west = MoveInOnedirection($k,$list,1,0,\@ret);
	   	    $cnt1000++ if($west < 1000);
		my $up = MoveInOnedirection($k,$list,2,1,\@ret);
	   	    $cnt1000++ if($up < 1000);
		my $down = MoveInOnedirection($k,$list,2,0,\@ret);
	   	    $cnt1000++ if($down < 1000);
		$staticvalues->{$k} = \@ret ; 
		#print "$x $y $z $north $south $east $west $up $down \n";

		if($cnt1000 > $CNT1000){
			$IGNOREALREADY++;
		}
	}
}
print "IGNOREALREADY = $IGNOREALREADY as they are > $CNT1000 \n";

if($RECALIBRATE){
print "recalibrating\n";
foreach my $k (@allPointList){
	if(! exists $tablefilled->{$k} && ! exists $boundary->{$k}){
	    my ($x,$y,$z) = geom_MakeCoordFromKey($k);
	    my $list = util_make_list($x,$y,$z);
		## todo remove values
		my $north = RecalibrateValues($k,$list,0,0,1);
		my $south = RecalibrateValues($k,$list,0,1,0);
		my $east = RecalibrateValues($k,$list,1,1,3);
		my $west = RecalibrateValues($k,$list,1,0,2);
		my $up =   RecalibrateValues($k,$list,2,1,5);
		my $down = RecalibrateValues($k,$list,2,0,4);
	}
}
}
else{
$recalibrated = $staticvalues;
}

print "Figuring out which are cavities - cnt1000 should be greater than $CNT1000\n";
foreach my $k (@allPointList){
	if(! exists $tablefilled->{$k} && ! exists $boundary->{$k}){
	   my @values = @{$recalibrated->{$k}};
	   my $cnt1000 = 0 ; 
	   foreach my $i (@values){
	   	    $cnt1000++ if($i < 1000);
	   }

	   if($cnt1000 > $CNT1000){
		   $cavity->{$k} = 1  ;
	   }
	   else{
		   	$, = " " ;
		   	print $ofh " not cavity @values \n";
		}
	}
}
}
