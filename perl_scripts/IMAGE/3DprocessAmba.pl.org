#!/usr/bin/perl -w 
use lib '/home/b/Bio/Code/perl_scripts/perl_scripts/';
use MyGeom;
use MyMagick;
use MyUtils;
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use PDB;
use POSIX ;
use Algorithm::Combinatorics qw(combinations) ;
use Math::Geometry ;
use Math::Geometry::Planar;
my $polygon = Math::Geometry::Planar->new; 
   #$contour = Math::Geometry::Planar->new; creates a new contour object;


#use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
$, = " ";
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$display,$listfile,$protein);
my ($figurerealstack,$resolution,$ChangeColorOfRest,$justborder,$justidentify,$csv,$contourcolor,$color,$distfromperiphery,@expressions,$from,$to,@colornms,$specifiedcolor);
my $verbose = 1 ;
my $NUMITERS = 100 ; 
my $DISTHRESH = 1 ; 
GetOptions(
            "from=s"=>\$from ,
            "to=s"=>\$to ,
            "ChangeColorOfRest"=>\$ChangeColorOfRest ,
            "csv=s"=>\$csv ,
            "protein=s"=>\$protein ,
            "contourcolor"=>\$contourcolor ,
            "display"=>\$display ,
            "infile=s"=>\$infile ,
            "specifiedcolor=s"=>\$specifiedcolor ,
            "listfile=s"=>\$listfile ,
            "color=s"=>\@colornms ,
            "outfile=s"=>\$outfile ,
            "justidentify=s"=>\$justidentify ,
            "expr=s"=>\@expressions,
            "dist=f"=>\$DISTHRESH ,
            "resolution=f"=>\$resolution ,
            "justborder"=>\$justborder ,
            "figurerealstack=s"=>\$figurerealstack ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
my $RATIO = 0.3 / $resolution ; 
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a specifiedcolor ile name => option -specifiedcolor ") if(!defined $specifiedcolor );
usage( "Need to give a listfile ile name => option -listfile ") if(!defined $listfile );
usage( "Need to give a color => option -color ") if(!@colornms);
my $ifh = util_read($infile);
#usage( "Need to give a from -option -from  ") if(!defined $from);
#usage( "Need to give a to -option -to  ") if(!defined $to);
$csv = "$specifiedcolor.rawdata";
my $ofhcsv = util_write($csv);
my $ofhperi = util_write("out.peri");
my $ofhcentre = util_write("out.centre");
my $ofhadp = util_write("out.addedpoints");
my $ofhvol = util_write("volume");

my @list= util_read_list_sentences($listfile);

my $colornms = util_make_table(\@colornms);

my $image = new MyMagick($infile);

my $w = $image->GetWidth();
my $h = $image->GetHeight();


my $obj = $image->GetObj();

my $RlistF = $image->GetRowInfoContour(0,$colornms);
my $RlistR = $image->GetRowInfoContour(1,$colornms);
#### ALL points are actually all points on the contour
my ($minX,$minY,$maxY,$maxX,$CminX,$CminY,$CmaxY,$CmaxX,$midX,$midY,$AREACALCBYPIXEL,$insidetable,@allpoints) = GetAllPoints($image,$RlistF,$RlistR);
my @allpointstemp = @allpoints  ;
#if(defined $justborder){
if(1){
    my @l = @allpoints  ;
    while(@l){
		 my $x = shift @l;
		 my $y = shift @l;
         $image->SetColorToPixelByColorName($x,$y,"red");
	}
	$image->Write("border.png",$display);
	#die ;
}



my ($maxDist,$minDist,$pointAX,$pointAY,$pointBX,$pointBY,$newmidX,$newmidY) = geom_GetABForEllipse($midX,$midY,\@allpoints);

print "W = $w H = $h maxdist = $maxDist minDist = $minDist \n";
my $ofhmax = util_write("max");
my $MM = int($maxDist);
print $ofhmax "$MM \n";



my $cnt = 0; 
foreach my $i (@list){
	$cnt++ ; 
	last if($i eq $infile);
}
my $MIDCNT = $cnt ;
my $MID = $list[$cnt];

my $junkcnt = 0 ;

my @threeDpoints ; 
my @newcolors ; 
#push @newcolors , "red";
push @newcolors , "blue";
my $newcolornms = util_make_table(\@newcolors);

$cnt = 0 ; 
my $foundblue = 0 ; 
my $startblue ;
my $endblue ;
my $volume = 0 ;
my $tableMinDist = {};

if(defined $ChangeColorOfRest){
    foreach my $i (@list){
	   ChangeColorOfRest($image,$h,$insidetable,"black");
	   my $nmchanges = "$i.colorchanged.png";
	   $image->Write($nmchanges,$display);
	}
	die ;
}

if(defined $figurerealstack){
	my $ofhnewlist = util_write($figurerealstack);
	my @fl ; 
	my @sl ; 
	my $N = @list -1  ;
	my $HALF = $N/2 ; 
	my $HH = $HALF - 1 ; 
    while($HH ge 0){
	   my $idx = $HH-- ;
		
		my $i = $list[$idx];
        my $newimage = new MyMagick($i);
	    ChangeColorOfRest($newimage,$h,$insidetable,"black");
        my (@blue) = $newimage->FindCoordsForColor("blue");
		my $NBLUE = @blue ;
		print "found blue $NBLUE \n";
	    if(!@blue){
			last ;
	    }
		else{
			push @fl, $i ;
		}
	}


    foreach my $idx ($HALF..$N){ 
		
		my $i = $list[$idx];
        my $newimage = new MyMagick($i);
	    ChangeColorOfRest($newimage,$h,$insidetable,"black");
        my (@blue) = $newimage->FindCoordsForColor("blue");
		my $NBLUE = @blue ;
		print "found blue $NBLUE \n";
	    if(!@blue){
			last ;
	    }
		else{
			push @sl, $i ;
		}
	}
	my @revfl = reverse @fl ;
	push @revfl, @sl ;
	foreach my $i (@revfl){
		print $ofhnewlist "$i\n";
	}
	my $NN = @revfl ;
	die " N = $N , NN = $NN \n";

}


foreach my $i (@list){

	my $dist = abs($MIDCNT - $cnt);
    my $newimage = new MyMagick($i);


    my $w = $newimage->GetWidth();
    my $h = $newimage->GetHeight();
	ChangeColorOfRest($newimage,$h,$insidetable,"black");
	#my $nmchanges = "$i.colorchanged.png";
	#$newimage->Write($nmchanges,$display);

	if(defined $justidentify){
           my (@l) = $newimage->FindCoordsForColorAmbaRed($justidentify);
		   if(@l){
		   	 print "Identified color $justidentify\n";
		   }
		   while(@l){
		       my $x = shift @l;
		       my $y = shift @l;
               $newimage->SetColorToPixelByColorName($x,$y,"white");
		   }
	       $newimage->Write("$i.identify$justidentify.png",$display);
		   next ;
	}

    my (@speccolor) = $newimage->FindCoordsForColorAmbaRed($specifiedcolor);
    my (@blue) = $newimage->FindCoordsForColor("blue");
	$volume = $volume + @blue + @speccolor ; 
	if(!$foundblue && @blue){
		$foundblue = 1 ; 
		$startblue = $cnt ; 
	}
	if($foundblue && !@blue){
		$endblue = $cnt ; 
		print "Ending since we dont see blue\n";
		last ; 
	}


	### Just see if red is gettting identified

	my $nspeccolor = @speccolor ;
    print "W = $w H = $h dist=$dist i=$i midcnt=$MIDCNT speccolor=$nspeccolor \n";
	my $found = 0 ; 
	my @LLL ;  ## just fill once
	while(@speccolor){
		my $x = shift @speccolor ;
		my $y = shift @speccolor ;
		my $str = geom_MakeKeyFromCoord($x,$y);
		my $threeDstr = geom_MakeKeyFromCoord($x,$y,$cnt);
		next if(! exists $insidetable->{$str});

		print "\t found something in this slice - and the marked area\n";

		if(!$found){
            my $F = $newimage->GetRowInfoContour(0,$newcolornms);
            my $R = $newimage->GetRowInfoContour(1,$newcolornms);
            my ($minX,$minY,$maxY,$maxX,$CminX,$CminY,$CmaxY,$CmaxX,$midX,$midY,$AREACALCBYPIXEL,$III,@XXX) = GetAllPoints($newimage,$F,$R);
			@LLL = @XXX ; ## just fill once
		    my @l = @LLL ;
	        while(@l){
		        my $X = shift @l;
		        my $Y = shift @l;
                $newimage->SetColorToPixelByColorName($X,$Y,"white");
		    }

		}

		my $MINDIST = 1000 ; 
		my @l = @LLL ;
	    while(@l){
		    my $X = shift @l;
		    my $Y = shift @l;
		    my $dist = geom_Distance_2D($X,$Y,$x,$y);
			$MINDIST = $dist if($dist < $MINDIST);
		}

		#if($MINDIST eq 0){
			#print "Ignoring point as it is on the peripher\n";
			#next;
		#}

		print $ofhperi "$threeDstr $MINDIST \n";
		$tableMinDist->{$threeDstr} = $MINDIST;

        $newimage->SetColorToPixelByColorName($x,$y,"black");
		$found = 1 ; 

		push @threeDpoints, $threeDstr ; 


	    #my $D = geom_Distance_2D($midX,$midY,$x,$y);
		#$D = $D * $SCALE ;
		#my $ZZ = sqrt($XX*$XX + $D*$D);
	}
	if($found){
	    $newimage->Write("found/$i.found$cnt.png",$display);
	}

	#$slice2points->{$cnt} = \@threeDpoints ;
	$cnt++ ;
}
$endblue = $cnt if(!defined $endblue);
print $ofhvol "volume = $volume\n";


my $threeDCount = @threeDpoints ;
print "threeDpoints = $threeDCount\n";
my @addedpoints = PrepocessPoints(\@threeDpoints);
foreach my $adp (@addedpoints){
	print $ofhadp "$adp\n";
}
push @threeDpoints, @addedpoints ;

#################################################
#### Calculate mid slice - and thus centre
#################################################
my $midblue = ($startblue + $endblue)/2 ;
print "SSSSSSSS midblue = $midblue \n";
my $midblueimage = $list[$midblue];
my $threeDstrCentre ;

### just to add a scope
if(1){
        my $IMAGE = new MyMagick($midblueimage);
	    ChangeColorOfRest($IMAGE,$h,$insidetable,"black");
        my $RlistF = $IMAGE->GetRowInfoContour(0,$newcolornms);
        my $RlistR = $IMAGE->GetRowInfoContour(1,$newcolornms);
        #### ALL points are actually all points on the contour
        my ($minX,$minY,$maxY,$maxX,$CminX,$CminY,$CmaxY,$CmaxX,$midX,$midY,$AREACALCBYPIXEL,$III,@l) = GetAllPoints($IMAGE,$RlistF,$RlistR);
		$threeDstrCentre = geom_MakeKeyFromCoord($midX,$midY,$midblue);
		print $ofhcentre "#centre = $threeDstrCentre, name of file = $midblueimage \n";
		print  "#centre = $threeDstrCentre, name of file = $midblueimage \n";
        while(@l){
		     my $x = shift @l;
		     my $y = shift @l;
             $IMAGE->SetColorToPixelByColorName($x,$y,"white");
	    }
        $IMAGE->SetColorToPixelByColorName($midX,$midY,"green");
	    $IMAGE->Write("centre.png",$display);

}

foreach my $p (@threeDpoints){
	     my ($X,$Y,$Z) = geom_MakeCoordFromKey($threeDstrCentre);
	     my ($x,$y,$z) = geom_MakeCoordFromKey($p);
            
		 $z = $RATIO * $z ; 
		 $Z = $RATIO * $Z ; 
		 my $d = geom_Distance_3D($x,$y,$z,$X,$Y,$Z);
		 print $ofhcentre "$p $d\n";
}

#################################################
#### Group points
#################################################
my $grpcnt = 0 ; 
my $grps = {} ; 
my $POINTSINGROUP = {};
my $doit = 1 ;
while($doit){
	my ($grp) = CreateGroup(\@threeDpoints,$POINTSINGROUP);
	$grps->{$grpcnt} = $grp ; 
	$doit = @{$grp};
	if($doit){
	    print "Created Grp $grpcnt with $doit items\n";
		my $grpname = "GROUP$grpcnt";
		print $ofhperi "$grpname ";
	    foreach my $x (@{$grp}){
		    print "\t $x \n";
			print $ofhperi " $x ";
	    }
		print $ofhperi "\n";
	}
	$grpcnt++ ;
}

print "Done with processing\n";
exit ; 

sub PrepocessPoints{
	my ($listP) =@_ ; 
	my $tab = {};
    foreach my $P (@{$listP}){
	         my ($X,$Y,$Z) = geom_MakeCoordFromKey($P);
		     my $str = geom_MakeKeyFromCoord($X,$Y);
			 $tab->{$str} = [] if(!defined $tab->{$str});
			 push @{$tab->{$str}} , $Z ; 
	}

	my @ADDEDPOINTS ;
	foreach my $k (keys %{$tab}){
	      my @l = @{$tab->{$k}} ; 
		  my @sl = sort { $a <=> $b } @l ; 

		  while(@sl>1){
		  	  my $i = shift @sl ;
			  my $j = $sl[0];
			  if(abs($i - $j) > 1 ){
			  	      print "Add $i + 1 to $k  \n";
	                  my ($X,$Y) = geom_MakeCoordFromKey($k);
		              my $origstr = geom_MakeKeyFromCoord($X,$Y,$i);
					  my $MINDIST = $tableMinDist->{$origstr};
		              my $str = geom_MakeKeyFromCoord($X,$Y,$i+1);
		              print $ofhperi "$str $MINDIST \n";
					  push @ADDEDPOINTS,$str ;

			  }
		  }

	}
	return @ADDEDPOINTS ;
}

sub CreateGroup{
	my ($listP,$pointsingroup) =@_ ; 
	my @pointsingroup ;
	my $added = 1 ;
	while($added){
		$added = 0 ; 
        foreach my $P (@{$listP}){
			next if(exists $pointsingroup->{$P});

			if(!@pointsingroup){
				push @pointsingroup, $P; 
				$pointsingroup->{$P} = 1 ;
			}
			else{
                foreach my $p (@pointsingroup){
				   next if($p eq $P);
	               my ($X,$Y,$Z) = geom_MakeCoordFromKey($P);
	               my ($x,$y,$z) = geom_MakeCoordFromKey($p);
            
		           $z = $RATIO * $z ; 
		           $Z = $RATIO * $Z ; 
		           my $d = geom_Distance_3D($x,$y,$z,$X,$Y,$Z);
        
	               #print "HH = $d\n";
		           if($d < $DISTHRESH){
	                   push @pointsingroup, $P ; 
					   $added = 1 ; 
					   $pointsingroup->{$P} = 1 ;
		           }
	            }
			 }
	    }
	}

    my $done = {};
	my @l ;
	foreach my $x (@pointsingroup){
		if(!exists $done->{$x}){
			push @l, $x ; 
			$done->{$x} = 1 ;
		}
	}

	return  (\@l);

}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

sub GetAllPoints{
	my ($IMAGE,$listF,$listR) = @_ ; 
	my $AREACALCBYPIXEL = 0 ;
    my $n = @{$listF} -1 ; 
	################################################################
	########## This gets the bounding coordinates ###################
	################################################################
    my ($minX,$minY,$maxY,$maxX) ;
    my ($CminX,$CminY,$CmaxY,$CmaxX) ;
	$CmaxX = 0 ;
	$CminX = 1000000 ;
	my ($prevStartX,$prevEndX) ; 
	my (@allpoints);
	
	my $Insidetable = {};
    foreach my $i (0..$n){
	    next if($listR->[$i]  == -1 || $listF->[$i]  == -1);
    
        my $start = $listF->[$i] ;
        my $end = $listR->[$i] ;
        my $y = $i+1 ;
		if(!defined $prevStartX){
                foreach my $XXX ($start..$end){
					push @allpoints, $XXX ;
					push @allpoints, $y ;
				}
		}
		if($start < $CminX){
			$CminX = $start ; 
			$CminY = $y;
		}
		if($end > $CmaxX){
			$CmaxX = $end ; 
			$CmaxY = $y;
		}

	    $minY = $y if(!defined $minY);
	    $maxY = $y ;
    
	    $minX = $start if(!defined $minX);
	    $maxX = $start ;
    
		my $x = $start ;
		push @allpoints, $x ;
		push @allpoints, $y ;

		$x = $end ;
		push @allpoints, $x ;
		push @allpoints, $y ;


		############################################
		#### this makes the line continuous #########
		############################################
		if(defined $prevStartX){
				my $A =$prevStartX > $start ?  $start : $prevStartX ;
				my $B =$prevStartX > $start ? $prevStartX :$start  ;
                foreach my $XXX ($A..$B){
					push @allpoints, $XXX ;
					push @allpoints, $y ;
				}
				$A =$prevEndX > $end ?  $end : $prevEndX ;
				$B =$prevEndX > $end ? $prevEndX :$end  ;
                foreach my $XXX ($A..$B){
					push @allpoints, $XXX ;
					push @allpoints, $y ;
				}
		}
		$prevStartX = $start;
		$prevEndX = $end;
		

        foreach my $x ($start..$end){
           $AREACALCBYPIXEL++;
	       my $str = geom_MakeKeyFromCoord($x,$y);
	       $Insidetable->{$str} = 1 ;
           #$image->SetColorToPixelByColorName($x,$y,"green");
	    }
    }
	
    
   ## there might be a gap in the last line - fill it up
	my @allpointstemp = @allpoints  ;
	my $MMMinx = 100000;
	my $MMMaxx = 0  ;
	while(@allpointstemp){
			my $x = shift @allpointstemp ;
			my $y = shift @allpointstemp ;
			if($y eq $maxY){
				if($x < $MMMinx){
					$MMMinx = $x ; 
				}
				if($x > $MMMaxx){
					$MMMaxx = $x ; 
				}
			}
	}
	$MMMaxx = $MMMaxx - 1;
	$MMMinx = $MMMinx + 1;
    foreach my $x ($MMMinx..$MMMaxx){
		push @allpoints, $x ;
		push @allpoints, $maxY ;
	}

	my @allpointstemp = @allpoints  ;
	while(@allpointstemp){
			my $x = shift @allpointstemp ;
			my $y = shift @allpointstemp ;
            #$image->SetColorToPixelByColorName($x,$y,"green");
	}
	#$image->Write("tttt.png",$display);
	#die ;

	#####################################################
	## There are two midpoints - which one to choose? 
	##################################################
    my $midX = ($minX + $maxX)/2 ;
    my $midY = ($minY + $maxY)/2 ;
    #my $midX = (($minX + $maxX)/2 + ($CminX + $CmaxX)/2)/2 ;
    #my $midY = (($minY + $maxY)/2 + ($CminY + $CmaxY)/2)/2 ;
    print "midx =  $midX = ($minX + $maxX)/2 ; \n";
    print "midy  $midY = ($minY + $maxY)/2 ; \n";


    #$image->SetColorToPixelByColorName($minX,$minY,"red");
    #$image->SetColorToPixelByColorName($maxX,$maxY,"red");


    #$image->SetColorToPixelByColorName($CminX,$CminY,"red");
    #$image->SetColorToPixelByColorName($CmaxX,$CmaxY,"red");


	my $N = @allpoints; 
	print "Countor has $N points for GetAllPoints\n";
	return ($minX,$minY,$maxY,$maxX,$CminX,$CminY,$CmaxY,$CmaxX,$midX,$midY,$AREACALCBYPIXEL,$Insidetable,@allpoints);
}

sub ChangeColorOfRest{
	my ($IMAGE,$height,$INSIDETABLE,$color) = @_ ; 
    my $info ;
    foreach my $i (1..$height){
        my $p = $IMAGE->GetRowStraightorReversed(0,$i);
        while(@{$p}){
	        my $l = shift @{$p};
		    my ($r,$g,$b,$x,$y) = @{$l} ;
		    my $str = geom_MakeKeyFromCoord($x,$y);
		    next if( exists $INSIDETABLE->{$str});
            $IMAGE->SetColorToPixelByColorName($x,$y,$color);
        }
	}
}
