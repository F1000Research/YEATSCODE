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
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
use Scalar::Util qw(looks_like_number);


use Math::Trig;




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($aalist,$cutoff,$exception,$ann,$config,$p1,$p2,$infile,$score,$ignorepro,$outfile,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;
my $DBCA ;
my $DBCB ;
my $DBCN ;

$cutoff = 0.012 ;
my $MINSAMPLE = 1 ;
$ignorepro = 0 ;

my ($verify,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "aalist=s"=>\$aalist ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "cutoff=f"=>\$cutoff ,
            "config=s"=>\$config,
            "score=s"=>\$score,
            "ignorepro=i"=>\$ignorepro,
            "radii=i"=>\$radii ,
            "exception=s"=>\$exception ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a listfile -option -list  ") if(!defined $listfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $ofh= util_append("$outfile");

usage( "Need to give a aalist => option -aalist ") if(!defined $aalist);
my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my ($aainfo,$grplist) = util_ParseAAGroups($aalist);

ConfigPDB_Init($config);

my $prevnorm ; 
my $proteincnt = 0 ; 
my @alldistances ; 
my @nm2print ; 
foreach my $protein (@list){
    my $finalscoreCA = 0 ; 
    my $finalscoreCB = 0 ; 
    my $finalscoreCN = 0 ; 
    my $pdCAscore = 0 ; 
    my $pdCBscore = 0 ; 
    my $finalscorestr = "";
    #ConfigPDB_Init($config);
    my $i = $protein ;
    my @proteins ; 
    push @proteins, $i ; 

	print "Checking $APBSDIR/$i/$i.pqr\n";
    next if( ! -e "$APBSDIR/$i/$i.pqr");
    next  if( ! -e "$APBSDIR/$i/pot1.dx.atompot" && ! -e "$APBSDIR/$i/pot1.dx.atompot");
    my @info = util_ReadPdbs($PDBDIR,$APBSDIR,1,@proteins) ; 
    my $info = shift @info ;
    my $pdb1 = $info->{PDBOBJ};
    my $pqr1 = $info->{PQR};
    my $pots1 = $info->{POTS};

    my $PWD = cwd;

    my $pdb = "$PDBDIR/$protein.pdb";
    my $pdb1 = new PDB();
    $pdb1->ReadPDB($pdb);
    
    
    my @res = $pdb1->GetResidues();
    my $N = @res;
    my $prevres ; 
    my $prevpd = 0  ; 
    my @potCA ; 
	my $bad = 0 ; 
	my $N = @res ; 

	my $cnt = 0 ;
	my $numbercomparedCA = 0 ;
	my $numbercomparedCB = 0 ;
	my $numbercomparedCN = 1 ;
	my @pdCA ; 
	my @pdCB ; 
	my @pdCN ; 
	my @NM ; 
	my $good = 1 ; 
	$proteincnt++;


my $colortable = {};
my @NP = qw(G P A V L I M F );
foreach my $k (@NP){
	$colortable->{$k} = "yellow";
}

my @pos = qw (H K R);
foreach my $k (@pos){
	$colortable->{$k} = "red";
}

my @neg = qw (E D); 
foreach my $k (@neg){
	$colortable->{$k} = "blue";
}

my @misc = qw(Q N Y W C S T);
foreach my $k (@misc){
	$colortable->{$k} = "green";
}




	my $scoreDistCA = 0 ;
	my $maxscoreDistCA = 0 ;
	my $numberScoreDistCA = 0 ;
	
	my $HELIXDIFF = 4 ;
	my @COORDS ;
	my @REACTIVEATOMSLIST ;
	my @lol ;
	$lol[0] = [];
	$lol[1] = [];
	$lol[2] = [];
	$lol[3] = [];
    my $done = {};
    my $initval = 90 ; 
    my $rad = 5 ; 
	my $loopcnt = 0 ;
	foreach my $r (@res){
	    die if($r->GetAtomStr() eq "HETATM");
	    die if($r->GetName() eq "HOH");
		my $single =  $r->PrintSingleLetter($pdb1);
		my $num = $r->GetResNum();
		my $idxcnt = $loopcnt % 4 ;

        my $val = $initval - $loopcnt * 100 ; 
        my $x = util_format_float($rad * cos(deg2rad($val)),1) + 0 ;
        my $y = util_format_float($rad * sin(deg2rad($val)),1) + 0 ;
        my $str = "$x.$y";
        while(exists $done->{$str}){
            #print "$str existed\n";
            $y = $y + 1 ;	
            $str = "$x.$y";
        }
        $done->{$str} = 1 ;
		my $color = $colortable->{$single};
        print $ofh "\\node[fill=$color!100,draw=blue,very thick] (FinalNode) at ($x, $y) {$single$num} ; \n" ;

		$loopcnt++;
		push @{$lol[$idxcnt]}, $r ;
	}


	my @seq ; 
	my $finalseq = "";
	foreach my $idx (0..3){
		my @l = @{$lol[$idx]};
		print "$idx --> ";
		my @STACK ;

		my $seq = "";
		foreach my $r (@l){
		   my $single =  $r->PrintSingleLetter($pdb1);
		   my $num = $r->GetResNum();
		   my $str = $single . $num . " ";
		   print $str ;
	       die "Expected aa $single" if(!defined $aainfo->{$single});
	       my $xt = $aainfo->{$single}; 
		   $seq = $seq . $xt ;


		     my $CA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"CA");
			 my $RAnm = $pdb1->GetReactiveAtomFromAtom($CA);
		     my $RA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),$RAnm);
			 my $X = MakeVectorFrom2Points_Atoms($CA,$RA);
			 push @STACK, $X ;
			 if(@STACK > 1){
			 	my $N = @STACK ;
				my $v1 = $STACK[$N -1 ];
				my $v2 = $STACK[0];
			    my $angle = util_format_float(geom_AngleBetweenTwoVectors($v1,$v2),1);
				#print " (a=$angle ) ";
			 }


		}
		push @seq, $seq ;
		$finalseq = $finalseq . $seq . "M";

		print "  \n";
	}

	foreach my $seq (@seq){
		print "S = $seq \n";
	}
	print $ofh ">$i\n";
	print $ofh "$finalseq \n";


	my $failedcnt = 0 ;
	my $DIFFDIST =0;
    while(@res){
		
		my @RES ;
		my $res = shift @res ;
		push @RES, $res ;
		last if(@res < $HELIXDIFF);
		foreach my $i (0..3){
			die if(!defined $res[$i]);
			push @RES, $res[$i];
		}
		my @SINGLELETTER ;
		foreach my $r (@RES){
		    push @SINGLELETTER, $r->PrintSingleLetter($pdb1);
		}
		my @NUM ;
		foreach my $r (@RES){
		    my $num = $r->GetResNum();
			push @NUM,$num ;
		}




		if($DIFFDIST > 0.5){
			warn "Sorry, does not seem to be a helix - large distance beween CA and O";
		}


		my @CA ;
		my @C ;
		my @N ;
		my @O ;
		my @CB ;
		my @RA ;
		foreach my $r (@RES){
		     my $CA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"CA");
		     my $C = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"C");
		     my $N = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"N");
		     my $O = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"O");
		     my $CB = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"CB");
			 my $RAnm = $pdb1->GetReactiveAtomFromAtom($CA);
		     my $RA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),$RAnm);


		     push @CA, $CA ;
		     push @C, $C ;
		     push @N, $N ;
		     push @O, $O ;
		     push @CB, $CB ;
		     push @RA, $RA ;
		}
	    push @REACTIVEATOMSLIST, $RA[4];

	    #my @coord1 = $O[0]->Coords();
		#my @coord2 = $N[4]->Coords();
		my $a1 = $O[0];
		my $a2 = $N[4];
		my $p1 = vector( $a1->x(), $a1->y(), $a1->z() );
		my $p2 = vector( $a2->x(), $a2->y(), $a2->z() );
		my $P = $p1 - $p2 ;
		push @COORDS, $P;


       #my $potO1 = util_GetPotForAtom($O1,$pqr1,$pots1) *1  ;
       #my $potN2 = util_GetPotForAtom($N2,$pqr1,$pots1)  *1 ;
	   #my $potdiff = $potO1 - $potN2 ;


        my $distO1N2 = util_format_float($pdb1->DistanceAtoms($O[0],$N[4]),1);
        my $distCA = util_format_float($pdb1->DistanceAtoms($CA[0],$CA[4]),1);
        my $distOO = util_format_float($pdb1->DistanceAtoms($O[0],$O[4]),1);
		my $NM = $SINGLELETTER[0] . $NUM[0] . "-". $SINGLELETTER[4] . $NUM[4];

		my $angleBetweenO0N4 = util_format_float(-1,1); ;
		if(@COORDS>1){
			my $NN = @COORDS;
			my $a = $COORDS[$NN-1]->norm;
			my $b = $COORDS[$NN-2]->norm;
			$angleBetweenO0N4 = geom_AngleBetweenTwoVectors($a,$b);
		}

		if(@REACTIVEATOMSLIST>2){
			my $NN = @REACTIVEATOMSLIST;
			my $a = $REACTIVEATOMSLIST[$NN-1];
			my $b = $REACTIVEATOMSLIST[$NN-2];
			my $c = $REACTIVEATOMSLIST[$NN-3];
			my $X = MakeVectorFrom2Points_Atoms($a,$b);
			my $Y = MakeVectorFrom2Points_Atoms($b,$c);
			#$angle = geom_AngleBetweenTwoVectors($X,$Y);
		}
	    my $anglebetweenRandCA = geom_AngleBetweenTwoVectors_4Atoms($RA[0],$CA[0],$RA[4],$CA[4]);
	    my $anglebetweenOandCA = geom_AngleBetweenTwoVectors_4Atoms($O[0],$CA[0],$O[4],$CA[4]);

        #my $sortNM = util_sortsingleString($NM);
        my $distTmp = util_format_float($pdb1->DistanceAtoms($RA[4],$RA[0]),1);
        my $distTmp1 = util_format_float($pdb1->DistanceAtoms($RA[1],$RA[0]),1);
		if($distO1N2 > 4){
			$failedcnt++;
		}
		print "nextinhelix=$distTmp nextinnuber=$distTmp1 nm = $NM: $distO1N2 $distCA $distOO $angleBetweenO0N4 $anglebetweenRandCA $anglebetweenOandCA \n";
		$DIFFDIST = abs($distCA - $distOO);
		if($failedcnt > 2){
			die "Sorry, does not seem to be a helix";
		}
    }



}
print "TIKZ o/p written to $outfile \n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}


