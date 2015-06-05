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
my $verbose = 0 ;
my $DBCA ;
my $DBCB ;
my $DBCN ;

$ignorepro = 0 ;
my $HELIXDIFF = 4 ;
my $HBONDDIST = 3.7 ;
my $DISTO1N2 = 4 ;

my ($verify,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "p1=s"=>\$p1 ,
            "aalist=s"=>\$aalist ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "protein=s"=>\$protein ,
            "maxdist=f"=>\$maxdist ,
            "config=s"=>\$config,
            "score=s"=>\$score,
            "ignorepro=i"=>\$ignorepro,
            "radii=i"=>\$radii ,
            "DISTO1N2=f"=>\$DISTO1N2 ,
            "verbose=i"=>\$verbose ,
            "exception=s"=>\$exception ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $ofh= util_write("$outfile");

usage( "Need to give a aalist => option -aalist ") if(!defined $aalist);

my ($aainfo,$grplist) = util_ParseAAGroups($aalist);

ConfigPDB_Init($config);
my $ofhAllNames= util_append("allnames");

my $prevnorm ; 
my $proteincnt = 0 ; 
my @alldistances ; 
my @nm2print ; 

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my $ofhstats = util_append("stats");

my $helixnumber = "0";
my @frags = $pdb1->BreakIntoContinuosFragments();

foreach my $frag (@frags){
   my @indices = IdentifyHelix(@{$frag});
   my @merge ; 
   my $prevend ; 
   while(@indices){
   		my $start = shift @indices;
   		my $end = shift @indices;
		my $s = $start->GetResNum();
		my $e = $end->GetResNum();

        if(defined $prevend && $s < $prevend){
			print "Merging helix $s to $e\n";
			my $N = @merge ;
			$merge[$N-1]= $end ;
			next ;
		}
		push @merge, $start ;
		push @merge, $end ;
		$prevend = $e ; 


   }
   while(@merge){
   		my $start = shift @merge;
   		my $end = shift @merge;
		my $s = $start->GetResNum();
		my $e = $end->GetResNum();
		if(abs($s - $e) <= 4){
			print "Helix is less than 4, so ignoring $s $e\n";
			next ;
		}
		my $tableofres = {};
		foreach my $i ($s..$e){
			$tableofres->{$i} = 1 ;
		}
	    my $PPP = "$protein.helix$helixnumber";
	    my $OOO = "$PPP.pdb";
		print $ofhAllNames "$PPP\n";
        my $ofhhelix= util_write($OOO);
		print "Writing helixnumber $helixnumber to file $OOO : residues $s to $e \n";
		$pdb1->ReadPDBAndWriteSubset($pdb,$tableofres,$ofhhelix);
		close($ofhhelix);
		$helixnumber++;
	}
}


sub IdentifyHelix{
    my (@res)=@_;
    my $N = @res;
	if($N < 4){
		return ; 
	}
	
	my @HELIXSTART  ; 
	my @HELIXEND  ; 
	my @indices ; 
    while(@res){
		
		my @RES ;
		my $res = shift @res ;
		my $num = $res->GetResNum();
		push @RES, $res ;
		if(@res < $HELIXDIFF){
			if(@HELIXSTART){
				my $start = $HELIXSTART[0];
				my $end = $HELIXEND[0];
		        my $s = $start->GetResNum();
		        my $e = $end->GetResNum();
				print "Found prob helix $s to $e\n" if($verbose);
				push @indices, $start ; 
				push @indices, $end ; 
			}
		   last ;
		}
		foreach my $i (0..3){
			die if(!defined $res[$i]);
			push @RES, $res[$i];
		}
		my ($IsHelixPair,$distO1N2,$DIFFDIST) = IsHelixPair($pdb1,$RES[0],$RES[4]);
		if($IsHelixPair){
			if(@HELIXSTART){
		        $HELIXEND[0] = $RES[4]	
			 }
			 else{
			 	push @HELIXSTART, $RES[0];
		        push @HELIXEND, $RES[4]	
			 }
		}
		else{
			if(@HELIXSTART){

				my $start = $HELIXSTART[0];
				my $end = $HELIXEND[0];
		        my $s = $start->GetResNum();
		        my $e = $end->GetResNum();
				print "Found prob helix $s to $e\n" if($verbose);
				push @indices, $start ; 
				push @indices, $end ; 
			}
			@HELIXSTART = ();
			@HELIXEND = ();


		}
		my $single = $res->PrintSingleLetter($pdb1);

	}
	return @indices ;

}

sub IsHelixPair{
     my ($pdb1,@RES) = @_ ; 
		my @SINGLELETTER ;
		my @NUM ;
		my $NM = "";
		foreach my $r (@RES){
		    push @SINGLELETTER, $r->PrintSingleLetter($pdb1);
		    my $num = $r->GetResNum();
			$NM = $NM . "-" . $r->GetName() . "$num";
			push @NUM,$num ;
		}

		my @CA ; my @C ; my @N ; my @O ; my @CB ; my @RA ;
		foreach my $r (@RES){
		     my $CA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"CA");
		     my $C = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"C");
		     my $N = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"N");
		     my $O = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"O");
		     my $CB = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),"CB");
			 my $RAnm = $pdb1->GetReactiveAtomFromAtom($CA);
		     my $RA = $pdb1->GetAtomFromResidueAndType($r->GetResNum(),$RAnm);


		     push @CA, $CA ; push @C, $C ; push @N, $N ; push @O, $O ; push @CB, $CB ; push @RA, $RA ;
		}


        my $distO1N2 = util_format_float($pdb1->DistanceAtoms($O[0],$N[1]),1);
        my $distO2N1 = util_format_float($pdb1->DistanceAtoms($O[1],$N[0]),1);
        my $distCA = util_format_float($pdb1->DistanceAtoms($CA[0],$CA[1]),1);
        my $distOO = util_format_float($pdb1->DistanceAtoms($O[0],$O[1]),1);
        my $distNN = util_format_float($pdb1->DistanceAtoms($N[0],$N[1]),1);
        my $distCC = util_format_float($pdb1->DistanceAtoms($C[0],$C[1]),1);

       my $angle = geom_AngleBetweenTwoVectors_4Atoms($CA[0],$O[1],$CA[1], $O[1]);

		my $DIFFDIST = util_format_float(abs($distCA - $distOO),1);
		my $IsHelixPair = 1 ;
		if($distO1N2 > $HBONDDIST || $DIFFDIST > 1 || ($distO1N2 > $DISTO1N2 && $DIFFDIST > 0.5)){
			$IsHelixPair = 0 ; 
	    }
		print $ofhstats "$IsHelixPair   $distO1N2   $DIFFDIST   $distCA $distOO $protein $NM $distO2N1 $angle $distNN $distCC\n" if($verbose);
		return ($IsHelixPair,$distO1N2,$DIFFDIST);

}





sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}


