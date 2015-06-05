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
my ($ispolar,$infile,$outfile,$doall,$initial,$which_tech,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 4 ;
my $maxdist = 10 ;
my $verbose = 0 ;

my $IGNOREBACKBONE = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "doall"=>\$doall ,
            "listfile=s"=>\$listfile ,
            "hetatm=s"=>\$hetatm ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "maxdist=i"=>\$maxdist ,
            "ispolar=i"=>\$ispolar ,
            "initial=f"=>\$initial ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
## for debug 
my $fhallinfo = util_write("$protein.allinfo");
usage( "Need to give a hetatm -option -hetatm  ") if(!defined $hetatm);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a ispolar pdb id -option -ispolar  ") if(!defined $ispolar);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

system("mkdir -p $hetatm");
my $DONEATOMS = {};

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my $newsvrpdb = "$hetatm/$protein.$hetatm.pdb";
my $ofhpdb = util_write($newsvrpdb);



my @NONPOLAR = qw(PRO ALA GLY ILE LEU MET VAL);
my $NONPOLAR = {};
foreach my $i (@NONPOLAR){
	$NONPOLAR->{$i}= 1;
}



my $log = "$hetatm/$protein.log";
my $ofhlog = util_write($log);

my $pseudoatom = new Atom();
$pseudoatom->SetIdx(10000);


my $RTouched = {} ;
my $RTouchedResAtom = {} ;
my $RTouchedDrugAtom = {} ;
my $WritePDB = 1 ; ## just write once

## initial rad and increment
my $R = defined $initial ? $initial : 2.5 ;
my $INCR = 0.1 ;

$, = " " ;
my ($N);
my (@t);
my ($found,$RESNUM);
while((keys %{$RTouched}) < $howmany){
	$R = $R + $INCR ;
	$N = keys %{$RTouched} ;
	@t = (keys %{$RTouched});
	print "rad = $R number found: $N Res numbers: @t ============ \n";
	($found,$RESNUM) = ProcessOneRadiiforHetatm($R,$WritePDB,$RTouched,$RTouchedResAtom,$RTouchedDrugAtom,$hetatm);
	if(!$found){
		print STDERR "Did not find hetatm $hetatm \n";
		die ;
	}
	$WritePDB = 0 ; 
	last if($R > $maxdist);
}

	$N = keys %{$RTouched} ;
	@t = (keys %{$RTouched});
print "rad = $R number found: $N Res numbers: @t ============ \n";

my $N = keys %{$RTouched} ;
print "Found $N residues in $R radius\n";


my $claspout = "$protein.$howmany.clasp.in";
my $fhclasp = util_write($claspout);
my $table = "$protein.$howmany.$ispolar.table.in";
my $contact = "$protein.$howmany.$ispolar.contact.in";
my $fhtable = util_write($table);
my $fhcontact = util_write($contact);
my $cnt = 0 ;
foreach my $key (sort {$RTouched->{$a} <=> $RTouched->{$b}} keys %{$RTouched}){
	   $cnt++; 
       last if($cnt > $howmany);
	   my ($res) = $pdb1->GetResidueIdx($key);
	   my $v = $RTouched->{$key};
	   $v =~ s/\s//g;
	   my $atomnm = $RTouchedDrugAtom->{$key};
	   my $type = $RTouchedResAtom->{$key};
	   my $nm = $res->GetName() . "$key";;
	   print $fhclasp "$nm ";
	   print $fhtable "$nm/$type/$atomnm/$v ";
	   print $fhcontact "$atomnm $RESNUM\n";
}
print $fhclasp "\n";
print $fhtable "\n";
print $fhcontact "\n";
print "Wrote in $claspout \n";

##################################################




sub ProcessOneRadiiforHetatm{
   my ($radii,$writepdb,$residuestouched,$residuestouchedResatom,$residuestouchedDrugatom,$hetatmnm) = @_ ;
   my $datafile = "$hetatmnm/$hetatmnm.$radii.dat";
   system("touch $datafile") if(!-e $datafile); 
   my $ofhdata = util_append($datafile);
   
   my @backboneatoms = qw( O N CA);
   my $backboneatomstable = util_make_table(\@backboneatoms);
   my $ifh = util_read($pdb);

   my $info = {};
   my $ONERES ; 
   my $foundhetam = 0 ; 
   while(<$ifh>){
     next if(/^\s*$/);
	 if(/^HETATM/ ){

	  my $LINE = $_ ;
	  my $len = length($LINE) ;
	  
	  my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($LINE);


	  if($resname =~ /$hetatmnm/i){
	  	  $foundhetam = 1 ;

		  # ensure only one bound molecule is processed
	      $ONERES = $resnum  if(!defined $ONERES);
	      next if($resnum ne $ONERES);

	  	  print $ofhpdb $_ if($writepdb);

	      $pseudoatom->SetCoords($x,$y,$z);

          my $list = util_make_list($pseudoatom);
	      my ($junk,$neigh)  = $pdb1->GetNeighbourHoodAtom($list,$radii);


		   my $sort = {} ;

		   ## Hydrogen bonds
		   my $Allowed = {};
		   $Allowed->{"ON"} = 1 ; $Allowed->{"NO"} = 1 ;
		   $Allowed->{"OH"} = 1 ; $Allowed->{"HO"} = 1 ;
		   $Allowed->{"OO"} = 1 ;
		   $Allowed->{"NN"} = 1 ;
		   $Allowed->{"NH"} = 1 ; $Allowed->{"HN"} = 1 ;
		   $Allowed->{"SH"} = 1 ; $Allowed->{"HS"} = 1 ;

		   my $sortedneigh = {};
		   foreach my $a (@{$neigh}){

		       my $d = $pdb1->DistanceAtoms($pseudoatom,$a);
			   $sortedneigh->{$a} = $d ;
		   }
		   my @sorted = sort { $sortedneigh->{$a} <=> $sortedneigh->{$b} } (@{$neigh});
		   foreach my $a (@sorted){
			   #$a->Print();

			   my $NNMM = $a->GetName();

		       my $d = $pdb1->DistanceAtoms($pseudoatom,$a);

			   #print "$NNMM $atomnm $d lllllllllll\n";

			   my $num = $a->GetResNum();
			   my $type = $a->GetType();
	           my ($res) = $pdb1->GetResidueIdx($num);
			   my $resnm = $res->GetName();


	           if($IGNOREBACKBONE && exists $backboneatomstable->{$type}){
			      next ;
			   }


			   my @typel = split "", $type ;
			   $atomnm =~ s/\s//g;
			   my ($atomfirst) = ($atomnm =~ /(.)/);
			   my ($typefirst) = ($type =~ /(.)/);
			   my $X = $atomfirst . $typefirst ;
				
			   ## for debug 
		   	   if(! exists $DONEATOMS->{$a}){
			       print $fhallinfo "$resnm$num$type $atomnm $X $d\n";
			   }
				$DONEATOMS->{$a} = 1 ;

			   if(!defined $doall){
			       if($ispolar){
			           next if(exists $NONPOLAR->{$resnm});
			           next if(!exists $Allowed->{$X});
			       }
			       else{
			           #next if(!exists $NONPOLAR->{$resnm});
			           next if(exists $Allowed->{$X});
			       }
			   }

	                    
			   if(! exists $residuestouched->{$num}){
			      $residuestouched->{$num} = util_format_float($d,1) ;
			      $residuestouchedResatom->{$num} = $type ;
			      $residuestouchedDrugatom->{$num} = $atomnm ;
			   }
			   my $nm = $a->GetName();
			   $a->Print() if($verbose);
			   print $ofhdata "$atomnm $protein $resnm $num $type $d \n";
			   print $ofhlog "$radii $atomnm $protein $resnm $num $type $d \n";
		   }
	  } ## matches the hetatm
	}
  }

  close($ifh);
  return ($foundhetam,$ONERES) ;
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
