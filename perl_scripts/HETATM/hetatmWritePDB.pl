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
my ($ispolar,$infile,$outfile,$write,$pdb2het,$het2pdb,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
my $NUMOFRES = 4 ;
GetOptions(
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "het2pdb=s"=>\$het2pdb ,
            "pdb2het=s"=>\$pdb2het ,
            "listfile=s"=>\$listfile ,
            "hetatm=s"=>\$hetatm ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "write=i"=>\$write ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a protein pdb id -option -listfile  ") if(!defined $listfile);
usage( "Need to give a pdb2het pdb id -option -pdb2het  ") if(!defined $pdb2het);
usage( "Need to give a het2pdb pdb id -option -het2pdb  ") if(!defined $het2pdb);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my @ignorewhilewriting = qw (HOH);
my @ignorewhilereading = qw (CE FE CD O PO4 FE2 CU NA MN ACN MG HG ZN CA BME HOH ACT CL SO4);

my $ofhdata = util_append($pdb2het);
my $ofhhet2pdb = util_read($het2pdb);

my $done = {};
while(<$ofhhet2pdb>){
     next if(/^\s*$/);
	 if(/NOHET/ || /YESHET/){
	 	my ($junk,@l) = split ; 
		foreach my $i (@l){
		    $done->{$i} = 1 ;
		}
	 }
}
close($ofhhet2pdb);

my @list= util_read_list_sentences($listfile);
foreach my $protein (@list){
	next if(exists $done->{$protein});
    WriteHetData($protein);
}

sub WriteHetData{
   
   my ($protein) = @_ ;
   my $pdb = "$PDBDIR/$protein.pdb";
   if(!-e $pdb){
       print STDERR "$PDBDIR/$protein.pdb does not exist\n";
	   return ;
   }
   my $ifh = util_read($pdb);

   my $info = {};
   my $CHAIN ; 
   my $donehets = {};
   while(<$ifh>){
     next if(/^\s*$/);
	 last if(/ENDMDL/);
	 if(/^HETATM/){

	      my $LINE = $_ ;
	     my $len = length($LINE) ;
	  
	     my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($LINE);
		 my $key = "$atomnm , $alt_loc , $resname , $chainId , $resnum ,";
		 next if(exists $donehets->{$key});
		 $donehets->{$key} = 1 ;


	     $CHAIN = $chainId if(!defined $CHAIN);
	     if($CHAIN ne $chainId){
	       print "$protein Allow only one chain here -found $CHAIN and $chainId\n" ;
	       return ;
	     }
	  
	     $info->{$resname} = {} if(!defined $info->{$resname}); 
	     $info->{$resname}->{$resnum} = 0 if(! defined  $info->{$resname}->{$resnum});
	     $info->{$resname}->{$resnum} = $info->{$resname}->{$resnum} + 1 ;
	   }
  }
  close($ifh);

  my $ignoretable = util_make_table(\@ignorewhilewriting);
  print $ofhdata "$protein ";
  foreach my $hetatm (keys %{$info}){
  		next if(exists $ignoretable->{$hetatm});
  	    my $N = keys %{$info->{$hetatm}} ; 
		my $nAtoms ;
	    foreach my $res (keys %{$info->{$hetatm}}){
			$nAtoms = $info->{$hetatm}->{$res} if(!defined $nAtoms);
			my $oldn = $info->{$hetatm}->{$res} ;
			warn "Different numbers for $hetatm: $nAtoms and $oldn"  if ($nAtoms ne $oldn && $verbose);

		}
		print $ofhdata "$hetatm:$N:$nAtoms  "  ;
  }
  print  $ofhdata "\n";

}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
