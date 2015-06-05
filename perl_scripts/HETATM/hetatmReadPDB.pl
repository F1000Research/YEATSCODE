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
my ($ispolar,$pdb2het,$outfile,$het2pdb,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
my $NUMOFRES = 4 ;
GetOptions(
            "protein=s"=>\$protein ,
            "pdb2het=s"=>\$pdb2het ,
            "listfile=s"=>\$listfile ,
            "hetatm=s"=>\$hetatm ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "verbose=i"=>\$verbose ,
            "het2pdb=s"=>\$het2pdb ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a protein pdb id -option -het2pdb  ") if(!defined $het2pdb);
usage( "Need to give a protein pdb id -option -pdb2het  ") if(!defined $pdb2het);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my @ignorewhilewriting = qw (HOH);
my @ignorewhilereading = qw (CE FE CD O PO4 FE2 CU NA MN ACN MG HG ZN CA BME HOH ACT CL SO4);

my $ofhsizezone = util_write("size1");
my $donesizeone = {};

my $ignorewhilereading = util_make_table(\@ignorewhilereading);

ReadHetData();
sub ReadHetData{
   
   my $ifh = util_read($pdb2het);
   my $ofhcollate = util_write($het2pdb);

   my $info = {};
   my $sizeinfo = {};

   my $NOHET = {};
   my $YESHET = {};
   my $done = {};
   while(<$ifh>){
      my (@l) = split ;
	  my $protein = shift @l ;
	  my $N = @l ;
	  if(!$N){
	  	 $NOHET->{$protein} = 1 ;
	  }
	  else{
	      foreach my $x (@l){
	  	    my ($nm,$num,$size) = split ":", $x ;

			#if($nm eq "U1"){
				#print "$nm $num $size\n";
			#}


			if(exists $ignorewhilereading->{$nm} || $size eq 1){
				print "Ignoring $nm,$num,$size \n" if($verbose);
				if( ! exists $ignorewhilereading->{$nm} && ! exists $donesizeone->{$nm}){
				     $donesizeone->{$nm} = 1 ;
				     print $ofhsizezone "$nm\n";
				}
	  	        $NOHET->{$protein} = 1 ;
				next ;
			}

			if(! defined $YESHET->{$protein}){
	  	          $YESHET->{$protein} = 1 ;
			}
			else{
	  	          $YESHET->{$protein} =  $YESHET->{$protein} + 1 ;
			}

		    #print "$protein $nm $num $size \n";
		    if(!exists $info->{$nm}){
		        $info->{$nm} = [] ;
		        $sizeinfo->{$nm} = $size;
		    }
			if(! exists $done->{$nm.$protein}){
				$done->{$nm.$protein} = 1;
		        push @{$info->{$nm}}, $protein ;
			}
	      }

		  if(exists $YESHET->{$protein}){
		  	delete $NOHET->{$protein} ;
		  }
	  }
  }
  $, = " ";


  print $ofhcollate "NOHET ";
  foreach my $key (keys %{$NOHET}){
  	   print $ofhcollate "$key ";
  }
  print $ofhcollate "\n";

  print $ofhcollate "YESHET ";
  foreach my $key (keys %{$YESHET}){
  	   my $val = $YESHET->{$key};
  	   print $ofhcollate "$key:$val ";
  }
  print $ofhcollate "\n";


  my @sorted = sort { $sizeinfo->{$a} <=> $sizeinfo->{$b}} (keys %{$info});
  foreach my $key (@sorted){
  	 my @l = @{$info->{$key}};
	 print $ofhcollate "$key $sizeinfo->{$key} @l \n";
	
  }

}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
