#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions,$config);
my $aminoacid = 1;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "config=s"=>\$config ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_append($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a config file name => option -config ") if(!defined $config);

my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


ConfigPDB_Init($config);
my $pdb1 = new PDB();


my ($string,$firstline) = util_readfasta($infile);


  my @ALL = ($string =~ /./g);
  my $LEN = @ALL;
  #my @l = qw (Asp Glu Gly Ala Pro Leu Val );
  my @l = qw(Gly Pro Ala Val Leu Ile Met Phe Tyr Trp His Lys Arg Gln Asn Glu Asp Cys Ser Thr) ;
  if(!$aminoacid){
      @l = qw(Gly Ala Cys Thr) ;
  }
  #y @l = qw(111 222 333 444 555 666 777 888 999  10  11  12  13  14  15  16  17  18  19  20) ;
  my $sum = 0; 
  print $ofh "$protein $LEN ";

  my @lll ;

  my $table = {};
  $table->{"P"} = 1 ;
  $table->{"K"} = 1 ;
  $table->{"V"} = 1 ;

  my $cntamino = 0 ;
  foreach my $i (@l){
  	my $one = $pdb1->GetSingleLetter($i);
    my @M = ($string =~ /$one/g);
    my $N = @M ; 
	if(exists $table->{$one}){
	   $cntamino = $cntamino+$N;	
	}
	$sum = $sum + $N ; 
	my $per = int(100*($N/$LEN));
	print $ofh " $per ";
	push @lll, $per ;
  } 
  if($aminoacid){
     print $ofh "\n";
  }
  else{
     my ($mean,$sd) = util_GetMeanSD(\@lll);
	 my $diff = ($lll[0] + $lll[2]) - ($lll[1] + $lll[3]) ;
     print  $ofh " $sd $diff \n";
  }

  if($aminoacid && !$cntamino){
  	    die "Something wrong - maybe you are using nt for amino acids\n";
  }

  #my $percent = (100 * $sum )/$LEN ; 
  #print $ofh "$protein $percent\n";

  # do stuff with $string

sub usage{
    my ($msg) = @_ ;
	    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
		    die ;
			}
