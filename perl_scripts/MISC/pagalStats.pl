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
my ($infile,$outfile,$which_tech,$queryfile,$protein);
my (@expressions);
my $cutoff = 0 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "queryfile=s"=>\$queryfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my ($INFODB) = ParsePagalStats($infile);

if(defined $queryfile){
    my ($INFILE) = @_ ;
    my $ifh = util_read($queryfile);
    my $info = {};
    while(<$ifh>){
         next if(/^\s*$/);
         chop ;
	     my ($pdb,$nm,$D,$PD) = split ; 
	     next if(!defined $PD || $PD eq "nan" );
	     next if(!defined $D || $D eq "nan" );
		 my $MEANPD = $INFODB->{$nm}->{MEANPD};
		 my $MEAND = $INFODB->{$nm}->{MEAND};
		 my $SDD = $INFODB->{$nm}->{SDD};
		 my $SDPD = $INFODB->{$nm}->{SDPD};

		 my $diffD = util_format_float(abs($MEAND-$D));
		 my $diffPD = util_format_float(abs($MEANPD-$PD));
		 print "$diffD $SDD ---  $diffPD $SDPD ------------$MEAND $D , $MEANPD $PD \n";
	    
    }
	
}


my $info = $INFODB ;

my $what = "NUM";
my @resultssorted = sort { $info->{$a}->{$what} <=> $info->{$b}->{$what} } (keys %{$info});
foreach my $nm (@resultssorted){
   my $WHAT = $info->{$nm}->{$what}  ;
   my $N = $info->{$nm}->{NUM}  ;
		 my $MEANPD = $info->{$nm}->{MEANPD};
		 my $MEAND = $info->{$nm}->{MEAND};
		 my $SDD = $info->{$nm}->{SDD};
		 my $SDPD = $info->{$nm}->{SDPD};
   if($N > $cutoff){
         print $ofh  "$nm NUM=$N $MEAND $SDD $MEANPD $SDPD  \n";
   }
}


sub ParsePagalStats{
    my ($INFILE) = @_ ;
    my $ifh = util_read($INFILE);
    my $info = {};
    while(<$ifh>){
         next if(/^\s*$/);
         chop ;
	     my ($pdb,$nm,$D,$PD) = split ; 
	     next if(!defined $PD || $PD eq "nan" );
	     next if(!defined $D || $D eq "nan" );
	     if(!defined $info->{$nm}){
		    $info->{$nm} = {} ;
	        $info->{$nm}->{D} = [];
	        $info->{$nm}->{PD} = [];
	    }
	    push @{$info->{$nm}->{D}}, $D ;
	    push @{$info->{$nm}->{PD}}, $PD ;
    }
   foreach my $nm (keys %{$info}){
	 my @lD  = @{$info->{$nm}->{D}} ; 
	 my @lPD  = @{$info->{$nm}->{PD}} ; 
	 my $N = @lD ;
	 $info->{$nm}->{NUM} = $N ; 
	 my ($meanD,$sdD);
	 my ($meanPD,$sdPD);
	 if($N eq 1){
	      ($meanD,$sdD) = ($lD[0],0);
	      ($meanPD,$sdPD) = ($lPD[0],0);
	      $info->{$nm}->{SDD} = 0 ;
	      $info->{$nm}->{SDPD} = 0 ;
	 }
	 else{
	      ($meanD,$sdD) = util_GetMeanSD(\@lD);
	      ($meanPD,$sdPD) = util_GetMeanSD(\@lPD);

	      $info->{$nm}->{SDD} = $sdD ;
	      $info->{$nm}->{SDPD} = $sdPD ;
	 }
	 $info->{$nm}->{MEAND} = $meanD ;
	 $info->{$nm}->{MEANPD} = $meanPD ;
    }
    
    close($ifh);
    return $info; 
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
