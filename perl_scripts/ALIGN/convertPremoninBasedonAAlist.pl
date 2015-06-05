#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
  use Time::HiRes qw( usleep ualarm gettimeofday tv_interval
   clock_gettime clock_getres  clock
   );

use PDB;
use Atom;
use Residue;
use ConfigPDB;


use POSIX qw(floor);
use Math::Combinatorics;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($pdb1,$pdb2,$infile,$outfile,$atomidx,$dontrunpymol);
my ($tag,$anndir,$listfile,$annotate,$anndir1,$anndir2,$aalist,$dist,$size,$findresidues,$maxresults,$inconf,$outconf,$resultfile,$checkself);
my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my $MINDIST = 2 ;
$, = "  ";
my $LEN = 10 ; 
my $verbose = 1 ;
GetOptions(
            "pdb1=s"=>\$pdb1 ,
            "anndir1=s"=>\$anndir1 ,
            "anndir2=s"=>\$anndir2 ,
            "listfile=s"=>\$listfile ,
            "checkself"=>\$checkself ,
            "dontrunpymol"=>\$dontrunpymol ,
            "findresidues"=>\$findresidues ,
            "atomidx=s"=>\$atomidx ,
            "tag=s"=>\$tag ,
            "resultfile=s"=>\$resultfile ,
            "maxresults=i"=>\$maxresults ,
            "dist=f"=>\$dist ,
            "len=i"=>\$LEN ,
            "size=i"=>\$size ,
            "infile=s"=>\$infile ,
            "aalist=s"=>\$aalist ,
            "grpconfig=s"=>\$grpconfig ,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a aalist => option -aalist ") if(!defined $aalist);
usage( "Need to give a anndir1 => option -anndir1 ") if(!defined $anndir1);
usage( "Need to give a anndir2 => option -anndir2 ") if(!defined $anndir2);
usage( "Need to give a size => option -size ") if(!defined $size);
usage( "Need to give a tag => option -tag ") if(!defined $tag);
usage( "Need to give a listfile => option -listfile ") if(!defined $listfile);
#usage( "Need to give a residue number") if(!defined $atomidx);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $ofh = util_write($outfile);

ConfigPDB_Init($grpconfig);

my ($aainfo) = util_ParseAAGroups($aalist);

my @list= util_read_list_sentences($listfile);

foreach my $i (@list){
   my $infile = "$anndir1/$i.$size.$tag.premon.in"; 
   if(! -e $infile){
   	  print "Warning: $infile does not exist\n";
	  next ;
   }
   my $ifh = util_read($infile);

	my @aainfolist ; 
	foreach my $IDX (1..4){
	     my $configtmp = "$anndir1/$i.$size.$tag.premon.in.config$IDX";
		 if(-e $configtmp){
		 	system ("cp $configtmp $anndir2");
		 }
	     my $aalisttmp = "$anndir1/$i.$size.$tag.premon.in.aalist$IDX";

		 my $aainfotmp = $aainfo ;
		 if(-e $aalisttmp){
               ($aainfotmp) = util_ParseAAGroups($aalisttmp);
			   	print "USing different AAgroups\n";
		 }
		 push @aainfolist, $aainfotmp;
	}

   print "iinfile = $infile\n" if($verbose);

   next if(! -e $infile || -z $infile);
   print $ofh "$i\n";
   my $Outfile = "$anndir2/$i.$size.$tag.premon.in"; 
   my $OFH = util_write($Outfile);
   while(<$ifh>){
		if(/^\s*STR /){
                my ($junk,$str) = split ;
				print "$str llll\n" if($verbose);
				my @l = split "", $str ;

				my @newl ;
				foreach my $id (1..4){
					push @newl, $l[$id -1] . $id ;
				}
				my $sets = {};
				my $IDX = 0 ;
				foreach my $K (@newl){
					my $aainfotmp = $aainfolist[$IDX];
					$IDX++;
					$sets->{$K} = [];

					my $KKK = shift @l;
				    my $ORIG = $aainfotmp->{$KKK} ; 
				    foreach my $grp (keys %{$aainfotmp}){
					    my $val = $aainfotmp->{$grp};
						if($val eq $ORIG){
					        print "$grp $val $IDX $K \n" if($verbose);
							push @{$sets->{$K}}, $grp ;
						}
				    }
				}

				my $A = shift @newl ;
				my $B = shift @newl ;
				my $C = shift @newl ;
				my $D = shift @newl ;
				#die " $A $B $C $D\n";
				my @A = @{$sets->{$A}};
				foreach my $a (@A){
				     my @B = @{$sets->{$B}};
				     foreach my $b (@B){
				          my @C = @{$sets->{$C}};
				          foreach my $c (@C){
				               my @D = @{$sets->{$D}};
				               foreach my $d (@D){
							   	    print $OFH "STR $a$b$c$d\n"; ;
							   }
						  }
				    }
				}
				die if(@l);


		}
		else{
			print $OFH $_ ; 
		}
  }
  close($ifh);
  close($OFH);
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;

}
