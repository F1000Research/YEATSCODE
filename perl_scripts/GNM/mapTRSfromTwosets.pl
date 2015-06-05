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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cutoff=f"=>\$cutoff ,
            "verbose=i"=>\$verbose ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}


my $info = {};

my $map = {};
my $maptwomatches = {};
my $mapall = {};
while(<$ifh>){
     next if(/^\s*$/);
     chop ;
	 my (@l) = split ; 
	 my $N = @l ;
	 my $trs1 = $l[0];
	 my $A = $trs1 ;
	 $A =~ s/.ORF.*//;
	 $mapall->{$A} = 1 ;


	 if($N eq 1){
	 }
	 else{
	     my $eval = $l[$N-1];
		 my $trs2 = $l[1];
		 #if($eval < 0.00000000000000001){
		 if($eval < $cutoff){
			 my $B = $trs2 ;
			 $B =~ s/.ORF.*//;
			 if(exists ($map->{$A})){
			 	my @l = @{$map->{$A}};
			 	my $old = $l[0];

				my $oldprefix = $old ;
				my $PPP = $B;
				($PPP =~ s/_.*//);
				($oldprefix =~ s/_.*//);

				if($PPP ne $oldprefix){
					$maptwomatches->{$A} = 1 ;
					
				}
			 }
			 else{
			 	$map->{$A} = [];
			 	push @{$map->{$A}}, $B ; 
			 }
		 }
		 else{
		 }
	 }

}


my $TRS2new = {};

my $CNTmappeduniquely = 0;
my $CNTtwomaps = 0;
foreach my $k (keys %{$map}){
		my @l = @{$map->{$k}} ;
		my $N = @l ;

		foreach my $v (@l){
		    $TRS2new->{$v} = [] if(!defined $TRS2new);
		    push @{$TRS2new->{$v}}, $k;
		}

	if(! exists $maptwomatches->{$k}){
		$CNTmappeduniquely++;
	}
	else{
		$CNTtwomaps++;
	}
}


foreach my $k (keys %{$TRS2new}){
	my @l = @{$TRS2new->{$k}};
	print $ofh "$k ";
	foreach my $l (@l){
	     print $ofh "$l ";
	}
	print $ofh "\n";
}
my $CNTnotmapped = 0 ;
my $ofhnomapped = util_write("notmapped");
foreach my $k (keys %{$mapall}){
	if(! exists $map->{$k}){
		$CNTnotmapped++;
		print $ofhnomapped "$k\n";
	}
}
print "$CNTnotmapped CNTnotmapped $CNTmappeduniquely CNTmappeduniquely $CNTtwomaps CNTtwomaps\n";


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
