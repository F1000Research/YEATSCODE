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
my $howmany = 100000 ;
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
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $info = {};
my $infoSubject2TRS = {};
my $infoStretch = {};
while(<$ifh>){
     next if(/^\s*$/);
     chop ;
	 my ($trs,$subject,$hasintrons,$qS,$qE,$qDIFF,$sS,$sE,$sDIFF)= split;
	 if(!$hasintrons){
          my $subjectcopy = $subject ; 
		  $subjectcopy =~ s/\..*//;

		  ### ensure we dont have both fwd and reverse
		  if(!exists $info->{$subjectcopy}){
		  	 $info->{$subjectcopy} = $subject ;
		  	 $infoSubject2TRS->{$subjectcopy} = $trs ;
		  	 $infoStretch->{$subjectcopy} = [];
		  }
		  else{
		  	  if($info->{$subjectcopy} ne $subject){
				    print "There is both fwd and reverse for $subjectcopy\n";
			   }
		  }

		  push @{$infoStretch->{$subjectcopy}}, $sS ;
		  push @{$infoStretch->{$subjectcopy}}, $sE ;
		 
		  push @{$infoStretch->{$subjectcopy}}, $qS ;
		  push @{$infoStretch->{$subjectcopy}}, $qE ;
	 }
}
close($ifh);


foreach my $subjectcopy (keys %{$infoStretch}){
	my @l = @{$infoStretch->{$subjectcopy}} ;
	my $subject = $info->{$subjectcopy};
	my $trs = $infoSubject2TRS->{$subjectcopy};

	my $strS = "";
	my $strQ = "";
	if(1){
	   ($strS) = util_readfasta($subject);
	   my $lenS = length($strS);
	   ($strQ) = util_readfasta("../../FASTADIR/$trs.ALL.1.fasta");
	   my $lenQ = length($strQ);
	}

	while(@l){
	    my $sS = shift @l	;
	    my $sE = shift @l	;
	    my $qS = shift @l	;
	    my $qE = shift @l	;

		my $totalSlen = length($strS);
        my $sliceFirst = util_extractSliceFromFasta($strS,1,$sS-1);
        my $sliceLast = util_extractSliceFromFasta($strS,$sE+1,$totalSlen);
        my $sliceS = util_extractSliceFromFasta($strS,$sS,$sE);
        my $sliceQ = util_extractSliceFromFasta($strQ,$qS,$qE);

		my $ofhnewfasta = util_write("NEWFASTA/$subjectcopy.ALL.1.fasta");
		print $ofhnewfasta ">$subjectcopy\n";
		print $ofhnewfasta "${sliceFirst}${sliceQ}${sliceLast}\n";
	    my $lenQ = length($sliceQ);
	    my $lenS = length($sliceS);
		print "$subject $sS $sE $qS $qE $lenQ $lenS  \n";
	}
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
