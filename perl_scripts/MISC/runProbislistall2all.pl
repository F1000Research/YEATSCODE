#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;
use ConfigPDB;
use MyGeom;
use BP;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile);
my ($config,$justone,$table,$save,$arg,$findmatch,$tmpfile,@expressions);
my $simi ;
my $verbose = 0 ;
my $verb;
my $postfix = "ALL.1.fasta";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "justone=s"=>\$justone ,
            "outfile=s"=>\$outfile ,
            "tmpfile=s"=>\$tmpfile ,
            "postfix=s"=>\$postfix ,
            "verbose"=>\$verbose,
            "table=s"=>\$table,
            "save=s"=>\$save,
            "expr=s"=>\@expressions,
            "simi=f"=>\$simi ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a listfile -option -tmpfile  ") if(!defined $tmpfile);


my $ifhOUT = util_read($outfile);
my $donedeal = {};
while(<$ifhOUT>){
	my ($x,$y,$v) = split ;
	my $str = "$x.$y";
	$donedeal->{$str} = 1 ;
}
close($ifhOUT);



 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP) = util_SetEnvVars();

my @pdbs= util_read_list_words($listfile);


my $diefinally = 0 ;
foreach my $p (@pdbs){
    my $f = "$PDBDIR/$p.pdb";
	if(! -e $f){
	   print "Need to have pdb... see output for list\n";
	   $diefinally = 1 ;
	}
}
die if($diefinally);


my $N = @pdbs ;
print "Starting ... There are $N\n";
my $CNT = 0 ;
while(@pdbs){
	my $pdb1 = shift @pdbs ;
	$CNT++;
	print "$CNT..";
    foreach my $pdb2 (@pdbs){
	   next if($pdb2 eq $pdb1);
	   my $str = "$pdb1.$pdb2";
	   next if(exists $donedeal->{$str});

       my $exec = "runprobis.csh $pdb1 $pdb2 $outfile\n";
	   system($exec);
    }
}


chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
