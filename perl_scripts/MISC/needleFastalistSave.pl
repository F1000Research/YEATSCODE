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
my ($config,$justone,$table,$save,$arg,$findmatch,$needleoutfile,@expressions);
my $simi ;
my $verbose = 0 ;
my $verb;
my $TEX2PDF = "/home/sandeepc/DATA/Tex2Pdf";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "justone=s"=>\$justone ,
            "outfile=s"=>\$outfile ,
            "needleoutfile=s"=>\$needleoutfile ,
            "arg=s"=>\$arg ,
            "config=s"=>\$config ,
            "verbose"=>\$verbose,
            "table=s"=>\$table,
            "save=s"=>\$save,
            "expr=s"=>\@expressions,
            "simi=f"=>\$simi ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a listfile -option -needleoutfile  ") if(!defined $needleoutfile);
usage( "Need to give a arg -option -arg  ") if(!defined $arg);
my $CNT = 0 ; 


$simi = 100 if(defined $table);
$verbose = 1 if(defined $verbose);

 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP) = util_SetEnvVars();

my @pdbs= util_read_list_words($listfile);


my $diefinally = 0 ;
foreach my $p (@pdbs){
    my $f = "$FASTADIR/$p.ALL.1.fasta";
	if(! -e $f){
	   print "Need to have fastafile... see output for list\n";
	   print $ofh "$f\n";
	   $diefinally = 1 ;
	}
}
die if($diefinally);


my $N = @pdbs ;
print "Starting ... There are $N\n";
my $CNT = 0 ;
while(@pdbs){
	my $pdb1 = shift @pdbs ;
    my $f1 = "$FASTADIR/$pdb1.ALL.1.fasta";
	$CNT++;
	print "$CNT..";
    foreach my $pdb2 (@pdbs){
       my $f2 = "$FASTADIR/$pdb2.ALL.1.fasta";
	   my ($iden,$simival) = util_NeedleFiles($needleoutfile,$f1,$f2,$arg);
	   print $ofh "$pdb1 $pdb2 $iden $simival\n";
    }
}
print "\n";


chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
