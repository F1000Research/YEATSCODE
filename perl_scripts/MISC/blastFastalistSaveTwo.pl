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
my ($config,$l1,$l2,$justone,$table,$save,$arg,$findmatch,$needleoutfile,@expressions);
my $simi ;
my $verbose = 0 ;
my $verb;
my $TEX2PDF = "/home/sandeepc/DATA/Tex2Pdf";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "l1=s"=>\$l1 ,
            "l2=s"=>\$l2 ,
            "justone=s"=>\$justone ,
            "outfile=s"=>\$outfile ,
            "needleoutfile=s"=>\$needleoutfile ,
            "config=s"=>\$config ,
            "verbose"=>\$verbose,
            "table=s"=>\$table,
            "save=s"=>\$save,
            "expr=s"=>\@expressions,
            "simi=f"=>\$simi ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a l1 -option -l1  ") if(!defined $l1);
usage( "Need to give a l2 -option -l2  ") if(!defined $l2);
usage( "Need to give a listfile -option -needleoutfile  ") if(!defined $needleoutfile);


my $ifhOUT = util_read($outfile);
my $donedeal = {};
while(<$ifhOUT>){
	my ($x,$y,$v) = split ;
	my $str = "$x.$y";
	$donedeal->{$str} = 1 ;
}
close($ifhOUT);
my $ofh = util_append($outfile);


$simi = 100 if(defined $table);
$verbose = 1 if(defined $verbose);

 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP) = util_SetEnvVars();

my @pdbs1= util_read_list_words($l1);
my @pdbs2= util_read_list_words($l2);

## check for files
{
my @pdbs = (@pdbs1,@pdbs2);

my $diefinally = 0 ;
foreach my $p (@pdbs){
    my $f = "$FASTADIR/$p.ALL.1.fasta";
	if(! -e $f){
	   print "Need to have fastafile... see output for list\n";
	   print $ofh "$f\n";
	   $diefinally = 1 ;
	}
}
my $N = @pdbs ;
print "Starting ... There are $N\n";
die if($diefinally);
}


my $CNT = 0 ;
while(@pdbs1){
	my $pdb1 = shift @pdbs1 ;
    my $f1 = "$FASTADIR/$pdb1.ALL.1.fasta";
	$CNT++;
	print "$CNT..";
    foreach my $pdb2 (@pdbs2){
	   my $str = "$pdb1.$pdb2";
	   next if(exists $donedeal->{$str});

       my $exec = "blastp -query FASTADIR/$pdb1.ALL.1.fasta -subject FASTADIR/$pdb2.ALL.1.fasta -out $needleoutfile";
	   system($exec);
	   my ($expect) = util_ParseBlastPW($needleoutfile);
	   print $ofh "$pdb1 $pdb2 $expect \n";
    }
}
print "\n";

sub util_ParseBlastPW{
	my ($infile) = @_ ;
	my $ifh = util_read($infile);
	while(<$ifh>){
		if(/Expect/){
			my ($val) = (/Expect = (.*),/);
			return $val ;
		}
	}
	close($ifh);
	return 1 ; 
}



chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
