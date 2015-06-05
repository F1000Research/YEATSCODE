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
my ($what,$filesgiven,$config,$p1,$p2,$justone,$table,$save,$arg,$findmatch,$needleoutfile,@expressions);
my $simi ;
my $verbose = 0 ;
my $verb;
my $postfix = "ALL.1.fasta";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "filesgiven"=>\$filesgiven ,
            "listfile=s"=>\$listfile ,
            "justone=s"=>\$justone ,
            "what=s"=>\$what ,
            "outfile=s"=>\$outfile ,
            "needleoutfile=s"=>\$needleoutfile ,
            "postfix=s"=>\$postfix ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "verbose"=>\$verbose,
            "table=s"=>\$table,
            "save=s"=>\$save,
            "expr=s"=>\@expressions,
            "simi=f"=>\$simi ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#$outfile = "tmp.blast" if(!defined $outfile);
usage( "Need to give a p1 -option -p1  ") if(!defined $p1);
usage( "Need to give a p2 -option -p2  ") if(!defined $p2);
usage( "Need to give a what -option -what  ") if(!defined $what);
usage( "Need to give a outfile -option -outfile  ") if(!defined $outfile);


my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my @listoffasta ;
push @listoffasta, $p1;
push @listoffasta, $p2;
my $diefinally = 0 ;
foreach my $p (@listoffasta){
    my $f = defined $filesgiven? $p: "$FASTADIR/$p.$postfix";
	if(! -e $f){
	   print "Need to have fastafile $f... see output for list\n";
	   $diefinally = 1 ;
	}
}
die if($diefinally);

my $a = defined $filesgiven? $p1 : $FASTADIR/$p1.$postfix ;
my $b = defined $filesgiven? $p2 : $FASTADIR/$p2.$postfix ;


die if($what ne "P" && $what ne "N");

my $whichblast = $what eq "P" ? "$SRC/BLAST/blastp":"$SRC/BLAST/blastn" ;
my $exec = "$whichblast -query $a -subject $b -out $outfile";
print "$exec\n";
system($exec);
my ($expect) = util_ParseBlastPW($outfile);
print "Wrote to $outfile. $expect is the first expectation value \n";
print "\n";




chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
