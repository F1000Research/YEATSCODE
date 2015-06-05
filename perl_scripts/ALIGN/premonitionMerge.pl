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
my ($p1,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "p1=s"=>\$p1 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a protein 1 id -option -p1  ") if(!defined $p1);
#usage( "Need to give a dist -option -dist  ") if(!defined $dist);
#usage( "Need to give a mlength -option -mlength  ") if(!defined $mlength);
#usage( "Need to give a resnum -option -resnum  ") if(!defined $resnum);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my $file1 = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($file1);

my @res = $pdb1->GetResidues();
my $N = @res;
my $sum = 0 ;

my $DB = {};
my $DBDONE = {};
my $done = {};
my $indices = {};
my $first ;
foreach my $res (@res){
     next  if($res->GetAtomStr() ne "ATOM");
     my $resnum = $res->GetResNum(); 
     my $i = "$p1.$resnum.singleout";
     my $ifh = util_read($i) or die;
     print "Merging $i\n";
     while(<$ifh>){
	 	chomp; 
       next if(/^\s*$/);
       if(/^\s*#/ && /mlength/){
	   	   if(!defined $first){
		   	   $first = $_ ; 
		   }
		   else{
		   	  die if( $first ne $_);
		   }
	   }
	   my @l  = split ; 
	   my $nm = shift @l ;
	   foreach my $index (@l){
		    if(!exists $DB->{$nm}){
		        $DB->{$nm} = [];	
		    }
		    if(! exists $DBDONE->{$nm}->{$index}){
			    $DBDONE->{$nm}->{$index} = 1 ;
		        push @{$DB->{$nm}}, $index; 
		    }
	      }
	   }
   close($ifh);

}


print $ofh "# $first \n";

foreach my $k (keys %{$DB}){
	 my @list = @{$DB->{$k}} ; 
	 print $ofh "$k @list \n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
