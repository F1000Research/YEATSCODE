#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$model,$listfile,$target);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "model=i"=>\$model ,
            "target=s"=>\$target ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a target pdb id -option -target  ") if(!defined $target);
usage( "Need to give a model pdb id -option -model  ") if(!defined $model);
my $ofh = util_write("$target.out");


my $distinfo = ParseDist($infile);
my $pdinfo = ParsePD($listfile);

my $scoreinfo = {};
my $max ; 
my @sorted = sort { $pdinfo->{$a} <=> $pdinfo->{$b} } (keys %{$pdinfo});
foreach my $k (@sorted){
	 my $val = $pdinfo->{$k} ; 
	 if(! defined $max){
			$max = $val ;
	 }
	 my $DVAL = $distinfo->{$k} ; 
	 if($DVAL > 0.011){
	     $scoreinfo->{$k} = 0 ; 
	 }
	 else{
	     $scoreinfo->{$k} = $max/$val ;
	 }
}


print $ofh "PFRMAT QA\n";
print $ofh "TARGET $target\n";
print $ofh "AUTHOR 3663-5056-2522\n";
print $ofh "METHOD Escapist and Proqud - http://f1000research.com/articles/2-243 and  http://f1000research.com/articles/2-211/v3\n";
print $ofh "MODEL $model\n";
print $ofh "QMODE 1\n";
@sorted = sort { $scoreinfo->{$b} <=> $scoreinfo->{$a} } (keys %{$pdinfo});
foreach my $k (@sorted){
	 my $v = util_format_float($scoreinfo->{$k},3) ; 
	print $ofh "$k $v \n";
}
print $ofh "END\n";

sub ParseDist{
	my ($Infile) = @_ ; 
    my $ifh = util_read($Infile);
	my $info = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     my (@l) = split ; 
		 my $nm = $l[0];
		 my $v = $l[3];
	    $info->{$nm} = $v ;
    }
    close($ifh);
	return $info ;
}

sub ParsePD{
	my ($Infile) = @_ ; 
    my $ifh = util_read($Infile);
	my $info = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     my (@l) = split ; 
		 my $nm = $l[21];
		 my $v = $l[5];
	    $info->{$nm} = $v ;
    }
    close($ifh);
	return $info ;
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
