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
my ($infile,$outfile,$which_tech,$select,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "select"=>\$select ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $info = {};
my @L ; 
while(<$ifh>){
     next if(/^\s*$/);
     next if(/PD/);
	 s/.*D //;
     chop ;
	 my (@l) = split ; 
	 if(!defined $info->{1}){
		$info->{1} = \@l;
	}
	else{
		push @L, \@l;
	}
}

my @first = @{$info->{1}};
my $N = @first ;
my $Nsquare = $N * $N ; 
my $doit = {};
$doit->{1} = 1 ;
$doit->{2} = 1 ;
$doit->{4} = 1 ;
foreach my $l (@L){
	my @l = @{$l};
    my $sum = 0 ; 
	my $CNT = 0 ; 

	my $MAX = 0 ; 
	foreach my $X (@first){
		$CNT++;
		my $Y = shift @l ;
		next if(defined $select &&  !defined $doit->{$CNT});

		my $diff = abs($X - $Y) ;
		my $diffsquare = $diff * $diff ; 
		$sum = $sum + $diffsquare ;
		#print "$diff = diff \n";
		$MAX = $diff if($diff > $MAX);
	}
    my $rmsd = util_format_float(sqrt($sum)/ $N,1);
    print "RMDS = $rmsd max = $MAX \n";
}



chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
