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
my ($infile,$outfile,$which_tech,$listfile,$protein,$what);
my (@expressions);
my $NCH = 2 ;
my $NLEN = 10 ;
my $HYDVALUE = 8 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "what=s"=>\$what ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "NCH=i"=>\$NCH ,
            "NLEN=i"=>\$NLEN ,
            "HYDVALUE=i"=>\$HYDVALUE ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a what file name => option -what ") if(!defined $what);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $info = {};
while(<$ifh>){
     next if(/^\s*$/);
     chop ;
	 my ($nm,$len,$hyd,$percent,$numcharge) = split ; 

	 next if($numcharge < $NCH);
	 next if($len < $NLEN);
	 next if($hyd < $HYDVALUE);

	 die "Give what as NEG or POS or HYD" if(!($what eq "POS" || $what eq "NEG" || $what eq "HYD"));
	 next if($what eq "POS" && $percent < 0.8);
	 next if($what eq "NEG" && $percent > 0.2);

	 $percent = 1 - $percent if($what eq "NEG");

	 my $val = $what eq "HYD"? $hyd: $hyd * $percent  + $numcharge*$percent ;
	 $info->{$_} = $val ;
}

foreach my $k (sort {$info->{$b} <=> $info->{$a}} keys %{$info}){
	my $v = $info->{$k} ;
	print $ofh "$k \n";
}
close($ifh);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
