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


usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
while(<$ifh>){
     next if(/^\s*$/);
     s/.*=N //;
	 my (@l) = split ; 
	 my $table = {};
	 my $realname = {};
	 foreach my $l (@l){
	    my ($P) = ($l =~ /C(\d+)_G/);
		if(! exists $table->{$P}){
			$table->{$P} = 1 ;
			$realname->{$l} = 1 ;
			#last ;
		}
	 }

	 if(keys %{$realname} > 1){
	 	foreach my $nm (keys %{$realname}){
			print "$nm ";
		}
		print "\n";
	 }

}
close($ifh);


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
