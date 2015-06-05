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
my ($codonfile,$fastafile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($aminoacid,$ignorefile,@expressions);
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "codonfile=s"=>\$codonfile ,
            "fastafile=s"=>\$fastafile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -codonfile ") if(!defined $codonfile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
usage( "Need to give a input file name => option -aminoacid ") if(!defined $aminoacid);




my ($info,$mapValues) = util_ReadCodonBiasFile($codonfile);


my $codontable = util_getCodonTable();



my ($retstr,$firstline) = util_readfasta($fastafile);

my @AAA ;
if($aminoacid){
	@AAA = ($retstr =~ /(.)/g);
}
else{
    my @ppp = ($retstr =~ /(...)/g);
	my $final = join "",@ppp ;
	die if ($final ne $retstr);
	foreach my $i (@ppp){
		die "Error: $i does not exist in codontable" if(! exists $codontable->{$i});
		my $aa = $codontable->{$i};
		push @AAA, $aa;
	 }
}

print $ofh ">junk\n";
foreach my $aa (@AAA){
	my $predtable = $info->{$aa};
    my $r = floor(1000*rand())+1;
	die "$r does not exist" if(!exists $info->{$aa}->{$r});
	my $v = $info->{$aa}->{$r};
	print $ofh "$v";
}
print $ofh "\n";



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
