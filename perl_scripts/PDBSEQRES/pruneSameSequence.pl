#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR,$listfile);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "howmany=i"=>\$howmany ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a protein -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein -option -cutofflength  ") if(!defined $cutofflength);
my $ofh = util_write($outfile);

print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($infile,0);

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;


my $done = {};
my $total = 0 ;
my $final = 0 ;
foreach my $protein (@list){


    if(!exists $mapChainedName2Name->{$protein}){
	    warn "protein $protein does not exist - have you specified the chain?\n";
		next ;
    }
    my $seq = $mapChainedName2Name->{$protein};
	my $len = length($seq);
	$total++;

	next if(exists $done->{$seq});
	next if($len < $cutofflength);

	$final++;
	print $ofh "$protein\n";
	$done->{$seq} = 1;

}
print STDERR "Output written in $outfile. $final out of $total\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
