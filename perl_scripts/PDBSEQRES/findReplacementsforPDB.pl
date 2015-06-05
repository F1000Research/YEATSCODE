#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR,$listfile,$hetpdb);
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
            "hetpdb=s"=>\$hetpdb ,
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
usage( "Need to give a output file name => option -hetpdb ") if(!defined $hetpdb);
usage( "Need to give a protein -option -protein  ") if(!defined $protein);
my $ofh = util_write($outfile);


print STDERR "Info: parsing file $infile - might take some time\n";

my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($infile,0);

foreach my $k (keys %{$mapChainedName2Name}){
		print $ofh "$k \n";
}

if(!exists $mapChainedName2Name->{$protein}){
	die "protein $protein does not exist - have you specified the chain?\n";
}

my $seq = $mapChainedName2Name->{$protein};

my @pdbswithseq  = @{$infoSeq2PDB->{$seq}} ;
my $N = @pdbswithseq ;
my $pdbstr = join ",", @pdbswithseq ;
print "There are $N proteins with the same seq as $protein $pdbstr \n";


my ($NOHET,$YESHET) = util_parseHETPDB($hetpdb);

foreach my $pdb (@pdbswithseq){
     if(exists $NOHET->{$pdb}){
	 	print "Replacement has $pdb has no het\n";
	 }
}


print STDERR "Output written in $outfile\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
