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
my $ignorelist ;
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
            "ignorelist=s"=>\$ignorelist ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -listfile ") if(!defined $listfile);
my $ofh = util_write($outfile);

my $ignlist = {};
if(defined $ignorelist){
   my @list= util_read_list_sentences($ignorelist);
   map { s/\s*//g ; $ignlist->{$_} = 1 ; } @list ;
}

print STDERR "Info: parsing file $infile - might take some time\n";

my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);
print STDERR "Info: parsing file $listfile - might take some time\n";


my $finalstr = "";
foreach my $dbseq (keys %{$infoSeq2PDB}){
	$finalstr = $finalstr . $dbseq ;
}


#print $ofh $finalstr , "\n";
my ($infolist,$infoSeq2PDBlist,$mapChainedName2Namelist) = util_parsePDBSEQRESNEW($listfile,0);

print "Matcin\n";

	 my $DONE = 0 ;
	 my $NOTDONE = 0 ;
     foreach my $queryseq (keys %{$infoSeq2PDBlist}){
         my @l1 = @{$infoSeq2PDBlist->{$queryseq}};
	     my $nm1 = $l1[0];
		 next if(exists $ignlist->{$nm1});
         if($finalstr =~ /$queryseq/){
	 	    print $ofh " $nm1 \n";
	 	    print " $nm1 \n";
			$DONE++;
	     }
		 else{
		 	$NOTDONE++;
		 }
	    print "Done = $DONE, not done = $NOTDONE \n";
	 }



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
