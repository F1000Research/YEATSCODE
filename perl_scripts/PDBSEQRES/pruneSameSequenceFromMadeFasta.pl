#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($writedata,$all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR,$length,$listfile,$ignorefile);
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
            "length=s"=>\$length ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
            "writedata=i"=>\$writedata ,
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
my $ofh = util_write($outfile);
my $ofhmapping = util_write("$outfile.mapping");

print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata);


my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}




my $done = {};
my $total = 0 ;
my $final = 0 ;
my $ofhlength ;
$ofhlength = util_open_or_append($length) if(defined $length);
foreach my $seq (keys %{$infoSeq2PDB}){
	my $len = length($seq);

	next if($len < $cutofflength);

	die if(exists $done->{$seq});

	my @l = @{$infoSeq2PDB->{$seq}};

	$total++;

	my $NN =@l ;
	my $CNTignored = 0 ;
	foreach my $l (@l){
		$l =~ s/.ORF.*//;
		$CNTignored++ if(exists $ignoretable->{$l});
	}

	if($CNTignored){
		die if($CNTignored ne $NN);
		next ;
	}

	#if(@l>1){
		foreach my $l (@l){
			print $ofhmapping "$l ";
		}
		print $ofhmapping "\n";

	#}


	my $protein = $l[0];

	$final++;
	print $ofh "$protein\n";
	if( defined $length){
	   print $ofhlength "$protein $len\n";
	}
	
	$done->{$seq} = 1;

}
print STDERR "Output written in $outfile. $final out of $total after ignorefile\n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
