#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($size,$DIR,$listfile,$ignorefile);
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
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
            "size=i"=>\$size ,
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
usage( "Need to give a protein -option -cutofflength  ") if(!defined $cutofflength);
usage( "Need to give a protein -option -size  ") if(!defined $size);
$outfile = $outfile . "$size";
my $ofh = util_write("MERGED/$outfile");
my $ofhmapping = util_write("MERGED/$outfile.mapping");

print "processing for size $size\n";
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);


my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}




my $done = {};
my $total = 0 ;
my $final = 0 ;

my $BEGIN = {};
my $END = {};

my $ENDSTR = "";
foreach my $seq (keys %{$infoSeq2PDB}){
	my $len = length($seq);
	die if(exists $done->{$seq});

	## first n and last n 
	my ($begin,$end) ;
	if($size eq 6){
	   ($begin) = ($seq =~/^(......)/);
	   ($end) = ($seq =~/(......)$/);
	}
	elsif($size eq 5){
	   ($begin) = ($seq =~/^(.....)/);
	   ($end) = ($seq =~/(.....)$/);
	}
	elsif($size eq 4){
	   ($begin) = ($seq =~/^(....)/);
	   ($end) = ($seq =~/(....)$/);
	}
	elsif($size eq 7){
	   ($begin) = ($seq =~/^(.......)/);
	   ($end) = ($seq =~/(.......)$/);
	}
	elsif($size eq 8){
	   ($begin) = ($seq =~/^(........)/);
	   ($end) = ($seq =~/(........)$/);
	}
	elsif($size eq 9){
	   ($begin) = ($seq =~/^(.........)/);
	   ($end) = ($seq =~/(.........)$/);
	}
	elsif($size eq 10){
	   ($begin) = ($seq =~/^(..........)/);
	   ($end) = ($seq =~/(..........)$/);
	}
	elsif($size eq 11){
	   ($begin) = ($seq =~/^(...........)/);
	   ($end) = ($seq =~/(...........)$/);
	}
	elsif($size eq 12){
	   ($begin) = ($seq =~/^(............)/);
	   ($end) = ($seq =~/(............)$/);
	}
	elsif($size eq 13){
	   ($begin) = ($seq =~/^(.............)/);
	   ($end) = ($seq =~/(.............)$/);
	}
	elsif($size eq 14){
	   ($begin) = ($seq =~/^(..............)/);
	   ($end) = ($seq =~/(..............)$/);
	}
	elsif($size eq 15){
	   ($begin) = ($seq =~/^(...............)/);
	   ($end) = ($seq =~/(...............)$/);
	}
	else{
		die "size can be only from 4 to 7";
	}
	$BEGIN->{$seq} = $begin;
	$END->{$seq} = $end;

	my @l = @{$infoSeq2PDB->{$seq}};
	die @l  if(@l>1);

	my $name = $l[0];

	$ENDSTR = $ENDSTR . $end . "X" . $name . "xx" ;

	my $N = @l ;
	$done->{$seq} = 1;

}

#print "$ENDSTR\n";
my $CNTfound = 0 ;

my $matches = {};
foreach my $b (keys %{$BEGIN}){
	my $A = $BEGIN->{$b};
	my @l1  = @{$infoSeq2PDB->{$b}} ;
	my $P1 = $l1[0];

	my (@l) = ($ENDSTR =~ /$A/);
	my $N = @l ;
	if(@l eq 1){
		my $CP = $ENDSTR ;
		$CP =~ s/.*xx$A/$A/;
		$CP =~ s/xx.*//;
		my ($x,$match) = ($CP =~ /.*($A).(.*)/);
		if($P1 ne $match){
		    print $ofh "$x $P1 $match\n";
		    $CNTfound++;
			$matches->{$P1} = $match ;
		}
	}
	next ;

	my $foundone = 0 ;
   foreach my $e (keys %{$END}){
	   my $B = $END->{$e};
   	   if($A eq $B){
	   	  my @l2  = @{$infoSeq2PDB->{$e}} ;
	   	  my $P2 = $l2[0];
		  next if($P1 eq $P2);

	   	  if(!$foundone){
		  	   print $ofh "$A $P1 ";
			   $foundone = 1 ;
		  }


		  print $ofh "$P2  ";
	   }
   }
   if($foundone){
   	  print $ofh "\n";
	  $CNTfound++;
   }
}



my $ofhchain = util_write("MERGED/$outfile.chain");
while((keys %{$matches})){
my $newmatches = {};
foreach my $m (keys %{$matches}){
	my $i = $matches->{$m} ;
	if(exists $matches->{$i}){
		my $j = $matches->{$i};
		print $ofhchain " chain $m $i $j \n";
		$newmatches->{$m} = $j;
	}
}
$matches = $newmatches ;
}


print "found matches on $CNTfound\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
