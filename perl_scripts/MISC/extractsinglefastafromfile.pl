#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;


my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -which_tech ") if(!defined $which_tech);
my $ifh = util_read($infile);
my $info = {};
my $ofh ;
while(<$ifh>){
	 next if(/^\s*$/);
     if(/^\s*>/){
	 	if(defined $ofh){
           print "Extracted file $which_tech from $infile\n";
		   exit ;
		}
	 	close($ofh) if(defined $ofh);
		my @l = split ;
		my $str =  $l[0];
		$str =~ s/>//;
		if($str eq $which_tech){
            $ofh = util_write("$str.ALL.1.fasta");
		    print $ofh $_ ;
		}
	 }
	 else{
	    #die if(!defined $ofh);
	    #die if(/>/);
		if(defined $ofh){
	 	   print $ofh $_ ;
		}
	 }
}
close($ifh);


sub util_read{
     my ($outfile)= @_;
	 my $fh = new FileHandle($outfile,O_RDONLY) or die " could not read file $outfile as $!";
	return $fh ;
}

sub util_write{
     my ($outfile)= @_;
	  die "not defined" if(!defined $outfile);
	unlink $outfile;
	my $fh = new FileHandle($outfile,O_CREAT|O_WRONLY) or die " could not write file $outfile as $!" ;
	#print "Writing to $outfile\n";
	return $fh ;
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
