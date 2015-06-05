#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($trs,$forR,$reverse,$infile,$cutoff,$outfile,$which_tech,$seperator,$listfile);
my (@idx);
$seperator = ",";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "cutoff=s"=>\$cutoff ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "forR=s"=>\$forR ,
            "seperator=s"=>\$seperator ,
            "idx=i"=>\@idx,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);

my $ofh = util_write($outfile);

my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();
my $ifh = util_read($infile);

my $rettable = util_R_extractSingleLine($ifh,"X");

my $CNTignored = 0 ;
foreach my $trs (keys %{$rettable}){
	my @X = @{$rettable->{$trs}} ;
	my ($max) = util_GetMeanSD(\@X);
   if($max < 100){
   	  $CNTignored++;
	  next ;
   }

   die "$trs not found " if(!@X);
   
   my $XSTR = "X = c ( " ;
   my $seperator = "";
   my $first = 1 ;
   while(@X){
	   my $x = shift @X ;
	   if($first){
		   $seperator = "";
		   $first = 0;
	   }
	   else{
		   $seperator = ",";
	   }
	   $XSTR = $XSTR .  "$seperator $x";
   }

   print $ofh "$XSTR ) \n";
   print $ofh "print (\"$trs\") \n";
   print $ofh "myboxplot <- boxplot(X) \n";
   print $ofh "myboxplot\$out \# it will print the values of the outliers  \n";
   print $ofh "quantile(X,.99) \n";

}

print "ifnore $CNTignored\n";

sub util_R_extractSingleLine{
   my ($ifh,$varname) = @_ ;
   my $info = {};
   while(<$ifh>){
	 	my @l = split ;
		my $id = shift @l ;
		$info->{$id} = \@l;
   }
   return $info ;
}


print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
