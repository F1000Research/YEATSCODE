#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Tech ;
use Actel ;
use MyUtils;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$samefile,$outfile,$which_str,$ignore,$with_what);
GetOptions(
            "which_str=s"=>\$which_str ,
            "ignore=s"=>\$ignore ,
            "with_what=s"=>\$with_what ,
            "infile=s"=>\$infile ,
            "outfile=s"=>\$outfile ,
            "samefile"=>\$samefile ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a with_what => option -with_what ") if(!defined $with_what);
usage( "Need to give a which_str => option -which_str ") if(!defined $which_str);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
unlink $outfile ;
my $ofh = new FileHandle($outfile,O_CREAT|O_WRONLY) or die ;
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = new FileHandle($infile,O_RDONLY) or die "Could not read $infile";
print STDERR "Writing to file $outfile\n";
while(<$ifh>){
     next if(/^\s*$/);
     next if(defined $ignore && /$ignore/);
	 ## matches exact 
     ### s/\b$which_str\b/$with_what/g;

	 ## matches anywhere 
     s/$which_str/$with_what/g;
     print  $ofh $_ ;
}
close($ifh);
chmod 0777, $outfile ;
if(defined $samefile){
    util_printAndDo("cp -f $outfile $infile");
}
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
