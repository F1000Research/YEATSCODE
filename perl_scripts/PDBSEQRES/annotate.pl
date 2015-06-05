#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyWeb;
use PDB ;
use ConfigPDB;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($promfile,$html,$all,$infile,$outfile,$scores,$or,$silent,$groupinfo,$pdb);
my ($addto,$matches,$DIR,$listfile,$query2ref,$mapping);
my $howmany = 600000 ; 
my ($WEIGHTEC) ;
my $threshhold = 2 ; 
my $cutofflength = 0 ; 
my $cutoffscore = 10000 ; 
my $isdummy = 0 ; 
my @types = (); 
my @ntypes = (); 
my @motifs = (); 
my $caption = "XXXXXXXXXXXXXXXXXXXXX";
my $header1 = "XXXXXXXXXXXXXXXXXXXXX";
my $header2 = "";
my $title = "CLASP Database";
my $ANNFILEDIST = 0 ;
GetOptions(
            "all"=>\$all ,
            "query2ref"=>\$query2ref ,
            "groupinfo"=>\$groupinfo ,
            "scores"=>\$scores ,
            "silent"=>\$silent ,
            "mapping"=>\$mapping ,
            "infile=s"=>\$infile ,
            "promidx=s"=>\$promfile ,
            "pdb=s"=>\$pdb ,
            "title=s"=>\$title ,

            "addto=s"=>\$addto ,
            "header1=s"=>\$header1 ,
            "header2=s"=>\$header2 ,

            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "howmany=i"=>\$howmany ,
            "matches=s"=>\$matches ,
            "anndist=i"=>\$ANNFILEDIST ,
            "isdummy=i"=>\$isdummy ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "cutoffscore=f"=>\$cutoffscore ,
            "type=s"=>\@types,
            "caption=s"=>\$caption,
            "ntype=s"=>\@ntypes,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);

my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my ($ANNDIR) = "ANNOTATE";
my ($MAPPINGFILE) = $ENV{MAPPINGFILE} or die ;
my ($infoMapping,$uniqueEC,$uniqueSP) = util_read_Mapping_PDB_2_SWISSPROT($MAPPINGFILE);
ConfigPDB_Init($grpconfig);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info) = util_parsePDBSEQRES($infile,0);



my $activesitesoverlap = {};
my $ifhlist = util_read($listfile);
my @listoffiles ; 
my @origPDB ; 
while(<$ifhlist>){
     next if(/^\s*$/);
     chop ;
	 my ($nm,$fname) = split ; 
	 push @listoffiles, $fname ; 
	 push @origPDB, $nm ; 
}

foreach my $LIST (@listoffiles) {


   my $origPDB = shift @origPDB ; 

   $outfile = $LIST . ".annotated"  ;
   my $ofh = util_write($outfile);
   
   my @list= util_read_list_sentences($LIST);
   my $list = {};

   if(defined $scores){
   	   my @l ;
       foreach my $i (@list) { 
   	     my ($name,$score) = split " ",$i   ;  
   		 next if($score > $cutoffscore);
   	     $list->{lc($name)} = $score ; 
   		 push @l, $name ; 
   	}; 
   
   	@list = @l ;
   }
   else{
       map { s/\s*//g ;  $list->{lc($_)} = 1 ; } @list ;
   }
   
   
   	my @p ;
   	my $doubleslashes = "";
   	my @tableheaders ;
   	my $cellentries = {} ;
   	$doubleslashes = "\\\\";
   	util_printTablePre($ofh,$caption);
   
   
   	my $cnt = 0 ;
   	my $nm2Orig = {} ;
   
   
   	## iterate over @list - the list of proteins
       map { 
   	    my $nm = uc($_);
		next if($nm eq $origPDB);
   	    my $orignmnm = uc($_);
   	    my $lc = lc($_);
   		my $ec ;
           if(defined $mapping){
               #my ($infoMapping,$uniqueEC,$uniqueSP) = util_read_Mapping_PDB_2_SWISSPROT($MAPPINGFILE);
               my $x = util_getECfromPDB($infoMapping,$lc);
   	        if(defined $x && @{$x} > 0){
   		        $ec = $x->[0];
   		        #print "WWWWWWW $lc $ec \n";
   	        }
           }
   
   
   		my ($tableResidues,$nresidues,$juuunk) ;
   		my $EXISTSINACTIVESITE = "No" ;
   		my $h = $info->{$lc} ; 
   		if(!defined $h){
   		    warn "did not fine $lc ";
   	    }
   		else{
   			($tableResidues,$nresidues,$juuunk) =  Config_ReadAnnfile($nm,$ANNDIR,0);
   		    $EXISTSINACTIVESITE = "Dont know active residues" if($nresidues == 0);
   			#print "$EXISTSINACTIVESITE \n";
   			my ($closetoactivesitePDB,$dist) =  Config_ReadAnnfileClose($pdb,$ANNDIR,0);
   			$ANNFILEDIST = $dist if($ANNFILEDIST eq 0);
   
   
   			###########    Get the first entry #############
   
   
   			if($h->{LEN} > $cutofflength && $cnt < $howmany ){
   				if(defined $scores){
   				    my $score = $list->{$lc};
   					$ec = "X" if(!defined $ec);
   		            print $ofh "$nm,$ec,$h->{LEN} &  $h->{FULLNM} & $score & $nresidues $doubleslashes \n";
   				}
   				else{
   		                  print $ofh "$nm &  $h->{LEN} &  $h->{FULLNM} &  $nresidues $doubleslashes \n";
   				}
   				$cnt++;
   			}
   		}
   	} @list ;
   
   
   
   	my $lc = lc($caption);
       my $h = $info->{$lc} ; 
   	 $h->{FULLNM}  =  "HHH" if(! defined  $h->{FULLNM} );
   	my $newcaption = $caption . " : " . $h->{FULLNM} ; 
   	util_printTablePost($ofh,$newcaption);
   
   
   if(defined $scores && defined $query2ref && defined $promfile){
   		my $lc = lc($pdb);
   		my $uc = uc($pdb);
   		my $FFFHHH = util_write("$uc.matches");
   		my $h = $info->{$lc} ; 
   		my $LEN = $h->{LEN};
   	my $promIndex = 0; 
   	my $moonIndex = 0; 
       #my ($infoMapping,$uniqueEC,$uniqueSP) = util_read_Mapping_PDB_2_SWISSPROT($MAPPINGFILE);
       my $x = util_getECfromPDB($infoMapping,$pdb);
   	if(defined $x && @{$x} > 0){
   		my $ec = $x->[0];
   		my $CNTOFMATCHES    = 0 ;
   		my $CNTOFNOTMATCHES = 0 ;
           foreach my $c (0..$cnt-1){
   	        my @l = @{$cellentries->{$c}} ;
   		    my $N = @l - 1 ; 
   		    my $score = $l[$N -1];
   		    my $last = $l[$N];
   		    my $nm = $l[0];
   		    my $orignmnm = $nm2Orig->{$nm} or die ; 
   		    next if($pdb eq $orignmnm);
   		    if($score < $threshhold){
                   my $sc = util_EC_CorrelatePDBS($infoMapping,$pdb,$orignmnm,$WEIGHTEC);
   			    if(defined $sc){
   			       $sc = $sc/$score ;
   		           if(!defined $activesitesoverlap->{$pdb.$orignmnm}){
   			           $moonIndex = $moonIndex + $sc ;
   				       $CNTOFNOTMATCHES++;
   		           }
   				   else{
   				       print " OVERKAP $pdb.$orignmnm  $activesitesoverlap->{$pdb.$orignmnm} \n";
   					   print $FFFHHH "$orignmnm\n";
   			           $promIndex = $promIndex + $sc ;
   					   $CNTOFMATCHES++;
   					}
   			    }
   		    }
           }
   		my $TOTALMATCHES = $CNTOFMATCHES + $CNTOFNOTMATCHES ;
   		my $append = util_append($promfile);
   		die if(!defined $append);
   	    print $append "$threshhold, $ANNFILEDIST, $promIndex, $moonIndex, $TOTALMATCHES , $CNTOFMATCHES , $CNTOFNOTMATCHES, $cnt ,$LEN  PROM score for $pdb with EC $ec  \n";
   
   	   }
   	   else{
   		   print "Could not get EC number for $pdb \n";
   	   }
   }



   print STDERR "Output written in $outfile\n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}



