
package MyUtils;
use Carp ;
use POSIX ;
require Exporter;
use Algorithm::Combinatorics qw(combinations) ;
use Math::NumberCruncher;
use Math::MatrixReal;  # Not required for pure vector math as above
use Math::Geometry ; 
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
#use Math::Trig;
#use Math::Trig ':radial';
no warnings 'redefine';
my $EPSILON = 0.01;
use MyPymol;

my $verbose = 0 ;

  local $SIG{__WARN__} = sub {};

@ISA = qw(Exporter);
@EXPORT = qw( util_execute_chop util_exit_if_doesnt_exist util_get_tag util_print_n_log util_get_user 
            util_get_cmdline util_num_lines util_percentage util_percentages util_printAndDo util_get_pwd util_log_base util_check_machine
	    util_get_tech util_printAndWrite2Script 
        util_get_time_from_string util_isDesignSuccessful util_parseTimelog util_diff util_round
		parseFPocketFile 
        util_print_mgc_home util_get_mgc_home util_write_list_in_file util_fullline 
		util_ceil util_floor
		util_getComplimentaryString
		util_get_max 
		util_findLength

util_R_extractIndex
		util_extractSliceFromFastaString
		util_extractSliceFromFasta
		util_open_or_append
		util_Blastout_isrev
		util_parsePDBSEQRESNEW
		util_ConcatTwoStringsWithCommonBeginAndEnd
		util_ReadPremonOut
		util_GetDisulphide
		util_ParsePremonIn
		util_IsPro util_IsMSE 
		util_ReadCodonBiasFile
		util_getCodonTable
		util_helixWheel
		util_ReadPremonOutAndGiveFasta
		util_IsResisdueThisType
        util_verilog2eqn util_runSis util_two2thepowerof util_assign_list2list util_emitGrid
        util_wait_on_lockfiles util_writePsp util_get_grid_cmd util_print_synp_file
		util_writelist2file
		util_readAPBSPotentialFromStart
		util_sortsingleString
		util_readfasta
		util_ReadAnnotateFile util_ReadAnnotateFileFRAGAL
		util_ParseBlastPW utils_parseBlastOut util_ParseWebBlast util_PARSEBLAST
		util_parseHolo2Apo
		vprint vprintheader vSetVerbose vIncrVerbose vDecrVerbose
        util_copy_table
		util_GetClosestAtoms_intwoPDBs
		util_ParseDSSP
		util_Ann2Simple
		util_AlignAndMatchRemaining
		util_GetMeanSD
		parseCastp
		util_AlignAndMatchRemainingAtoms
util_AreResiduesContinuous
util_FindRmsd
util_FindRmsdAllAtoms
util_ProcessSingleLine
util_parseHETPDB
util_ProcessSingleLineAtomlist

util_RunHelixWheel
util_ParseAAGroups
util_GetPotForAtom
util_GetPotDiffForAtoms
util_GetPotDiffForResidues
util_ReadPdbs
util_WriteClustalAln
		util_WriteFastaFromAtoms
		util_WriteFastaFromResidueNumbers
		
        util_print_vars util_uniq2
        util_SAVEDIR util_INITNAMES
        util_read util_append util_write util_pick_random_from_list
		util_make_list util_make_table
		util_GetFasta
		util_wget util_get_pdb 
		util_makeCSH
		util_parse_pdbseqres
		util_readPDB
		util_read_list_words
		util_read_list_sentences
		util_enter_maxcutoff
        util_is_integer util_is_float util_EnterName util_EnterNumber util_EnterTwoNumbers util_EnterTwoNames 
		util_AddBeforeEach
		util_ignoreIfDistanceIsLessthan util_ignoreIfAnyAreEqual util_ignoreIfAnyAreFurtherThanCutoff
		util_split_list_numbers
		util_read_Mapping_PDB_2_SWISSPROT util_filter_basedon_EC util_getPDBID_basedon_SP util_getECfromPDB
		util_ReadLine
		util_ParseBlocks util_ParseBlockForString
		util_ProcessRowAndColumnsForMean
		util_SetEnvVars
		util_ScoreDistance util_ScorePD
		util_format_float
		util_mysubstrDontStripSpace util_mysubstr
		util_printResult
		util_printTablePre util_printTablePost util_printTableLine util_PrintMeanAndSD
		util_Annotate util_ReadAnnfile
		util_table_print
		util_pick_n_random_from_list
		util_ExtractSliceFromFasta
		util_GetFastaFiles 
		util_maketablefromfile 

util_NeedleSeq
util_NeedleFiles
util_NeedlePDBNamesFromSeq
util_NeedlePDBNamesFromFASTADIR
		util_readAPBSPotential
		util_usage util_CmdLine
		util_parsePDBSEQRES
		util_SortTwoStrings
		util_GetEnv
		util_IsZero
		util_Banner util_PrintInfo

		util_getTmpFile
		util_PrintOutConf

ParseAPBSResult
		util_printHtmlHeader util_printHtmlEnd util_HtmlizeLine
		util_HtmlTableHead util_HtmlTableEnd util_HtmlTableCell
		util_MakeLink
		util_EC_CreateLevels util_EC_AddPDB util_EC_CorrelatePDBS
		util_PairwiseDiff
util_readPeptideInfo
util_getECID

util_ExtractOneIdxFromFile
	    );

use strict ;
use FileHandle ;
use Getopt::Long;


my $havetokeepthispostive = 13 ;

sub util_SetEnvVars{
   my @vars = qw ( RESULTDIR PDBDIR FASTADIR APBSDIR FPOCKET SRC MATCH3D ANNDIR UNIPROT PREMONITION HELIXDIR DSSP CONFIGGRP BLASTOUT BLASTDB);
   my @ret ;
   print STDERR "=============\n" if($verbose);
   foreach my $var (@vars){
       my $v  = $ENV{$var} or die "$var not set";
       #print  STDERR " $var = $v\n " ;
	   push @ret, $v ;
   }
   print STDERR  "\n===============\n" if($verbose);
   return @ret ;
}

sub util_readPDB{
    my ($pdb1) = @_ ;
	my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
    my $origpdb = $pdb1 ;
    $pdb1 = "$PDBDIR/$pdb1.pdb";

    my $pdb = new PDB();
    #$pdb->SetLogFile($ofh);
    $pdb->ReadPDB($pdb1);
    return $pdb ;
}



sub util_GetEnv{
	my ($l) = @_ ;
    my $ret = $ENV{$l} or die "Need to set environment variable $l"; 
	return $ret ;
}

sub util_make_list{
	my (@l) = @_ ;
	return \@l ;
}
sub util_make_table{
	my ($l) = @_ ;
	my $t = {};
	map {$t->{$_} = 1 ;} @{$l} ;
	return $t ;
}


sub util_wait_on_lockfiles{
    my ($ofh,$lockfiles) = @_ ;
    my @lockfiles = @{$lockfiles};
    print $ofh "echo Waiting for lock files to disappear\n";
    print $ofh "while(1) \n";
    print $ofh "sleep 5\n";
    my $str = "";
    map { $str.= " -e $_ || " ; } (@lockfiles);
    $str .= " -e $lockfiles[0]" ;
    print $ofh "if( $str)  then \n";
    print $ofh "else \n";
    print $ofh "break \n";
    print $ofh "endif \n";
    print $ofh "end \n";
}

sub util_assign_list2list{
     my ($list1,$list2,$assignment) = @_ ;
     my @list1 = @{$list1};
     my @list2 = @{$list2};
     map { $assignment->{$_} = [] if(!defined $assignment->{$_}) ; } @list1;
     my $num = @list1 ;
     my $cnt = 0 ;
     foreach my $l2 (@list2){
         my $l1 = $list1[$cnt];
         push @{$assignment->{$l1}} , $l2;
         $cnt++;
         $cnt = 0 if($cnt == $num);
     }
     #return $assignment;
}


sub util_two2thepowerof{
     my ($number) = @_;
     my $ret = 1; 
     foreach my $i (1..$number){
        $ret = 2*$ret;
     }
    return $ret ;
}


sub util_fullline{
    my ($fh,$firstline,$delim) = @_ ;
    if($firstline =~ /\s*$delim\s*$/){
        return ($firstline);
    }
    my $retline = $firstline ;
    while(<$fh>){
         chop ;
         $retline .= $_;
         if(/\s*$delim\s*$/){
             return ($retline);
         }
    }
    undef ;
}


sub util_print_vars{
     my ($dirname,$ofph) = @_ ;
     my @envs = qw(PERLLIB );
     map { print $ofph  "#$_ = $ENV{$_}\n"; } @envs ;
}


sub util_write_list_in_file{
     my ($filenm) = shift @_ ;
     my (@list) = @_ ;
     unlink $filenm ;
     my $fh = new FileHandle($filenm,O_CREAT|O_WRONLY) or croak " $! $filenm" ;
     map { print $fh "$_\n"; } @list ;
}

sub util_round {
    my($number) = shift;
    return int($number + .5);
}

sub util_get_pwd{
   my $PWD = getcwd ;
   return $PWD ;
}


sub util_execute_chop{
  my ($exec,$fname) = @_;
  my  $ret = ` $exec $fname `;
  chop $ret ;
  $ret ;
}
sub util_printAndDo{
   my ($what,$dry) = @_ ;
   my $comment = defined $dry ? " Will run ( this is dry run ) " : "Running" ;
   print STDERR "$comment $what ...\n";
   system($what) if(!defined $dry);
}
sub util_printAndWrite2Script{
   my ($what,$fh) = @_ ;
   ## 
   #$what =~ s/\"/\\\"/g ;
   croak " undefined file handle " if(!defined $fh);
   print STDERR "Running $what ...\n";
   print $fh "$what ; ";
}


sub util_exit_if_doesnt_exist{
  my ($fname) = @_;
  croak "File $fname does not exist. Quitting " if(!-e $fname);
}

sub util_get_tag{
  return "" if(!-e "CVS/Tag");
  my $fh = new FileHandle("CVS/Tag",O_RDONLY) or croak ;
  my $tag = <$fh>;
  $tag =~ s/^T//;
  $tag ;
}

sub util_print_n_log {
    my ($ofh,$msg) = @_ ;
	print $msg ;
	print $ofh $msg ;
}

sub util_get_user{
   my $user = `whoami` ;
   chop $user ;
   $user ;
}
sub util_get_cmdline{
    my ($exec , $list ) = @_ ;
	map { $exec = $exec . " $_ " ; } @$list ;
	$exec ;

}
sub util_num_lines{
  my ($file) = @_ ;
  return 0  if(!-e  $file);
  my $fh = new FileHandle($file,O_RDONLY) or croak ;
  my $num = 0 ; 
  while(<$fh>){ $num++ ; }
  $num ;
}

sub util_percentage{
   my ($a,$b,$justval) = @_ ;  
   croak "a is undefined" if(!defined $a);
   croak "b is undefined" if(!defined $b);
   return $a if(defined $justval);
   $a = 1 if($a eq 0 || !defined $a);
   $b = 1 if($b eq 0 || !defined $b);
   my $percent = ($a - $b)/$a ;  # Changing the diff so that we are in sync with harness eqn
   #my $percent = 1 - ($a/$b);
   #my $percent = ($a/$b)-1 ; 
   #int($percent*100) . "%" ;
   int($percent*100);
}

sub util_diff{
   my ($a,$b,$justval) = @_ ;  
   return $a if(defined $justval);
   $a = 1 if($a eq 0 || !defined $a);
   $b = 1 if($b eq 0 || !defined $b);
   my $percent = ($a - $b); 
}

sub util_log_base {
    my ($base, $value) = @_;
    if($value < $base){
        return 1;
    }
     my $val =  log($value)/log($base);
     my $finalval =ceil($val);
     return $finalval;
}

sub util_ceil{
    my ($value) = @_;
     my $finalval =ceil($value);
	 return $finalval;
}
sub util_floor{
    my ($value) = @_;
     my $finalval =floor($value);
	 return $finalval;
}

sub util_read_list_numbers{
    my @list = ();
    my ($file) = @_ ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " Error for file: $file $!" ;
    while(<$fp>){
         next if(/^\s*$/);
		 s/\s*//g;
         chomp ;
         croak " just expect a number " if(!/^[+-]?\d+(\.\d+)?$/);
         my ($num) = $_ ;
         push @list, $num ;
    }
    return @list ;
}

sub util_maketablefromfile{
    my ($file) = @_ ;
    my $fp = util_read($file);
	my $table ;
    while(<$fp>){
         next if(/^\s*$/);
         next if(/^\s*#/);
         my ($a,$b) = split ; 
         $table->{$a}= $b ;
    }
    return $table ;

}


sub util_read_list_words{
    my @list = ();
    my ($file) = @_ ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " $!" ;
    while(<$fp>){
         next if(/^\s*$/);
         next if(/^\s*#/);
         chop ;
         my @l = split ; 
         push @list, @l ;
    }
    return @list ;

}
sub util_read_list_sentences{
    my @list = ();
    my ($file) = @_ ;
	my $PWD = getcwd ;
    my $fp = new FileHandle($file,O_RDONLY) or croak " util_read_list_sentences :  $! $file $PWD" ;
    while(<$fp>){
         next if(/^\s*$/);
         chop ;
         push @list, $_ ;
    }
    return @list ;

}
sub util_get_abs_path{
    my ($fname) = @_ ; 
    my $pwd = util_get_pwd();
    return $fname if($fname =~ /^\s*\// || $fname =~ /^\s*\~/);
    my $abspath = $pwd . "/" . $fname ;
    return $abspath ;
}



sub util_get_time_from_string{
 my ($str) = @_ ;

 my ($a,$b,$c,$d) = ($str =~ /(\d+)\.(\d+)u\s*(\d+)\.(\d+)s/);
 return $a + $c ;

}

sub util_isDesignSuccessful{
    my ($logname) = @_ ;
    my $fp = new FileHandle($logname,O_RDONLY) or croak " $! $logname" ;
    my $foundSynthesized = "FAILED";
    while(<$fp>){
        if(/Finished synthesizing design/){
            $foundSynthesized = "PASSED";
        }
    }
    return $foundSynthesized;
}

sub util_parseTimelog {
   my ($logname) = @_ ;
   my $fp = new FileHandle($logname,O_RDONLY) or croak " $logname $!" ;
   my $lastline;
   while(<$fp>){
     $lastline = $_;
   }
   return util_get_time_from_string($lastline);
}

sub util_open_or_append{
     my ($outfile)= @_;
	 die "Blank file name " if($outfile =~ /^\s*$/);
	 my $fh ;
	 if(! -e $outfile){
         $fh = util_write($outfile);
	 }
	 else{
         $fh = util_append($outfile);
	 }
     return $fh ;
}

sub util_append{
     my ($outfile)= @_;
	 die "Blank file name " if($outfile =~ /^\s*$/);
     #unlink $outfile;
     my $fh = new FileHandle($outfile,O_WRONLY|O_APPEND) or croak " could not write file $outfile as $!" ;
     return $fh ;
}

sub util_write{
     my ($outfile)= @_;
	 croak "not defined" if(!defined $outfile);
     unlink $outfile;
     my $fh = new FileHandle($outfile,O_CREAT|O_WRONLY) or croak " could not write file $outfile as $!" ;
	 #print "Writing to $outfile\n";
     return $fh ;
}

sub util_writelist2file{
     my ($outfile,$list)= @_;
     unlink $outfile;
     my $fh = new FileHandle($outfile,O_CREAT|O_WRONLY) or croak " could not write file $outfile as $!" ;
	 foreach my $x (@{$list}){
	 	print $fh "$x\n";
	 }
     return $fh ;
}

sub util_read{
     my ($outfile)= @_;
     my $fh = new FileHandle($outfile,O_RDONLY) or croak " could not read file $outfile as $!";
     return $fh ;
}

sub util_wget{
     my ($file)= @_;
     util_printAndDo("wget --no-proxy $file");
}

sub util_makeCSH{
     my ($ofh)= @_;
     print $ofh "#!/bin/csh -f\n";
}

sub util_get_pdb{
     my ($dir,$file)= @_;
	 my $fname = $dir . "/" . $file . ".pdb" ;
	 if(-e $fname){
	 }
	 else{
	 	util_wget($fname); 
	 }
	 return $fname ; 
}
                                                                                                                                                             
sub util_pick_random_from_list{
    my ($list) = @_ ;
    my @temporaries = @{$list};
    my $num = @temporaries ;
    my $r = floor($num*rand());      
    my $operator = $temporaries[$r] or croak ;
    return $operator ;
}

sub util_pick_n_random_from_list{
    my ($list,$n) = @_ ;
	my @l ; 
	my $done ; 
	my $N = @{$list} ; 
	foreach my $i (1..$n){
		my $over = 0 ; 
		while (!$over){
		    my $x = util_pick_random_from_list($list);
			#print "$x $i $N\n";
			if(! exists $done->{$x}){
		        push @l, $x ; 
			    $over =  1 ; 
				$done->{$x} = 1 ; 
			}
		}
	}
	return \@l ; 
}




sub util_enter_maxcutoff{
	print STDERR " Enter max cutoff number \n";
	my $qnum = <> ;
	chop $qnum ; 

	my $ret = 1000000 ;
	$ret = $qnum if($qnum !~ /^\s*$/);
	return $ret ;
}
sub util_is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

sub util_is_float {
      defined $_[0] && $_[0] =~ /^[+-]?\d+(\.\d+)?$/;
}

sub util_EnterName{
	my ($default) = @_ ; 
	if(defined $default){
		print "Press enter to choose default name $default\n";
	}
	my $qname = <> ;
	chomp $qname ; 
	if(defined $default){
		if($qname =~ /^\s*$/){
			return $default ; 
		}
	}
	my @l = split " ",$qname; 
	if(@l != 1){
		print "Warning: Need a single name. Try again\n";
		return util_EnterName(); 
	}
	return $l[0]; 
}

sub util_EnterNumber{
	my $qnum = <> ;
	chomp $qnum ; 
	my @l = split " ",$qnum; 
	if(@l != 1 || !util_is_integer($l[0])){
		print "Warning: Need a single number. Try again\n";
		return util_EnterNumber(); 
	}

	return $l[0]; 
}
sub util_EnterTwoNumbers{
	my $qnum = <> ;
	chomp $qnum ; 
	my @l = split " ",$qnum; 
	if(@l != 2 || !util_is_integer($l[0]) || !util_is_integer($l[1])){
		print "Warning: Need 2 numbers. Try again\n";
		return util_EnterTwoNumbers(); 
	}

	return @l ;
}
sub util_EnterTwoNames{
	my $qname = <> ;
	chomp $qname ; 
	my @l = split " ",$qname; 
	if(@l != 2){
		print "Warning: Need 2 Names. Try again\n";
		return util_EnterTwoNames(); 
	}

	return @l ;
}


sub util_AddBeforeEach{
	my ($expr,@l) = @_ ; 
	my $str  = join " -$expr ",@l ;
	$str  = " -$expr $str";
	return $str ; 
}

sub util_ignoreIfDistanceIsLessthan{
	my ($num,@l) = @_ ; 
    my $iter = combinations(\@l, 2);
    while (my $c = $iter->next) {
        my @combo = @{$c} ; 
	    if(abs($combo[0] - $combo[1]) < $num){
			return 1 ; 
		}
	}
	return 0 ;
}
sub util_ignoreIfAnyAreEqual{
	my (@l) = @_ ; 
    my $iter = combinations(\@l, 2);
    while (my $c = $iter->next) {
        my @combo = @{$c} ; 
	    if($combo[0] == $combo[1]){
			return 1 ; 
		}
	}
	return 0 ;
}

sub util_split_list_numbers{
	my ($l,$point) = @_ ; 
	my @l = @{$l} ; 

	my @lte ;
	my @gt ;
	foreach my $e (@l){
		if($e <= $point){
			push @lte, $point ;  
		}
		else{
			push @gt, $point ;  
		}
	}
	my $n1 = @lte ;
	my $n2 = @gt ;
	return (\@lte, \@gt,$n1,$n2);
}


sub util_read_Mapping_PDB_2_SWISSPROT{
	my ($infile) = @_;
    my $info = {};
    my $ifh = util_read($infile);
	print "Reading mapping $infile\n";

    #ignore first two lines
    <$ifh>;
    <$ifh>;
	my $uniqueEC ; 
	my $uniqueSP ; 
    while(<$ifh>){
         next if(/^\s*$/);
	     s/|//g;
	     s/,//g;
	     s/,//g;
	     my ($nm,$j1,$chainid,$j2,$n1,$j3,$n2,$j4,$swissprot,$jjj,$ec) = split ;
		 #print " $swissprot $ec\n";
    
		     if($chainid eq "A"){
	            if(!defined $info->{$nm}){
		           $info->{$nm} = {} ;
	           }
	  	         #print "$nm $chainid $swissprot  \n";
	             $info->{$nm}->{SWISSPROT} = $swissprot ;
	             $info->{$nm}->{EC} = [] if(!defined $info->{$nm}->{EC});
	             $info->{SWISSPROT2EC}->{$swissprot} = [] if(!defined $info->{SWISSPROT2EC}->{$swissprot}) ;
	             if($ec ne "0.0.0.0"){
	                push @{$info->{$nm}->{EC}}, $ec ; 
					push @{$info->{SWISSPROT2EC}->{$swissprot}},$ec ;
				    $uniqueSP->{$swissprot} = 1 ; 
				    $uniqueEC->{$ec} = 1 ; 
				 }
			}
    }
    close($ifh);
	return ($info,$uniqueEC,$uniqueSP) ; 
}

sub util_getECfromPDB{
	my($info,$nm) = @_ ; 
    $nm = lc($nm);
	my $ec =  $info->{$nm}->{EC} ; 
    return undef if(!defined $ec);
	return $ec ; 
}


sub util_filter_basedon_EC{
	my ($info,$uniqueEC,$uniqueSP,$list,$ofh,$ignore0000) = @_ ; 
	my @list = @{$list};
    my $N = @list  ;
    my $ecdone = {};
    my $spdone = {};
    my $pdbTable = {};
    my $cnt = 0 ;
    my $ignored = util_write("ignored");
    foreach my $UCPDBID  (@list){

       my $pdbid = lc($UCPDBID);
	   if(!defined $info->{$pdbid}){
	   	  print $ignored "PDB $pdbid not found in mapping \n";
		  next ;
	   }
       my @ec = @{$info->{$pdbid}->{EC}};
	   my $sp = $info->{$pdbid}->{SWISSPROT};
	   my $N = @ec ; 
	   if($N > 1){
	   	  print $ignored "PDB $pdbid ignored as it has more than 1 :$N ECS\n";
		  next ;
	   }
	   if(@ec == 0){
	        my @ECCC = @{$info->{SWISSPROT2EC}->{$sp}}; 
			if(@ECCC == 0){
				print $ignored "Ignoring pdb $pdbid with swiss $sp as there are no ECS\n";
				next ; 
			}
	   		print "Did not find EC number for $pdbid\n";
			@ec = @ECCC ;
       }
   	   print $ofh "$UCPDBID\n" if(defined $ofh);
	   $spdone->{$sp} = $pdbid ; 

	   my $added = 0 ; 
	   print "did not find EC for $pdbid \n" if(@ec == 0);
	   foreach my $ec (@ec){
	   		  #if(exists $ecdone->{$ec}){
			  	  #print $ignored "$ec exists already for $pdbid $ecdone->{$ec} \n";
				  #next ;
			  #}
			  $added = 1 ;
	          $ecdone->{$ec} = [] if(!defined $ecdone->{$ec});
	          push @{$ecdone->{$ec}}, $pdbid ; 
		}

	    $pdbTable->{$pdbid} = $sp ; 
	    $cnt++ if($added);

   }
   print "Wrote ignored PDBS in ignored\n";

   return ($ecdone,$spdone,$pdbTable,$cnt); 
}


sub util_getPDBID_basedon_SP{
	my ($info,$uniqueEC,$list,$id) = @_ ; 
    my ($ecdone,$spdone,$pdbTable,$cnt)= util_filter_basedon_EC($info,$uniqueEC,$list);
	if(exists $pdbTable->{$id}){
		return $id ; 
	}
	else{
       my $sp = $info->{$id}->{SWISSPROT};
	   if(exists $spdone->{$id}){
		   return $spdone->{$id}  ;  
	   }
	   else{
		   return undef ; 
	   }
	}
}


###########################################################################
# Parses a file and returns blocks between start and end
###########################################################################

sub util_ParseBlocks{
	my ($infile,$start,$end) = @_;
    my @blocks ;
    my $ifh = util_read($infile);

    while(<$ifh>){
         if(/$start/){
		 	my @block ;
			push @block,$_;
			while(!/$end/){
                $_ = <$ifh>;
                if(!$_){
					print "NULL - hnece lasting\n";
					 last ;
				}
                if(/$start/){
                    @block = ();
                    push @block,$_;
                }
                else{
                    push @block,$_;
                }

			}
			push @blocks, \@block ;
		 }
    }
    close($ifh);
	return \@blocks ;
}

sub util_ParseBlockForString{
	my ($block,$str) = @_;
    my @lines ;
    foreach my $i (@{$block}){
		$_ = $i ;
         if(/$str/){
			push @lines,$_;
		 }
    }
	return \@lines ;
}



#############################################################
## Calculate the mean and SD of a matrix, columnwise 
## Also takes a range, and tells how many lie outside that range.
## So this range ideally makes sense for one column only
#############################################################
sub util_ProcessRowAndColumnsForMean{

   my ($rows,$top,$low)= @_ ;
   my @rows = @{$rows};
   my $nrows = @rows - 1;
   my @cols ; 


    #print "Number of rows -1 = $nrows\n";
	my $once = 1 ;
    foreach my $i (0..$nrows){
	    my $row = $rows[$i];
	    my @row = @{$row};
	    my $ncol = @row - 1 ;
        if($once){
		   print "Number of column -1 = $ncol\n";
		   $once = 0 ;
		}
        foreach my $j (0..$ncol){
		    if(!defined $cols[$j]){
			    my @l = ();
			    $cols[$j] = \@l;
		    }
		    push @{$cols[$j]}, $row[$j];
	    }
    }
    
	my @means ; 
	my @sds ; 
	my $NN = @cols ;
    foreach my $container (@cols){
        my $mean = Math::NumberCruncher::Mean($container) or warn "Mean not found" ;
        my $sd = Math::NumberCruncher::StandardDeviation($container) or warn "sd not found" ;
		next if(!defined ($mean && $sd));
		push @means, util_format_float($mean,1) ;
		push @sds, util_format_float($sd,1) ;
        my $incnt = 0 ;
        my $outcnt = 0 ;
	    foreach my $i (@{$container}){
		    $i = -($i);
		    if($i > $low && $i < $top){
			    $incnt++;
		    }
		else{
			    $outcnt++;
		    }
	    }
	    print "mean = $mean sd = $sd \n";
	    print "incnt = $incnt outcnt = $outcnt \n";
    }
	return (\@means,\@sds);

}


sub util_format_float{
	my ($d,$v) = @_; 
	croak if(!defined $d);
	if (defined $v && $v eq 3) {return sprintf("%8.3f", $d);}
	if (defined $v && $v eq 1) {return sprintf("%8.1f", $d);}

	return sprintf("%8.3f", $d); 
}

sub util_mysubstr {
	my ($len,$str,$start,$end) = @_ ; 
	if($len < $start ){
		die "dying $str" if(!($str =~ /TER/));
	    return "" ;
	}
	my $diff = $end == 1 ? $end :  $end - $start + 1 ; 
	$start-- ; 
	#print "($str,$start,$diff)\n";
	my $s =  substr($str,$start,$diff);
	$s =~ s/\s*//g ; 
	return $s ; 
}

sub util_mysubstrDontStripSpace {
	my ($len,$str,$start,$end) = @_ ; 
	if($len < $start ){
		die "dying $str" if(!($str =~ /TER/));
	    return "" ;
	}
	my $diff = $end == 1 ? $end :  $end - $start + 1 ; 
	$start-- ; 
	return substr($str,$start,$diff);
}


sub ParseAPBSResult{
   my ($size,$infile,$listfile) = @_ ;
   my $ifh = util_read($infile);
   my $list = {};
   if(defined $listfile){
       my @list= util_read_list_sentences($listfile);
       map { $list->{$_} = 1 ; } @list ;
   }
   

   my $blocks = util_ParseBlocks($infile,"Starting read","Ending read");
   my @blocks = @{$blocks};

   my $len = @blocks ; 
   print STDERR "\%There were $len blocks\n";




    my @diffs = ();
    foreach my $block (@blocks){

	chomp $block ;
	    my $lines = util_ParseBlockForString($block,"Resultfile is");
	    if(@{$lines} == 0){
		    die ;
	    }
	    my $line = shift @{$lines}; 
	    my ($pdb)  = ($line =~ /\/(....)\.pdb.out/);
	    ($pdb)  = ($line =~ /(....)\.pdb.out/) if(!defined $pdb);
	    if(defined $listfile && defined $pdb ){
		     if(!exists $list->{$pdb}){
			    print "nexting\n";
				#return undef ;
			 	next;
			 }
	    }
    
    
	    $lines = util_ParseBlockForString($block,"potential");
	    if(@{$lines} != $size){
			print "nexting as size $size doestn match\n";
				return undef ;
		    next ;
	    }
    
       my @vals  ; 
       my $somethingwrong = 0 ; 
       while($line = shift @{$lines}){
		    my ($val) = ($line =~ /potential\s*=\s*(.*)/);
		    if(!defined $val || $val =~ /^\s*$/){
			    $somethingwrong = 1 ;
			    last ; 
		    }
		    #print $ofh "$val ";
		    push @vals, $val;
	    }
    
	    next if($somethingwrong);
    
    
        my $iter = combinations(\@vals, 2);
        my @diff  ;
        while (my $c = $iter->next) {
                my @combo = @{$c} ; 
			    my ($a,$b) = @combo ; 
			    my $d= $a - $b ;
			    #print " $d= $a - $b ; \n";
			    $d=util_format_float($d,1);
			    push @diff , $d; 
	    }
	    push @diffs , \@diff ;
    }
    
    #my ($means,$sds) = util_ProcessRowAndColumnsForMean(\@diffs,250,150);
	#util_PrintMeanAndSD($ofh,$means,$sds);


    return \@diffs ;
}

sub util_PrintMeanAndSD{
	my ($ofh,$means,$sds) = @_ ;
    #print $ofh "\\rowcolor{orange} \n";
	print $ofh "Mean & " ;
	util_printTableLine($ofh,$means);
    #print $ofh "\\rowcolor{orange} \n";
	print $ofh "SD & ";
	util_printTableLine($ofh,$sds);
}


sub util_printResult{
    my ($ofh,$bestScore,$bestList,$cnt) = @_ ; 
	carp "Best score not defined" if(! defined $bestScore);
	print $ofh  "#RESULT $cnt  SCORE - $bestScore\n";
	print  $ofh "# ";
    foreach my $atom (@{$bestList}){
          my ($res,$num,$type) = split "/", $atom ;
	      print $ofh  "-  $atom ";
    }
	print  $ofh "\n";
}


sub util_printTablePre{
    my ($ofh,$caption) = @_ ; 
	$caption = "XXX" if(!defined $caption);
    print $ofh "\\begin{center} \n";
    print $ofh "\\begin{table*} \n";
   print $ofh "\\caption { $caption }  \n";
    #print $ofh "\\rowcolors{1}{tableShade}{white} \n";
    
    print $ofh "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c            } \n";
	#print $ofh "\\rowcolor{red!50}  \n";
}

sub util_printTablePost{
    my ($ofh,$caption) = @_ ; 
	$caption = "XXX" if(!defined $caption);

   #print $ofh "\\rowcolor{green!50}  \n";
   #print $ofh "\\rowcolor{orange!50}  \n";
   print $ofh "\\end{tabular}  \\label{label} \n";
   #print $ofh "\\caption {\\bf $caption }  \n";
   print $ofh "\\end{table*} \n";
   print $ofh "\\end{center} \n";
}

sub util_printTableLine{
	my ($ofh,$list) = @_ ; 
	my @l ; 
	map { 
	   if(util_is_float($_)){
	   	  $_ = util_format_float($_,1);
	   }
	   push @l, $_ ; 
    } @{$list} ;

	my $str = join " & ", @l;
	print $ofh " $str \\\\ \n";
	#print $ofh "\\hline \n";
}

sub util_readAPBSPotentialFromStart{
	my ($protein,$apbsdir) = @_ ; 
    my $pqrfile = "$apbsdir/$protein/$protein.pqr";
    my $pqr = new PDB();
    $pqr->ReadPDB($pqrfile,"hackforpqr");
    my @pots = ();
    my $potential = "$apbsdir/$protein/pot1.dx.atompot";
    util_readAPBSPotential(\@pots,$potential);
	return ($pqr,\@pots);
}

sub util_readAPBSPotential{
	my ($pots,$potential) = @_ ; 
    my $ifh = util_read($potential);
	print "util_readAPBSPotential : Reading file $potential \n";
    while(<$ifh>){
         next if(/^\s*$/);
	     chop ;
	     my @l = split ",",$_;
	     #print "$l[3] \n";
	     push @{$pots}, $l[3] ;
    }
	close($ifh);
	
}

sub util_CmdLine{
    my ($nm,$var) = @_ ;
	#$var = ${$nm} ;
    die "Error: In command line parsing. Needed command line ==> $nm" if(!defined $var);
}

sub util_parsePDBSEQRESNEW{
	my ($infile,$all,$writedata,$ignorerepeats) = @_ ; 
    my $ifh = util_read($infile);
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

	my $seenNm ; 
	my $info  = {}; 
	my $infoSeq2PDB = {} ; 
	my $mapChainedName2Name = {} ; 
	my $uniqueseqs = {} ; 
	
	my $nm ;
	my $CNT = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
	     if(/>/){
	        ($nm) = />(.*)\s*/;
			$nm = uc($nm);
			die if(!defined $nm); 
			warn "Repeat $nm " if(exists $info->{$nm});
			$CNT++;
			$info->{$nm} = "" ;
			next ;
		}


		## add to seq
		die if(/>/);
		if(!defined $_){
		    die "$nm $CNT " ;
		}
		chomp;
		die $_ if(!defined $nm);
		$info->{$nm} = $info->{$nm} . $_ ;

	 }

	 foreach my $k (keys %{$info}){
	 	my $seq = $info->{$k};
		my $length = length($seq);
		#print "$length $k \n";
		$infoSeq2PDB->{$seq} = [] if(!defined $infoSeq2PDB->{$seq}) ;
		push @{$infoSeq2PDB->{$seq}}, $k ;
	 }

	 my $N = (keys %{$infoSeq2PDB});
	 print "There were $CNT seqeunces to start with, and there are $N unique \n";

	 if(defined $writedata && $writedata){
         my $ofh = util_write("list.unique");

		 my $SORTING = {};
         foreach my $k (keys %{$infoSeq2PDB}){
		 	my $len = length($k);
			$SORTING->{$k} = $len ;
		 }


		 my @sl = sort { $SORTING->{$b} <=> $SORTING->{$a}} (keys %{$SORTING});
		 my $CNTwrite = 0 ;
         foreach my $k (@sl){
	        my @pdbswithseq  = sort @{$infoSeq2PDB->{$k}} ;
		 	my $len = length($k);
			#if($CNTwrite < 3 ){
			if(1){
			   #print "$len \n";
		       foreach my $pdb (@pdbswithseq){
				   ($pdb) = split " ", $pdb ;
			       #print "KKKKKKKK $pdb\n";
			       my $fastanm = "$FASTADIR/$pdb.ALL.1.fasta";
				   $CNTwrite++;
			       my $FH = util_write($fastanm);
                  print $FH "\>$pdb\n";
                  print $FH "$k\n";
		          close($FH);
       
		       }
		       my $v = $pdbswithseq[0];
	 	       print $ofh "$v\n";
			}
	     }
		 print "Wrote $CNTwrite finally\n";
             close($ofh);


         $ofh = util_write("list.allpdbschained");
         foreach my $k (keys %{$mapChainedName2Name}){
	 	    print $ofh "$k\n";
	     }
         close($ofh);

         $ofh = util_write("list.allpdbs");
         foreach my $k (keys %{$info}){
	 	    $k = uc($k);
	 	    print $ofh "$k\n";
	     }
         close($ofh);
	 }


	 return ($info,$infoSeq2PDB,$mapChainedName2Name) ;
}


sub util_parsePDBSEQRESOLLLLLL{
	my ($infile,$all,$writedata) = @_ ; 
    my $ifh = util_read($infile);
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

	my $seenNm ; 
	my $info  = {}; 
	my $infoSeq2PDB = {} ; 
	my $mapChainedName2Name = {} ; 
	my $uniqueseqs = {} ; 
	
    while(<$ifh>){
         next if(/^\s*$/);
	     if(/^\s*>/){
	        my ($nm,$chain,$type,$len,$fullnm) = parseSingleLinePDBSEQRES($_);
			if(!defined $type){
				($nm) = (/>(.*)\s+/);
				$type = "protein";
				$len = "50";
				$fullnm = $nm ;
				$chain = "A";
			}
			next if($type ne "protein");

			next if(!defined $nm); 


	        my $seq = <$ifh> ;
		    chomp $seq ; 
			next if(!length($seq));

			my $chainednm = uc($nm . $chain);
			$mapChainedName2Name->{$chainednm} = $seq ;

			if(! exists $infoSeq2PDB->{$seq}){
				$infoSeq2PDB->{$seq} = [];
			}
			push @{$infoSeq2PDB->{$seq}}, $chainednm ; 

			$uniqueseqs->{$seq} = $chainednm ;


			next if(exists $seenNm->{$nm}); 
			$seenNm->{$nm} = 1 ;


		    $info->{$nm} = {};
		    $info->{$nm}->{TYPE} = $type ;
		    $info->{$nm}->{LEN} = $len ;
		    $info->{$nm}->{FULLNM} = $fullnm ;
		    $info->{$nm}->{SEQ} = $seq ;

		    $info->{$nm}->{NM} = $nm ;


	 	}
	 }

	 if(defined $writedata && $writedata){
         my $ofh = util_write("list.unique");
         foreach my $k (keys %{$infoSeq2PDB}){
	        my @pdbswithseq  = sort @{$infoSeq2PDB->{$k}} ;
		    foreach my $pdb (@pdbswithseq){
				($pdb) = split " ", $pdb ;
				#print "KKKKKKKK $pdb\n";
			    my $fastanm = "$FASTADIR/$pdb.ALL.1.fasta";
			    my $FH = util_write($fastanm);
               print $FH "\>$pdb\n";
               print $FH "$k\n";
		       close($FH);
    
		    }
		    my $v = $pdbswithseq[0];
	 	    print $ofh "$v\n";
	     }
             close($ofh);


         $ofh = util_write("list.allpdbschained");
         foreach my $k (keys %{$mapChainedName2Name}){
	 	    print $ofh "$k\n";
	     }
         close($ofh);

         $ofh = util_write("list.allpdbs");
         foreach my $k (keys %{$info}){
	 	    $k = uc($k);
	 	    print $ofh "$k\n";
	     }
         close($ofh);
	 }


	 return ($info,$infoSeq2PDB,$mapChainedName2Name) ;
}
  



sub parseSingleLinePDBSEQRES{
	my ($line) = @_ ; 
	my ($nm,$chain,$type,$len,$fullnm) = ($line =~ /^.(....).(.)\s*mol:(\w+)\s*length:(\d+)\s*(.*)/);
	return ($nm,$chain,$type,$len,$fullnm) ;

}
sub util_Annotate{
	my ($file) = @_ ; 
	my ($outfile) = "$file". ".annotate";
	util_printAndDo("annotate.pl -in ~/pdb_seqres.txt -lis $file -out $outfile -cutoff 100 -anndis 10 ");

}

sub util_IsZero{
	my ($num) = @_ ; 
	if(abs($num) < $EPSILON){
		return 1 ; 
	}
	else{
		return 0 ; 
	}

}

sub util_Banner{
	my ($str) = @_ ; 
	print STDERR "===============================================\n";
	print STDERR "$str\n";
	print STDERR "===============================================\n";

}
sub util_PrintInfo{
	my ($str) = @_ ; 
	print STDERR "Info: $str\n";

}


sub util_getTmpFile{
    my $time = int(time() *rand());
    my $tmpfile =  "sandeeptmp.$time";
	return $tmpfile ;
}

sub util_percentages{
	my (@values) = @_ ;
	my $sum = 0 ; 
	#print " KKKKKKKKKKKKK " , @values , "\n";
	foreach my $v (@values){
	   $sum = $sum + $v ; 
	}
	my @l ; 
	foreach my $v (@values){
	   push @l, util_format_float(($v*100)/$sum,3) ; 
	}
	return @l ;
}

sub util_printHtmlHeader{
	my ($ofh,$header1,$header2) = @_ ;
    print $ofh "<html> \n";
    print $ofh "<h1>$header1</h1> \n";
    print $ofh "<body> \n";

	if(defined $header2){
    print $ofh "<html> \n";
    print $ofh "<h2>$header2</h2> \n";
    print $ofh "<body> \n";
	}
}



sub util_printHtmlEnd{
	my ($ofh) = @_ ;
    print $ofh "</body> \n";
    print $ofh "</html> \n";
}

sub util_HtmlizeLine{
	my ($line) = @_ ;
	chomp $line ; 
	#return $line  . "<br />" ;
	return $line ;

}

sub util_HtmlTableHead{
	my ($ofh,@headers) = @_ ;
    print $ofh "<table border=\"1\" cellpadding=\"5\" cellspacing=\"5\" width=\"100%\">\n";
    print $ofh "<tr>\n";
	foreach my $h (@headers){
        print $ofh "<th>$h</th>";
	}
    print $ofh "\n";
    print $ofh "</tr>\n";

}


sub util_HtmlTableEnd{
	my ($ofh) = @_ ;
	print $ofh "</table>\n";
}



sub util_HtmlTableCell{
	my ($str) = @_ ;
	return "<td>$str</td>";
}



sub util_table_print{
	my ($table) = @_ ;
	my $N = keys %{$table};
	print "There are $N entries \n";
	print "===========\n";
	foreach my $k (sort { $a <=> $b} keys %{$table}){
		print "Tableprint $k $table->{$k} \n";
	}
}


sub util_MakeLink{
	my ($nm,$link) = @_ ;
    return  "<a href=\"$link\"> $nm</a>";
}

sub util_EC_CorrelatePDBS{
	my ($info,$a,$b,$VALUE) = @_ ; 
	print "util_EC_CorrelatePDBS for $a $b\n";
	my $x = util_getECfromPDB($info,$a);
	my $y = util_getECfromPDB($info,$b);
	if(!defined $x || !defined $y){
		print "Could not get EC for $a \n" if(!defined $x);
		print "Could not get EC for $b \n" if(!defined $y);;
		return undef;
	}

	my @ec1 = @{$x};
	my @ec2 = @{$y};
	my $ec1 = $ec1[0];
	my $ec2 = $ec2[0];
	if(!defined $ec1 || !defined $ec2){
		print "Could not get EC for $a \n" if(!defined $ec1);
		print "Could not get EC for $b \n" if(!defined $ec2);;
		return undef;
	}
	print "util_EC_CorrelatePDBS for $ec1 $ec2\n";

	$ec1 =~ s/\./YYY/g;
	$ec2 =~ s/\./YYY/g;
	my @l1 = split "YYY", $ec1 ;
	my @l2 = split "YYY", $ec2 ;
   
   my $N1 = @l1 -1 ;
   my $N2 = @l2 -1 ;
   my $score = abs($N1 - $N2);

   foreach my $n (0..$N1){
			my $v1 = $l1[$n] ; 
			my $v2 = $l2[$n] ; 
			if($v1 != $v2){
				$score = $score + $VALUE/($n+1) ; 
			    print " Found mismatch $n $score $v1 $v2\n";
				last ;
			}
   }
   return $score ; 
}

sub util_EC_AddPDB{
	my ($info,$pdb) = @_ ; 
	my $ec = util_getECfromPDB($info,$pdb);
	return $ec if(!defined $ec);
	my @l = split "\.", $ec ;
    
   my $N = @l -1 ;
   my $obj = {} ;
   foreach my $n (0..$N){
			my $v = $l[$n] ; 
		    $obj->{$n}->{$v} = [] if(!exists $obj->{$n}->{$v}) ;
			push @{$obj->{$n}->{$v}}, $pdb ; 
   }
}


sub util_EC_CreateLevels{
	my ($ecdone) = @_ ; 
    my $MAXLEVEL = 0 ;
    my $obj = {} ; 
    foreach my $ec (sort keys %{$ecdone}){
	    my @l = split "\.", $ec ;
	    $MAXLEVEL = @l if(@l > $MAXLEVEL);
    
	    my $N = @l -1 ;
	    foreach my $n (0..$N){
		    my $obj->{$n} = [];
	    }
    }
	return ($obj,$MAXLEVEL) ;
}


sub util_readPeptideInfo{
   my ($info,$nm,$infile) = @_ ;
   $info->{$nm} = {};
   die "ERROR: util_readPeptideInfo $infile does not exist" if(! -e $infile);
   my $ifh = util_read($infile);
   while(<$ifh>){
        next if(/^\s*$/);

        if(/Average Residue Weight/){
           my ($charge) = (/Charge\s*=\s*(.*)/) or die;
           $info->{$nm}->{CHARGE} = $charge ;
       }
	   if(/Residues = /){
	   	    my ($nres) = (/Residues = (\d+)/) or die;
			$info->{$nm}->{NRES} = $nres ;
	   }
        if(/^(Molecular|Basic|Acidic|Polar)/i){
		   my (@l) = split ;
		   my $N = @l -1 ;
           $info->{$nm}->{$l[0]} = $l[$N]; ;
		}
   }
   if(defined $info->{$nm}->{Acidic} && defined $info->{$nm}->{Basic}){
       $info->{$nm}->{AcidBasic} = $info->{$nm}->{Acidic} + $info->{$nm}->{Basic} ;
       $info->{$nm}->{PAB} = $info->{$nm}->{Acidic} + $info->{$nm}->{Basic} + $info->{$nm}->{Polar} ;
   }
   close($ifh);
}


sub util_ReadLine{
		   my ($LINE) = @_ ;
		   my $len = length($LINE);
           my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z );

           $atomstr = util_mysubstr($len,$LINE,1 ,  6);
           $serialnum = util_mysubstr($len,$LINE,7 , 11);


           $atomnm = util_mysubstrDontStripSpace($len,$LINE,13 , 16);
           $alt_loc = util_mysubstr($len,$LINE,17,1);
           $resname = util_mysubstr($len,$LINE,18 , 20);
           $chainId = util_mysubstr($len,$LINE,22,1);
           $resnum = util_mysubstr($len,$LINE,23 , 26);
           $codeforinsertion = util_mysubstr($len,$LINE,27,1);
           $x = util_mysubstr($len,$LINE,31 , 38);
           $y = util_mysubstr($len,$LINE,39 , 46);
           $z = util_mysubstr($len,$LINE,47 , 54);
		   return ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ); 
}
sub util_getECID{
	my ($ec,$level) = @_ ; 
	$ec =~ s/\./YYY/g ; 
	my @l = split "YYY",$ec ; 
	my $str = "";
	my $first = 1 ; 
	foreach my $i (1..$level){
		if(!$first){
			$str = $str . ".";
		}
		$str = $str . $l[$i -1];
		$first = 0 ; 
	}
	return $str ; 
}

sub util_uniq2 {
	my ($list) = @_ ;
    my %seen = ();
    my @r = ();
    foreach my $a (@{$list}) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
    }
    return @r;
}
sub vSetVerbose{
	my ($val) = @_ ;
	my $ret = $ENV{VERBOSE};
	if(defined $ret){
		$ENV{VERBOSE} = $val ;
	}
}
sub vIncrVerbose{
	my $ret = $ENV{VERBOSE};
	if(defined $ret){
		$ENV{VERBOSE} = $ENV{VERBOSE} + 1 ;
	}
}
sub vDecrVerbose{
	my $ret = $ENV{VERBOSE};
	if(defined $ret){
		$ENV{VERBOSE} = $ENV{VERBOSE} - 1 ;
	}
}

sub vprint{
	my ($str) = @_ ;
	my $ret = $ENV{VERBOSE};
	if(defined $ret && $ret > 0 ){
		my $tab = "";
		foreach my $i (1..$ret){
			$tab = $tab. "\t";
		}
	    print STDERR "$tab$str\n";
	}
}
sub vprintheader{
	my ($str) = @_ ;
	my $ret = $ENV{VERBOSE};
	if(defined $ret && $ret > 0){
		my $tab = "";
		foreach my $i (1..$ret){
			$tab = $tab. "\t";
		}
	    print STDERR "\n\n======================================================================\n";
	    print STDERR "$tab$str\n";
	    print STDERR "======================================================================\n";
	}
}
sub util_copy_table{
	my ($table) = @_ ; 
	my $ret = {};
	foreach my $k (keys %{$table}){
		$ret->{$k} = $table->{$k}; 
	}
	return $ret; 
}
sub util_GetFastaFiles{
	my ($FASTADIR,@pdbs) = @_ ; 
    my $mapinfo ;
    my $fhnotfound = util_write("list.notfound");
    my @files = ();
    foreach my $i (@pdbs){
	    $i = uc($i);
        $i =~ /^\s*$/;
        my @f = <$FASTADIR/$i*fasta>;
	    my $j = $f[0];
	    if(!defined $j){
		    warn "Did not find fasta file for pdb $i in dir $FASTADIR \n";
		    print $fhnotfound "$i\n";
            #push @ignored, $i;
		    next ;
	    }
	    print " pushing $j \n" ;
        push @files, $j;
	    $mapinfo->{$j} = $i ;
    }
	return ($mapinfo,@files) ; 
}

sub util_WriteFastaFromAtoms{
    my ($pdb,$allatoms,$fastafh,$origpdb) = @_ ; 
    my @allresidues ;
    my $allres  = {};
    foreach my $atom (@{$allatoms}){
		my $nm = $atom->GetName();
		next if($nm =~ /HOH/);
	    my $resnum = $atom->GetResNum();
	    my ($res) = $pdb->GetResidueIdx($resnum);
	    push @allresidues ,$res ; 
	    $allres->{$res->GetResNum()} = $res ; 
    }
    my $allresiduesN = @allresidues ;
    my $s2 = "";
    my $s1 = "";
    foreach my $i (sort  { $a <=> $b }  keys %{$allres}){
	    my $r = $allres->{$i} ; 
	    my $s = $r->PrintSingleLetter($pdb);
	    $s1 = $s1 .  "$s$i," ;
	    $s2 = $s2 . "$s" ;
    }
    print $fastafh "\>$origpdb.$s1;\n";
    print $fastafh "$s2\n";
}


sub util_WriteFastaFromResidueNumbers{
    my ($pdb,$resnumbers,$fastafh,$origpdb) = @_ ; 
    my @allresidues ;
    my $allres  = {};
    foreach my $resnum (@{$resnumbers}){
	    my ($res) = $pdb->GetResidueIdx($resnum);
	    push @allresidues ,$res ; 
	    $allres->{$res->GetResNum()} = $res ; 
    }
    my $allresiduesN = @allresidues ;
    my $s2 = "";
    my $s1 = "";
    foreach my $i (sort  { $a <=> $b }  keys %{$allres}){
	    my $r = $allres->{$i} ; 
	    my $s = $r->PrintSingleLetter($pdb);
	    $s1 = $s1 .  "$s$i," ;
	    $s2 = $s2 . "$s" ;
    }
    print $fastafh "\>$origpdb.$s1;\n";
    print $fastafh "$s2\n";
}

sub util_GetPotForAtom{
	my ($a,$pqr,$pots) = @_ ;
	croak "not defined" if(!defined $a);
	my $number = $a->GetResNum();
	#print STDERR "kkkkkkkkk $number \n";
	my $atomnm = $a->GetType();
    my ($aPqr) = $pqr->GetAtomFromResidueAndType($number,$atomnm) or croak "No PQR?" ;
    croak "could not find $number $atomnm" if(!defined $aPqr);

    my ($x,$y,$z) = $a->Coords();
    my ($x1,$y1,$z1) = $aPqr->Coords();
	if(1 && !util_IsZero($x-$x1+$y -$y1+$z-$z1)){
			# this will mismatch as rotation is done
             #warn "Warning: $x,$y,$z $x1,$y1,$z1 do not match" ; 
	}

    #my ($i1) = $a->GetIdx();
    my ($i2) = $aPqr->GetIdx();
    #imp -1 
	my $NPots = @{$pots};
	#$a->Print();
    my $pot = $pots->[$i2-1] or die "Expected to find potential for residue number $number i2=$i2 and $NPots";
	return $pot ; 
}

sub util_ReadPdbs{
	my ($PDBDIR,$APBSDIR,$readpotential,@P) = @_ ; 
	my @ret ; 
	die "Expected at least one protein" if(!@P);
	foreach my $p1 (@P){
         my $file1 = "$PDBDIR/$p1.pdb";
         my $PPP = new PDB();
         $PPP->ReadPDB($file1);
	     print STDERR "Reading $p1\n";

	     my $pqrfile = "$APBSDIR/$p1/$p1.pqr";
	     my $pqr = new PDB();
	     my @pots = ();
		 if($readpotential){
		      #$pqr->SetReadCharge(1);
	          $pqr->ReadPDB($pqrfile,"hack");
	          #$pqr->ReadPDB($pqrfile);
	          my $potential = "$APBSDIR/$p1/pot1.dx.atompot";
	          util_readAPBSPotential(\@pots,$potential);
		 }

         

		 my $info->{PDBNAME} = $p1 ;
		 $info->{PDBOBJ} = $PPP ;
		 $info->{PQR} = $pqr ;
		 $info->{POTS} = \@pots ;
		 push @ret , $info ; 
	}
	return @ret ;
}

sub util_WriteClustalAln{
	my ($protein1,$matches,$matchedProteins,$alnfh,$map,$annMap) = @_ ; 
	my @matchedProteins = @{$matchedProteins} ; 
	my $N = @matchedProteins ; 
	my $halfN = $N/2 ;
	my $quarterN = $halfN/2 ;
	my @matches = @{$matches} ; 
	my $first = 1 ; 
	my $EXTEND = {};
	while(@matches){
		my $protein2 = shift @matchedProteins ; 
		my $MATCH = shift @matches ; 
		if($first){
              print $alnfh "CLUSTAL 2.1 multiple sequence alignment\n\n\n";
			  my $pdb1 = $map->{$protein1} ;
	          my $mappedname = exists $annMap->{$protein1} ? $protein1. "." . $annMap->{$protein1} : $protein1 ;
		      my $XXX =  sprintf ( "%-12s",  $mappedname); 
              print $alnfh "$XXX";
			  my $CNT = 0 ; 
			  foreach my $resnum (sort {$a <=> $b}  keys %{$MATCH}){
			  	 $CNT++ ; 
				 $EXTEND->{$CNT} = {};
				 $EXTEND->{$CNT}->{CNT} = 1 ; 
				 $EXTEND->{$CNT}->{MATCHES} = {};
	             my ($res) = $pdb1->GetResidueIdx($resnum);
	             my $s = $res->PrintSingleLetter($pdb1);
				 $EXTEND->{$CNT}->{MATCHES}->{$s} =  1 ; 
                 print $alnfh "$s";
			  }
              print $alnfh "\n";
			  $first = 0 ;
		}


		my $pdb2 = $map->{$protein2} ;
	    my $mappedname = exists $annMap->{$protein2} ? $protein2. ".". $annMap->{$protein2} : $protein2 ;
		my $XXX =  sprintf ( "%-12s",  $mappedname); 
        print $alnfh "$XXX";
		my $CNT = 0 ; 
		foreach my $k (sort {$a <=> $b}  keys %{$MATCH}){
			 $CNT++ ; 
			 my $resnum = $MATCH->{$k} ;
			 my $s = "-";
			 if($resnum ne "-"){
				  $EXTEND->{$CNT}->{CNT} =  $EXTEND->{$CNT}->{CNT} + 1 ; 
	              my ($res) = $pdb2->GetResidueIdx($resnum);
	              $s = $res->PrintSingleLetter($pdb2);
				  $EXTEND->{$CNT}->{MATCHES}->{$s} =  1 ; 
			 }
             print $alnfh "$s";
		}
        print $alnfh "\n";
	}

	foreach my $CNT (sort {$a <=> $b}  keys %{$EXTEND}){
		my $cnt = $EXTEND->{$CNT}->{CNT} ;
		if($cnt >= (3*$quarterN)){
		   my $res = $EXTEND->{$CNT}->{MATCHES} ;
		   my @keys = (keys %{$res});
		   my $str = join "|" , @keys ;
		   $str = "(" . $str . ")";
		   print "$CNT $cnt $str ========= \n";
		}
	}

	
}

sub util_GetClosestAtoms_intwoPDBs{
    my ($p1,$p2,$PDBDIR,$maxdist) = @_ ; 
    my $file1 = "$PDBDIR/$p1.pdb";
    my $file2 = "$PDBDIR/$p2.pdb";
    my $pdb1 = new PDB();
    $pdb1->ReadPDB($file1);
    my $pdb2 = new PDB();
    $pdb2->ReadPDB($file2);
    
    my @reslist = $pdb1->GetResidues();
    
    my @atoms ;
    while(@reslist){
        my $r = shift @reslist ;
        my @aaa = $r->GetAtoms();
    	push @atoms, @aaa;
    }
    
    my @resultsall ; 
    my $DONE ;
    foreach my $atom (@atoms){
         my $list = util_make_list($atom);
         my ($junk,$neighatoms)  = $pdb2->GetNeighbourHoodAtom($list,$maxdist);
    	 foreach my $r (@{$neighatoms}){
    	       my $atomstr = $r->GetAtomStr();
    		   next if($atomstr eq "HETATM");
    
    	        my $d = $atom->Distance($r) ;
    			my $nm = $atom->GetName() . " -> " . $r->GetName();
    			my $info = {};
    	        $info->{NAME} = $nm;
    	        $info->{SCORE} = $d;
    	        push @resultsall, $info ;
    	 }
    
    }
    
    my @resultssorted = sort { $a->{SCORE} <=> $b->{SCORE} } @resultsall ;
    my $CNT = 0 ; 
    foreach my $r (@resultssorted){
    	$CNT++;
    	my $nm = $r->{NAME};
    	my $score = $r->{SCORE};
    	print STDERR "NAME = $nm score = $score \n";
    	last if($CNT eq 10);
    }
    return @resultssorted ;
}


sub util_Ann2Simple{
	my ($infile,$outfile) = @_; 
    my $ofh ;
    $ofh = util_write($outfile) if(defined $outfile);
    my $ifh = util_read($infile);
	my @retlist  ; 
    while(<$ifh>){
         next if(/^\s*$/);
         chop ;
	     if(/^POINTS/){
	 	    s/POINTS//g;
		    s/\// /g;
		    my @l = split " ", $_ ; 
		    while(@l){
			    my $a = shift @l ;
			    my $b = shift @l ;
			    #$a =~ s/\s*//g;
			    #$b =~ s/\s*//g;
			    print  $ofh "$b$a " if(defined $ofh);
				push @retlist, $a ;
		    }
		    print $ofh "\n" if(defined $ofh);
	     }
    }
	return @retlist ;
}

sub util_GetPotDiffForAtoms{
     my ($pdb1,$pqr1,$pots1,$a,$b) = @_ ; 
     my $pota = util_GetPotForAtom($a,$pqr1,$pots1) ;
     my $potb = util_GetPotForAtom($b,$pqr1,$pots1) ;
	 my $diff = $pota - $potb ;
	 return $diff;
}

sub util_GetPotDiffForResidues{
    my ($pdb1,$pqr1,$pots1,$res1,$res2,$what) = @_ ; 
	my $a = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),$what);
	my $b = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),$what);
	my $diff = util_GetPotDiffForAtoms($pdb1,$pqr1,$pots1,$a,$b);
}

sub util_ProcessSingleLine{
     my ($pdb1,$pqr1,$pots1,$line,$isnew) = @_ ; 
	 my $atomlist = $pdb1->ParseResultLine($line,$isnew);
	 return util_ProcessSingleLineAtomlist($pdb1,$pqr1,$pots1,$atomlist);
}

sub util_ProcessSingleLineAtomlist{
     my ($pdb1,$pqr1,$pots1,$atomlist) = @_ ; 
	 my @names ; 
	 foreach my $a (@{$atomlist}){
        my $pot = util_GetPotForAtom($a,$pqr1,$pots1) ;
		push @names, $a->GetName();
	 }
	 my $name = join ",",@names ; ;
	 my @dist = @{$pdb1->DistanceInGivenSetOfAtoms($atomlist)};
	 my @pots = @{$pdb1->PotInGivenSetOfAtoms($atomlist,$pqr1,$pots1)};
	 return (\@dist,\@pots,$name);
}

sub util_FindRmsd{
	my ($pdb1,$pdb2) = @_ ;
    
    my @res = $pdb1->GetResidues();
    my $N = @res;
    my $cnt = 0 ;
    my $sum = 0 ;
    my $cntmatch = 0 ; 
    foreach my $res (@res){
        next if($res->GetAtomStr() ne "ATOM");
	    my $resnumA = $res->GetResNum() ;
	    my $resnumB = $res->GetResNum() ;
        my $CAatom1 = $pdb1->GetAtomFromResidueAndType($resnumA,"CA");
        my $CAatom2 = $pdb2->GetAtomFromResidueAndType($resnumB,"CA");
	    my $d = util_format_float($pdb1->DistanceAtoms($CAatom2,$CAatom1),1);
        $cnt++;
	    $sum = $sum + $d * $d ; 
    }
    
    my $rmsd = util_format_float(sqrt($sum/$cnt),3) ; 
    print  " $rmsd $cnt\n";
    return $rmsd ;
}

sub util_GetMeanSD{
    my ($container) =@_ ;
	my $mean = util_format_float(Math::NumberCruncher::Mean($container),0) or warn "Mean not found" ;
	my $sd = util_format_float(Math::NumberCruncher::StandardDeviation($container),0) or warn "sd not found" ;
	return ($mean,$sd);
}

sub util_FindRmsdAllAtoms{
	my ($pdb1,$pdb2) = @_ ;
    
    my @atoms = $pdb1->GetAtoms();
    my $N = @atoms;
    my $cnt = 0 ;
    my $sum = 0 ;
    my $cntmatch = 0 ; 
    foreach my $atom1 (@atoms){
	    my $type = $atom1->GetType();
	    my $resnum = $atom1->GetResNum() ;
        my $atom2 = $pdb2->GetAtomFromResidueAndType($resnum,$type);
		next if (!defined $atom2);
	    my $d = util_format_float($pdb1->distanceAtoms($atom1,$atom2),1);
        $cnt++;
	    $sum = $sum + $d * $d ; 
    }
    
    my $rmsd = util_format_float(sqrt($sum/$cnt),3) ; 
    print  " $rmsd $cnt\n";
    return $rmsd ;
}
sub util_sortsingleString{
   my $s = shift ;
   my @sl = split "", $s ;
   my @XX = sort @sl ; 
   
   my $rev = join "", @XX ;
   return $rev ; 
}

sub util_AreResiduesContinuous{
	my ($pdb1,$res1,$res2,$eitherispro) = @_ ;
	my $CA1 = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),"CA");
	my $CA2 = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),"CA");
    my $d1 = util_format_float($pdb1->DistanceAtoms($CA1,$CA2),1);
    my $diff1 = abs ($d1 - 3.8);
	if($eitherispro){
        my $XXX = abs ($d1 - 3);
		$diff1 = $XXX if($XXX < $diff1);
	}
    return (0,$d1,$diff1) if($diff1 > 0.3);
    return (1,$d1,$diff1) ;
}

sub util_IsPro{
	my ($res) = @_ ;
	my $name = $res->GetName();
	return 1 if($name eq "PRO");
	return 0 ; 
}

sub util_IsResisdueThisType{
	my ($res,$type) = @_ ;
	my $name = $res->GetName();
	return 1 if($name eq "$type");
	return 0 ; 
}
sub util_IsMSE{
	my ($res) = @_ ;
	my $name = $res->GetName();
	return 1 if($name eq "MSE");
	return 0 ; 
}
sub util_readfasta{
    my ($infile) =@_ ; 
    my $str = "";
    my $ifh = util_read($infile);
	my $firstline ; 
    while(<$ifh>){
         if(/^\s*>/){
	 	    $firstline = $_ ; 
		    next ; 
	     }
         next if(/^\s*$/);
	     chop ;
	     $str = $str . $_ ; 
    }
    $str =~ s/\s*//g;
    return ($str,$firstline) ; 
}

sub util_ExtractSliceFromFasta{
	my ($ofh,$infile,$start,$end) = @_ ; 
    my ($str,$firstline) = util_readfasta($infile);

    chop $firstline ;
    my $time = int(time() *rand());
    print "$time";
	$time = "$start.$end";
    
    print $ofh ">$time.$firstline  \n";
    
    my @l = split "", $str ; 
    my $N = @l ;
    print "$N sll \n";
    foreach my $i ($start..$end){
	    $i = $i - 1 ; 
	    print $ofh "$l[$i]";
    }
    print $ofh "\n";
}


## THIS NEEDS TO BE FIXED FOR FRAGMENTCOMPARE.PL ###
sub util_ReadAnnotateFile{
	my ($infile,$DIFF) = @_ ; 
	$DIFF = 50 if(!defined $DIFF);
    my $ifh = util_read($infile);
    my $info = {};
    while(<$ifh>){
         chop; 
         next if(/^\s*$/);
	     if(/^\s*(Repeat|Region)/){
	         my (@l) = split ; 
		     my $start = $l[1];
		     my $end = $l[3];
		     my $diff = abs($start - $end);
		     next if($diff > $DIFF);
    
		     my $anno = $_ ;
		     print "$start $end $anno \n";
		     $info->{$start} = $anno ;
		     $info->{$end} = $anno ;
	     }
	     if(/^\s*Modified/){
	         my (@l) = split ; 
		     my $start = $l[2];
		     my $anno = $_ ;
		     $info->{$start} = $anno ;
	     }
	 	    
    }
	return $info ;
}

sub util_ReadAnnotateFileFRAGAL{
	my ($infile,$DIFF) = @_ ; 
	$DIFF = 50 if(!defined $DIFF);
    my $ifh = util_read($infile);
    my $inforepeatregion = {};
    my $infomodified = {};
	my $CNT = 0; 
    while(<$ifh>){
         chop; 
         next if(/^\s*$/);
		 $CNT++;
	     if(/^\s*(Repeat|Region)/){
	         my (@l) = split ; 
		     my $start = $l[1];
		     my $end = $l[3];
		     my $diff = abs($start - $end);
		     if($diff > $DIFF){
			 	print "DIFF $diff is freater than $DIFF\n";
				next ;
			 }
    
		     my $anno = $_ ;
		     $inforepeatregion->{$CNT}->{START} = $start;
		     $inforepeatregion->{$CNT}->{END} = $end;
	     }
	     if(/^\s*Modified/){
	         my (@l) = split ; 
		     my $start = $l[2];
		     my $anno = $_ ;
		     $infomodified->{$CNT}->{START} = $start ;
	     }
	 	    
    }
	return ($inforepeatregion,$infomodified) ;
}


sub util_SortTwoStrings{
	my ($a,$b) = @_ ; 
	if($a lt $b){
		return ($a, $b);
	}
	else{
		return ($b, $a);
	}
    
}

sub util_ParseAAGroups{
	my ($in) = @_ ; 
    my $ifh = util_read($in);
	my $grp = 0 ; 
	my $info ={};
	my $map3to1 ={};
	my $map1to3 ={};
	my @grps ;
	while(<$ifh>){
        next if(/^\s*$/);
		#print ; 
	    my (@l) = split ; 
		$grp++;
		my $NM ; 
	    while(@l){
	 	   my $single = shift @l ; 
		   if(!defined $NM){
		       $NM = $single ;
			   push @grps, $single ;
		   }
		   #$info->{$single} = $grp ; 
		   $info->{$single} = $NM ; 
		   #print "$single \n";
	 	   shift @l ; 
	 	   shift @l ; 
	 	   my $three = shift @l ; 

		   $three =~ s/\(//;
		   $three =~ s/\)//;
		   $three = uc($three);
		   $map1to3->{$single} = $three ;
		   $map3to1->{$three} = $single ;
	    }
	}
	return ($info,\@grps,$map3to1,$map1to3);
}

sub util_PairwiseDiff{
	my ($l1,$l2) = @_ ; 
	my @l1 = @{$l1};
	my @l2 = @{$l2};
	my $N1 = @l1 ;
	my $N2 = @l2 ;
	die "$N1 $N2 not equal" if($N1 ne $N2);
	my $diff = 0 ;
	while(@l1){
	   my $a = shift @l1 ;
	   my $b = shift @l2 ;
	   $diff = $diff + abs ($a-$b);
	}
	return ($diff);
}

sub util_parseHETPDB{
	my ($infile) = @_ ;
    my $NOHET = {};
    my $YESHET = {};
    my $HET2PDB = {};
    my $HET2PDBSIZE = {};
    my $ifh = util_read($infile);
    while(<$ifh>){
         if(/^\s*NOHET/){
	          my ($nm,@v) = split ; 
		      $NOHET = util_make_table(\@v);
	     }
         elsif(/^\s*YESHET/){
	          my ($nm,@v) = split ; 
			  foreach my $v (@v){
			  	   my ($protein,$num) = split ":", $v ;
				   $YESHET->{$protein} = $num ;
			  }
	    }
		else{
	          my ($nm,$size,@v) = split ; 
		      $HET2PDB->{$nm} = (\@v);
		      $HET2PDBSIZE->{$nm} = $size ;  
		}
    }
	return ($NOHET,$YESHET,$HET2PDB,$HET2PDBSIZE);
}

sub util_parseHolo2Apo{
	my ($infile) = @_ ;
    my $info = {};
    my $ifh = util_read($infile);
    while(<$ifh>){
		my ($a,$b,$num) = split ;
		$info->{$a} = $b ; 
    }
	return $info ;
}

sub util_ParsePremonIn{
	my ($infile,$errfh) = @_ ;
    my $ifh = util_read($infile);
    my $info = {};
	my ($statusStr,$statusD,$statusPD,$statusREALSTRING,$resultstr);
	while(<$ifh>){
		print if($verbose);
		if(/^\s*STR /){
			$statusStr = 1 ;
			my ($junk,$str) = split ;
			$info->{STR} = [] if(!defined  $info->{STR});
			push @{$info->{STR}}, $str ;
		}
		if(/^\s*D /){
			$statusD = 1 ;
			my ($junk,@v) = split ;
			$info->{D} = \@v 
		}
		if(/^\s*PD /){
			$statusPD = 1 ;
			my ($junk,@v) = split ;
			$info->{PD} = \@v ;
			print "@v ......\n";
		}
		if(/^\s*REALSTRING /){
			$statusREALSTRING = 1 ;
			my ($junk,@v) = split ;
			$info->{REALSTRING} = \@v 
		}
		if(/^\s*RESULTSTYLE /){
			$resultstr = 1 ;
			my ($junk,@v) = split ;
			$info->{RESULTSTYLE} = \@v 
		}
	}
	my $status = 1 ; 
	if(!defined ($statusStr && $statusPD && $statusD && $statusREALSTRING && $resultstr)){
	    print $errfh "Could not read $infile\n" ; 
		die ;
		$status = 0 ;
	}
	close($ifh);
	return ($status,$info);
}

sub util_ReadPremonOut{
   my ($premon,$size) = @_ ;
   my $info = {};
   my $ifh = util_read($premon);
   while(<$ifh>){
          next if(/^\s*$/);
          next if(/^\s*#/);
	      my @l  = split ; 
	      my $nm = shift @l ;
	      my $len = length($nm);
	      die "len = $len $nm " if($len ne $size);
	      $info->{$nm} = \@l ;
   }
   close($ifh);
   return $info ;
}

sub util_ReadPremonOutAndGiveFasta{
   my ($premon,$size) = @_ ;
   my $info = {};
   my $ifh = util_read($premon);
   while(<$ifh>){
          next if(/^\s*$/);
          next if(/^\s*#/);
	      my @l  = split ; 
	      my $nm = shift @l ;
	      my (@aacodes) = split "", $nm ;
	      my $len = length($nm);
	      die "len = $len $nm " if($len ne $size);
		  foreach my $match (@l){
			    $match =~ s/\./ /g;
			    my (@numbers) = split " ", $match ;
				foreach my $idx (0..$size-1){
					my $aa = $aacodes[$idx];
					my $num = $numbers[$idx];
					$info->{$num} = $aa ;
				}
				
		  }  
   }

   my $fasta = "";
   foreach my $i (sort {$a <=> $b} keys %{$info}){
		my $s = $info->{$i};
		$fasta = $fasta . $s ;
   }
   close($ifh);

   return $fasta ;
}

## onus is onto caller to delete seq file
sub util_GetFasta{
    my ($pdb1,$name) = @_ ; 
    my $pdb = new PDB();
    $pdb->ReadPDB($pdb1);

    my $outfile = util_getTmpFile();
    my $ofh = util_write($outfile);
    my $seq = $pdb->WriteFasta("$name",$ofh);
    print STDERR "Wrote fasta in $outfile\n" if($verbose);

	close($ofh);
    return ($seq,$outfile) ;
}

sub util_ScorePD{
	my ($l1,$l2) = @_ ;
	my @l1 = @{$l1};
	my @l2 = @{$l2};
	my $N = @l1 ;
	die if($N ne @l2);
	my @listofdiffs ;
	foreach my $i (0..$N-1){
		my $p1 = $l1[$i];
		my $p2 = $l2[$i];
	    my ($realdiff,$thr);
		my $absdiff = util_format_float(abs($p1 - $p2),3);
    
		my $ABSREF = abs($p1);
		my $ABSQUERY = abs($p2);
		$thr = 150 ;
		    
		if($ABSREF > 400){
		         $thr = 170  ;
		}
		elsif($ABSREF > 300){
		         $thr = 150  ;
		}
		elsif ($ABSREF > 200){
		         $thr = 125  ;
		}
		#$thr = 75   if($ABSREF > 100 && $ABSREF < 131 && $ABSQUERY < $ABSREF);
		$thr = 75   if($ABSREF > 100 || $ABSQUERY > 100);
    
    
		if(abs($p1) < 99.9 && abs($p2) < 99.9){
			$thr = "bothbelow100";
			$realdiff = 0 ; 
		 }
		 else{
			  $realdiff = $absdiff < $thr ? 0 : $absdiff;
		 }
		 push @listofdiffs, $realdiff ;

	}

    my @resultssorted = sort { $a <=> $b } @listofdiffs ; 
	my $indexformax = @resultssorted ;
	my $maxdiff = int($resultssorted[$indexformax - 1]);

	return ($maxdiff);
}

sub util_ScoreDistance{
	my ($l1,$l2,$MAXALLOWEDDISTDEV) = @_ ;
	my @l1 = @{$l1};
	my @l2 = @{$l2};
	my $N = @l1 ;
	die if($N ne @l2);
	my $sumsq =  0 ;
	my $maxdiff  = 0 ;
	foreach my $i (0..$N-1){
		my $a = $l1[$i];
		my $b = $l2[$i];

		my $diff = $a - $b;
		my $diffsq = $diff * $diff ;

		my $absdiff = util_format_float(abs($a - $b),1);
		$maxdiff = $absdiff if($absdiff > $maxdiff);

		$sumsq = $sumsq + $diffsq ;
	}
	my $rmsd = util_format_float(sqrt($sumsq/$N),1);
	if($maxdiff  > $MAXALLOWEDDISTDEV){
		$rmsd = $rmsd + $maxdiff ;
	}
	return ($rmsd,$maxdiff);
}

sub parseFPocketFile{
    my ($protein,$pdb1,$addtopdb,$infile) = @_ ;
    my $info = {};
    my $infoperatom = {};
    my $ifh = util_read($infile);
	my $volume  = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
         chop ;

		 if(/Real volume/){
		     ($volume) = (/Real volume \(approximation\)\s*:(.*)/);
		 }
    
	    if(/^ATOM/){
	       my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($_);

		   $atomnm =~ s/\s*//g;
	       my $I = $resname . "/" . $resnum ;
	       my $ATOM = $I . "/" . $atomnm ;
	       $info->{$resnum} = $I ; 
		   $pdb1->AddPocket($resnum) if($addtopdb);
		   $ATOM =~ s/\s*//g;
	       $infoperatom->{$ATOM} = 1  ; 
		   #print "$ATOM \n";
	    }
    }
	return ($info, $infoperatom,$volume);
}

sub parseCastp{
    my ($infile,$subtract) = @_ ;
    my $info = {};
    my $infoperatom = {};
    my $ifh = util_read($infile);
	my $volume  = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
         chop ;

	    if(/^ATOM/){
	       my ($atomstr , $serialnum , $atomnm , $alt_loc , $resname , $chainId , $resnum , $codeforinsertion , $x , $y , $z ) = util_ReadLine($_);
		   my @l = split ;
		   my $N = @l ;
		   my $POCKNUMER = $l[$N-2];

		   $atomnm =~ s/\s*//g;
	       my $I = $resname . "/" . $resnum ;
	       my $ATOM = $I . "/" . $atomnm ;
		   if(!defined $info->{$POCKNUMER}){
		        $info->{$POCKNUMER} = {};
		   }
	       $info->{$POCKNUMER}->{$resnum} = $I ; 

		   #$ATOM =~ s/\s*//g;
	       #$infoperatom->{$ATOM} = 1  ; 
		   #print "$ATOM \n";
	    }
    }

	my $biggestnumber ;
	foreach my $POCKNUMER ( sort { $a <=> $b} keys %{$info}){
		my @res = (keys %{$info->{$POCKNUMER}});
		my $N = @res ;
		#print "$POCKNUMER has $N res - castp \n";
		$biggestnumber = $POCKNUMER ;

	}
	my $infobiggest = $info->{$biggestnumber - $subtract};


	return ($info,$infobiggest);
}

sub util_AlignAndMatchRemaining{

   my ($p1,$p2,$infile,$maxdist) = @_ ;
   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();

   
   my $file1 = "$PDBDIR/$p1.pdb";
   my $file2 = "$PDBDIR/$p2.pdb";
   my $pdb1 = new PDB();
   $pdb1->ReadPDB($file1);
   my $pdb2 = new PDB();
   $pdb2->ReadPDB($file2);

   my ($atoms1,$atoms2) = pymolin_getResultsLine($infile,$pdb1,$pdb2);
   
   return util_AlignAndMatchRemainingAtoms($pdb1,$pdb2,$atoms1,$atoms2,$maxdist);
   	
   
}

sub util_AlignAndMatchRemainingAtoms{

   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
   
   my ($pdb1,$pdb2,$atoms1,$atoms2,$maxdist) = @_ ;

   #print "Aligning geom_Align3PointsToXYPlane \n";
   my ($done1,$remainingAtoms1) = MyGeom::geom_Align3PointsToXYPlane($pdb1,$atoms1,$verbose);
   my ($done2,$remainingAtoms2) = MyGeom::geom_Align3PointsToXYPlane($pdb2,$atoms2,$verbose);
   
   
   my $i = shift(@{$remainingAtoms1}) or die "Need an extra residue";
   my $JJ = shift(@{$remainingAtoms2}) or die "Need an extra residue";
   my $NM2 = $JJ->GetName();
   
   	my @tmp1 = (@{$done1}, $i);
   	$i->Print() if($verbose);
   	my @atomlist ;
   	push @atomlist, $i ;
   	my ($results,$combined) = $pdb2->GetNeighbourHoodAtom(\@atomlist,$maxdist);
   	print STDERR  "Atoms close to this one\n" if($verbose);
   	my $sort ;
       foreach my $j (@{$combined}){
   		my $NM = $j->GetName();
   		if($NM eq $NM2){
		   my $d = util_format_float($i->Distance($j),1) ;
   		   print "$d $NM ;;;;;;;;;;\n" if($verbose);
   		   return $d; ;
   		}
   	}
   
   	return 1000 ;
   	
   
}

sub util_ParseDSSP{
   my ($protein,$pdb1,$dsspinfile,$what,$force) = @_ ;
   print "Parsing DSSP for $protein for $what\n";

   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();

   my $finaloutfile = "$HELIXDIR/$protein.$what" . "LIST";
   if(!defined $force && -e $finaloutfile && ! -z $finaloutfile){
   	   print "$finaloutfile already exist and not zero, so dssp files have already been processed\n";
       return $finaloutfile ;
   }

   system ("touch HTH") if(! -e "HTH");
   my $ofhhth = util_append("HTH");

   ## make the dssp file is not there 
   if(defined $force ||  ! -e $dsspinfile){
	   print "\tWriting DSSP $dsspinfile\n";
   	   system("mkdssp -i $PDBDIR/$protein.pdb -o $DSSP/$protein.dssp");
   }
   my $ofhPDBHELIX = util_write($finaloutfile);
   
   my $ifh = util_read($dsspinfile);
   while(<$ifh>){
        last if(/#  RESIDUE AA STRUCTURE BP1 BP2/);
   }
   
   my @arr = ();
   foreach my $i (0...10000){
   	push @arr, 0 ;
   }
   
   my $final = 0 ;
   while(<$ifh>){
   	    next if(/^\s*$/);
        chop ;
   	     my (@l) = split ; 
   	     my $n = $l[1];
   	     $final= $n;
   	     my $h = $l[4];
		 if($what eq "HELIX"){
   	         if($h eq "H" || $h eq "G" || $h eq "I"){
   	 	        $arr[$n] = 1 ;
   	         }
		 }
		 else{
   	         if($h eq "B" || $h eq "E"){
   	 	        $arr[$n] = 1 ;
   	         }
		 }
   }
   $final++;
   #print "Setting 2 value to $final\n";
   $arr[$final]= 2;
   
   my $idx = 0 ;
   my $prevend = 0 ;
   my $prevstart = 0 ;
   my $start = 0 ;
   my $helixnumber = 0 ;
   while($idx < 10000){
      ($idx,$start) = __GetNextStretch(\@arr,$idx,$helixnumber,$protein,$pdb1,$ofhPDBHELIX,$what);
      last if($idx eq -1);
      my $diff = $start - $prevend ;
      if($diff < 5 && $helixnumber){
      	   my $prevhnumber = $helixnumber -1 ;
   	       my $QQQ = "$protein.${what}$prevhnumber";
   	       my $PPP = "$protein.${what}$helixnumber";
           print $ofhhth "$QQQ $PPP $protein $diff $prevstart $idx   \n";
      }
   
      $prevend = $idx; 
      $prevstart = $start; 
      $idx++;
      $helixnumber++;
   
   }
   print "\tSplit in $helixnumber of $what\n";
   
   close($ifh);
   return $finaloutfile ;
}
   
sub __GetNextStretch{
   	my ($l,$idx,$hnumber,$protein,$pdb1,$ofhPDBHELIX,$what) = @_;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR) = util_SetEnvVars();
   	my @l = @{$l};
   	my @ret ; 
   	my $start;
   	my $end;
   	foreach my $i ($idx...10000){
   		my $x = $l[$i];
   	    return -1 if($x eq 2);
   		if($x eq 1){
   			$start = $i;
   			last ;
   		}
   	}
   	return -1 if(!defined $start) ;
   
   	foreach my $i ($start...10000){
   		my $x = $l[$i];
   		if($x eq 0 || $x eq 2){
   			$end = $i -1 ;
   			last ;
   		}
   	}
   	    die "start = $start " if(!defined $end);
   		my $tableofres = {};
   		foreach my $i ($start..$end){
   			$tableofres->{$i} = 1 ;
   		}
   	    my $PPP = "$HELIXDIR/$protein.${what}$hnumber";
   	    my $OOO = "$PPP.pdb";
   		print $ofhPDBHELIX "$protein.${what}$hnumber\n";
        my $ofhhelix= util_write($OOO);
   		print "Writing hnumber $hnumber to file $OOO : residues $start to $end \n" if($verbose);
        my $pdb = "$PDBDIR/$protein.pdb";
   		$pdb1->ReadPDBAndWriteSubset($pdb,$tableofres,$ofhhelix);
   		close($ofhhelix);
   	    print "idx=$idx start=$start end=$end \n" if($verbose);
		if($what eq "HELIX"){
		   $pdb1->AddHelix($hnumber,$start,$end);
		}
		else{
		   $pdb1->AddBeta($hnumber,$start,$end);
		}
   	return ($end,$start) ;
}

sub util_helixWheel{
   my ($protein,$pdb1,$donefh,$ofhcircular,$writetex,$justseqeunce) = @_ ; 
   my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();
   my $tableATOM = {};
   $tableATOM->{LYS} = "NZ";
   $tableATOM->{ARG} = "NH1";
   $tableATOM->{HIS} = "ND1";
   $tableATOM->{SER} = "OG";
   $tableATOM->{THR} = "OG1";
   $tableATOM->{CYS} = "SG";
   $tableATOM->{ASN} = "OD1";
   $tableATOM->{GLN} = "NE2";
   $tableATOM->{TYR} = "CB";
   $tableATOM->{PHE} = "CB";
   $tableATOM->{TRP} = "CB";
   $tableATOM->{ASP} = "OD1";
   $tableATOM->{GLU} = "OE1";
   $tableATOM->{GLY} = "CA";
   $tableATOM->{VAL} = "CB";
   $tableATOM->{ALA} = "CB";
   $tableATOM->{LEU} = "CB";
   $tableATOM->{MET} = "CB";
   $tableATOM->{ILE} = "CB";
   $tableATOM->{PRO} = "CB";

   my $outfile = "$HELIXDIR/$protein.wheel.tikz.tex";

   my $ofhpymol= util_write("pymol.helix.in");
   system ("touch $HELIXDIR/commands");
   my $ofhcommands = util_append("$HELIXDIR/commands");

    
    my $colortable = {};
    my $value = {};
    my $chargedtable = {};
    my @NP = qw(G P A V L I M );
    foreach my $k (@NP){
    	$colortable->{$k} = "red";
    	$value->{$k} = 100 ;
    }
    
    my @pos = qw (H K R);
    foreach my $k (@pos){
    	$chargedtable->{$k} = 1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 100 ;
    }
    
    my @neg = qw (E D); 
    foreach my $k (@neg){
    	$chargedtable->{$k} = -1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 60 ;
    }
    my @amide = qw(Q N );
    foreach my $k (@amide){
    	$colortable->{$k} = "blue";
    	$value->{$k} = 10 ;
    }
    
    my @misc = qw(F Y W C S T);
    foreach my $k (@misc){
    	$colortable->{$k} = "red";
    	$value->{$k} = 75 ;
    }
    
    my $HYDROVAL = {};
    ## from "Computer programs to identify and classify amphipathic a helical domains"
    $HYDROVAL->{MET} = 0.975;
    $HYDROVAL->{ILE} = 0.913;
    $HYDROVAL->{LEU} = 0.852;
    $HYDROVAL->{VAL} = 0.811;
    $HYDROVAL->{CYS} = 0.689;
    $HYDROVAL->{ALA} = 0.607;
    $HYDROVAL->{THR} = 0.525;
    $HYDROVAL->{GLY} = 0.484;
    $HYDROVAL->{SER} = 0.402;
    $HYDROVAL->{HIS} = 0.333;
    $HYDROVAL->{PRO} = 0.239;
    $HYDROVAL->{PHE} = 1.036;
    $HYDROVAL->{TRP} = 0.668;
    $HYDROVAL->{TYR} = 0.137;
    $HYDROVAL->{GLN} = -0.558;
    $HYDROVAL->{ASN} = -0.701;
    $HYDROVAL->{GLU} = -1.396;
    $HYDROVAL->{LYS} = -1.518;
    $HYDROVAL->{ASP} = -1.600;
    $HYDROVAL->{ARG} = -2.233;
    
    my $hydroval = {};
    foreach my $k (keys %{$HYDROVAL}){
    	my $single = $pdb1->GetSingleLetter($k);
    	$hydroval->{$single} = $HYDROVAL->{$k};
    }
    
    
    
   my $ofh ;
   $ofh = util_write("$outfile") if($writetex);
    print $ofh "% for protein $protein \n" if($writetex);
        
        my @res = $pdb1->GetResidues();
        my $N = @res;
    	my $cnt = 0 ;
    
    	my @singles ;
    	my @numbers ;
    	if(! defined $justseqeunce){
    	     foreach my $r (@res){
    	         die if($r->GetAtomStr() eq "HETATM");
    	         die if($r->GetName() eq "HOH");
    		     my $single =  $r->PrintSingleLetter($pdb1);
    		     push @singles, $single ;
    		     my $num = $r->GetResNum();
    		     push @numbers, $num ;
    	     }
    	}
    	else{
            my $ifh = util_read($justseqeunce);
    		my $sequence = "";
            while(<$ifh>){
                 next if(/^\s*$/);
                 next if(/^\s*>/);
                 chop ;
    			 $sequence = $sequence . $_ ;
            }
    		$sequence =~ s/\s*//g ;
    		@singles = split "", $sequence ;
    		my $N = @singles;
    		foreach my $i (1..$N){
    			push @numbers, $i ; 
    		}
            close($ifh);
    	}
    
    
    	#my @lol ;
    	#$lol[0] = []; $lol[1] = []; $lol[2] = []; $lol[3] = [];
        my $done = {};
        my $initval = -270 ; 
        my $rad = 5 ; 
    	my $loopcnt = 0 ;
    
    	my @centre = qw (0 0 0);
    	my $finalveccharged; 
    	my $finalvechydro; 
    
    	my $TOTALCHARGED = 0 ;
    	my $TOTALPOS = 0 ;
    
    
    	my $bins = {};
    	my $TOTALNUM = 0 ; 
    	my $circularseq = {};
    	my $circularval = {};
    	while(@singles){
    		$TOTALNUM++;
    		my $single =  shift @singles ;
    		my $num = shift @numbers ;
    
    
    		my $idxcnt = $loopcnt % 4 ;
    		if(!defined $bins->{$idxcnt}){
    			$bins->{$idxcnt} = [];
    		}
    
    
    		my $x = $pdb1->{SINGLE2THREE}->{$single} ;
    		my $what = $tableATOM->{$x};
    		#print "$x/$num/$what $idxcnt $x\n";
    		my $SSSS = $single . $num;
    
    		push @{$bins->{$idxcnt}}, "$x/$num/$what" ;
    
            my $val = $initval - $loopcnt * 100 ; 
    		
    		my $cval = (deg2rad($val ));
    		while( exists $circularval->{$cval}){
    			$cval = $cval + 0.01;
    		}
    		$circularseq->{$SSSS} = $cval ; 
    		$circularval->{$cval} = 1 ;
    
            $x = util_format_float($rad * cos(deg2rad($val)),1) + 0 ;
            my $y = util_format_float($rad * sin(deg2rad($val)),1) + 0 ;
            my $str = "$x.$y";
    
    		my @thispoint ; 
    		push @thispoint, $x ;
    		push @thispoint, $y ;
    		push @thispoint, 0 ;
    
    		    my $VEC = MyGeom::MakeVectorFrom2Points(\@centre,\@thispoint) ;
    		    my $tmpvec = $VEC->norm * $hydroval->{$single} ;
    		    if(!defined $finalvechydro){
    			    $finalvechydro = $tmpvec ;
    		    }
    		    else{
    			    $finalvechydro = $finalvechydro + $tmpvec ;
    		    }
    
    		if(exists $chargedtable->{$single}){
    			$TOTALCHARGED++;
    			my $factor = $chargedtable->{$single} ;
    			if($factor eq 1){
    				$TOTALPOS++;
    			}
    		    #my $tmpvec = $VEC->norm * $hydroval->{$single} *$factor;
    		    my $tmpvec = $VEC->norm * $factor;
    		    if(!defined $finalveccharged){
    			    $finalveccharged = $tmpvec ;
    		    }
    		    else{
    			    $finalveccharged = $finalveccharged + $tmpvec ;
    		    }
    		}
    
            while(exists $done->{$str}){
                #print "$str existed\n";
                $y = $y + 1 ;	
                $str = "$x.$y";
            }
            $done->{$str} = 1 ;
    		my $color = $colortable->{$single};
    		my $VAL = $value->{$single};
    		#my $VAL = 100 ;
            print $ofh "\\node[fill=$color!$VAL,draw=blue,very thick] (FinalNode) at ($x, $y) {$single$num} ; \n" if($writetex);
    		print $ofhpymol "select A, /PDBB//A/$num/CA\n";
    		print $ofhpymol "# $single\n";
    		print $ofhpymol "color $color, A\n";
    
    
    		$loopcnt++;
    		#push @{$lol[$idxcnt]}, $r ;
    	}
    
        my $fastacircular = "$protein.circular.fasta";
        print $ofhcircular ">$fastacircular\n";
    	foreach my $k (sort { $circularseq->{$a} <=> $circularseq->{$b}} keys %{$circularseq}){
    		my $v = $circularseq->{$k};
    		my ($firstletter) = ($k =~ /(.)/);
    		print $ofhcircular "$firstletter";
    
    	}
    		print $ofhcircular "\n";
    
    
    	#### FOR CLASP #########
        my ($REALPDB) = ($protein =~ /(....)/) ;
    	$REALPDB = $REALPDB . "A";
    	foreach my $k (sort  keys %{$bins}){
    		my $l = $bins->{$k} ;
    		my $CONCATSTR = "";
    		my $cntconcat = 0 ;
    		foreach my $i (@{$l}){
    			$cntconcat++;
    			$CONCATSTR = $CONCATSTR . " $i ";
    		}
    		if($cntconcat > 1){
    		    print $ofhcommands "echo $CONCATSTR > ! IN \n";
    		    print $ofhcommands "printPairwiseOneAfterAnother.pl -out out -c \$CONFIGGRP -ra 222 -pr $REALPDB -in IN -helix $protein \n";
    		}
    	}
    
    my $Lcharged = 0  ;
    my $cx = 0 ;
    my $cy = 0 ;
    if(defined $finalveccharged){
        $Lcharged = util_format_float($finalveccharged->length,1) ;
        ($cx) = $finalveccharged->x * 1;
        ($cy) = $finalveccharged->y * 1;
    }
    my $Lhydro = util_format_float($finalvechydro->length,1) ;
    
    my ($hx) = $finalvechydro->x * 1;
    my ($hy) = $finalvechydro->y * 1;
    print $ofh "\\node[fill=black!10,draw=yellow,very thick] (HPH) at ($hx, $hy) {\$\\mu_{h}\$} ; \n" if($writetex);
    #print $ofh "\\node[fill=white!100,draw=blue,very thick] (Ch) at ($cx, $cy) {Ch} ; \n" ;
    ;
    print $ofh "\\node[fill=blue!100,draw=blue] (centre) at (0, 0) {} ;\n" if($writetex);
    
      print $ofh "\\begin{scope}[every path/.style=line]\n" if($writetex);
       print $ofh "\\path          (centre) -- node [near end] {$Lhydro} (HPH);\n" if($writetex);
       #print $ofh "\\path          (centre) -- node [near end] {$Lcharged} (Ch);\n" if($writetex);
    print $ofh "\\end{scope}\n" if($writetex);
    
    my $PERCENTPOS = -1 ;
    if($TOTALCHARGED){
        $PERCENTPOS = util_format_float($TOTALPOS/$TOTALCHARGED,1);
    }
    
    print "TIKZ o/p written to $outfile Lcharged = $Lcharged and Lhydro = $Lhydro PERCENTPOS = $PERCENTPOS TOTALNUM = $TOTALNUM \n" if($verbose);
    
    
    #print $ofhdata "$protein $N $Lhydro $Lcharged $PERCENTPOS \n";
    print $donefh "$protein $N $Lhydro $PERCENTPOS $TOTALCHARGED \n";
    close($ofh) if($writetex);
    system("unlink in ; ln -s $outfile in");
}

sub util_RunHelixWheel{
	my ($protein,$helixlist,$force) = @_;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();

	my $ifh = util_read($helixlist);
	my $donefile = "$HELIXDIR/$protein.HELIXVALUES";
	my $circularfasta = "$HELIXDIR/$protein.HELIXCIRCULARFASTA";
	if(!defined $force && -e $donefile){
		 print "$donefile exists, so the helix values have already been extracted\n";
	     return ;
	}
	print "Writing HELIXVALUES to $donefile\n";
	my $donefh = util_write($donefile);
	my $ofhcircular = util_write($circularfasta);
	while(<$ifh>){
		next if(/^\s*$/);
		chomp ;
		my ($f) = split ;
        my $pdb = "$HELIXDIR/$f.pdb";

        my $pdb1 = new PDB();
        $pdb1->ReadPDB($pdb);
		my $writetex = 0 ;
        util_helixWheel($f,$pdb1,$donefh,$ofhcircular,$writetex);

	}
	
}

sub util_GetDisulphide{
	my ($protein,$pdb1,$plsprint) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();
    my @res = $pdb1->GetResidues();
    my $N = @res;

	my $CNT = 0 ;
    while(@res){
		
		my $res1 = shift @res ;
		foreach my $res2(@res){
		next if($res1->GetAtomStr() eq "HETATM");
		next if($res2->GetAtomStr() eq "HETATM");
		next if($res1->GetName() eq "HOH");
		next if($res2->GetName() eq "HOH");

		my $bothcys = util_IsResisdueThisType($res1,"CYS") && util_IsResisdueThisType($res2,"CYS") ;
		next if(!$bothcys);


		my $num1 = $res1->GetResNum();
		my $num2 = $res2->GetResNum();
		my $CA1 = $pdb1->GetAtomFromResidueAndType($res1->GetResNum(),"SG");
		my $CA2 = $pdb1->GetAtomFromResidueAndType($res2->GetResNum(),"SG");

        my $distCA = util_format_float($pdb1->DistanceAtoms($CA1,$CA2),1);
		if($distCA < 2.4){
		       print "dist = $distCA $num1 $num2 \n" if($plsprint);
			   $CNT++;
			   $pdb1->AddDisulphide($num1,$num2);
		}

         }
	}


	return $CNT;

}

sub util_NeedleFiles{
	my ($outfile,$f1,$f2,$arg) = @_ ;
	my $end = " -endopen 1 -endweight 1 -endextend 0.2";
	my $execinit = ParseArgFile($arg);
	#print "INIT = $execinit \n";
    #$exec = "needle -gapopen 25  -gapex 0.5 -ou $outfile $file1 $file2  2> /dev/null ";
    my $exec = $execinit . "-ou $outfile $f1 $f2  2> /dev/null "; 
	system("touch $outfile");
	system($exec);

   my ($iden,$simi); 
   my $ifh = util_read($outfile);
   while(<$ifh>){
        next if(/^\s*$/);
		if(/Identity/){
			($iden) = (/.*\((.*)\)/);
			$iden =~ s/\%//;
		}
		if(/Similarity/){
			($simi) = (/.*\((.*)\)/);
			$simi =~ s/\%//;
		}
		
   }
	return ($iden,$simi);
}

sub util_NeedleSeq{
	my ($outfile,$seq1,$seq2,$arg,$pdb1name,$pdb2name) = @_ ;
		
	my $f1 = (util_getTmpFile());
	my $f2 = (util_getTmpFile());
	my $ofh1 = util_write($f1);
	my $ofh2 = util_write($f2);

    print $ofh1 ">$pdb1name";
    print $ofh1 "$seq1\n";

    print $ofh2 ">$pdb2name\n";
    print $ofh2 "$seq2\n";

	close($ofh1);
	close($ofh2);

	return util_NeedleFiles($outfile,$f1,$f2,$arg);
}

sub util_NeedlePDBNamesFromSeq{
	my ($outfile,$pdb1name,$pdb2name,$arg) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();
    my $pdb1 = "$PDBDIR/$pdb1name.pdb";
    my $pdb2 = "$PDBDIR/$pdb2name.pdb";

	my ($seq1,$f1) = util_GetFasta($pdb1,$pdb1name);
	my ($seq2,$f2) = util_GetFasta($pdb2,$pdb2name);

	
	my $N1 = length($seq1);
	my $N2 = length($seq2);
	if(!$N1 || !$N2){
	    print "empty sequence in one $N1 or $N2\n";
	    return (-1,-1)  ;
	}



   my ($iden,$simi)=	 util_NeedleFiles($outfile,$f1,$f2,$arg,$pdb1name,$pdb2name);
   unlink($f1);
   unlink($f2);
   return ($iden,$simi);
}



sub util_NeedlePDBNamesFromFASTADIR{
	my ($outfile,$pdb1name,$pdb2name,$arg) = @_ ;
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP) = util_SetEnvVars();
    my $f1 = "$FASTADIR/$pdb1name.ALL.1.fasta";
    my $f2 = "$FASTADIR/$pdb2name.ALL.1.fasta";

	if(! -e $f1 || ! -e $f2){
		die "$f1 or $f2 does not exist \n";
	}

	return util_NeedleFiles($outfile,$f1,$f2,$arg);
}

sub ParseArgFile{
	my ($arg) = @_ ;
	my $exec  = "needle " ;
    my $ifh = util_read($arg);
    while(<$ifh>){
         next if(/^\s*$/);
	     my ($nm,$junk) = split ; 
		 $exec = $exec . "-". $nm . " $junk " ;
    }
	return $exec ; 
}

sub util_ParseBlastPW{
	my ($infile) = @_ ;
	my $ifh = util_read($infile);
	while(<$ifh>){
		if(/Expect/){
			my ($val) = (/Expect = (.*),/);
			if(!defined $val){
			($val) = (/Expect = (.*)/);
			}
			return $val ;
		}
	}
	close($ifh);
	return -1 ; 
}
sub utils_parseBlastOut{
	my ($infile,$cutoff) = @_ ;
    my $ifh = util_read($infile);
    my @l ; 
    my $push = 0 ; 
    while(<$ifh>){
         next if(/^\s*$/);
	     chomp ;
	     if(/Sequences producing significant alignments/){
		    $push = 1 ; 
		    next ;
	     }
	     if($push && !/(chromosome|hypothe|unnamed|uncharacterized)/i){
	         push @l, $_ 
	     }
    }
    my $v = $l[0] ;
    die "$infile problem" if(!defined $v);
    
	if(defined $cutoff){
         my @lll = split " ", $v ;
         my $N = @lll ;
         my $evalue = $lll[$N-1];
         if($evalue > $cutoff){
	         print "Warning; $infile has evalue $evalue more than cutoff $cutoff\n";
         }
	}
	return $v ;
}
sub util_Blastout_isrev{
   my ($infile) = @_ ;
   my $ifh = util_read($infile);
   #my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();

   my $info = {};
   my ($subject);
   my $isrev = 0 ;
   while(<$ifh>){
        next if(/^\s*$/);
	    chomp ;
        if(/>/){
	 	   s/>//;
		   ($subject) = split ;
	    }
	    elsif(/Sbjct/){
	 	   my (@l) = split ;
		   my $s = $l[1];
		   my $e = $l[3];
		   if($s > $e){
		   	 $isrev = 1 ;
		   }
		   last ;
	    }
   }
   return ($subject,$isrev);
}


sub util_ConcatTwoStringsWithCommonBeginAndEnd{
	my ($A,$B,$N) = @_;

	foreach my $i (5..$N){
		my $Ax = util_findLength($A,$i,"end");
		my $Bx = util_findLength($B,$i,"begin");
		if($Ax eq $Bx){
			$A =~ s/$Ax//;
			my $ret = $A . $B ;
			return $ret ;
		}
	}
	return -1 ;

}

sub util_findLength{
	my ($seq,$size,$what) = @_ ;

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
	elsif($size eq 16){
	   ($begin) = ($seq =~/^(................)/);
	   ($end) = ($seq =~/(................)$/);
	}

	elsif($size eq 17){
	   ($begin) = ($seq =~/^(.................)/);
	   ($end) = ($seq =~/(.................)$/);
	}
	elsif($size eq 18){
	($begin) = ($seq =~/^(..................)/); 
	($end) = ($seq =~/(..................)$/); 
	}
elsif($size eq 19){
($begin) = ($seq =~/^(...................)/); 
($end) = ($seq =~/(...................)$/); 
}
elsif($size eq 20){
($begin) = ($seq =~/^(....................)/); 
($end) = ($seq =~/(....................)$/); 
}
elsif($size eq 21){
($begin) = ($seq =~/^(.....................)/); 
($end) = ($seq =~/(.....................)$/); 
}
elsif($size eq 22){
($begin) = ($seq =~/^(......................)/); 
($end) = ($seq =~/(......................)$/); 
}
elsif($size eq 23){
($begin) = ($seq =~/^(.......................)/); 
($end) = ($seq =~/(.......................)$/); 
}
elsif($size eq 24){
($begin) = ($seq =~/^(........................)/); 
($end) = ($seq =~/(........................)$/); 
}
elsif($size eq 25){
($begin) = ($seq =~/^(.........................)/); 
($end) = ($seq =~/(.........................)$/); 
}
elsif($size eq 26){
($begin) = ($seq =~/^(..........................)/); 
($end) = ($seq =~/(..........................)$/); 
}
elsif($size eq 27){
($begin) = ($seq =~/^(...........................)/); 
($end) = ($seq =~/(...........................)$/); 
}
elsif($size eq 28){
($begin) = ($seq =~/^(............................)/); 
($end) = ($seq =~/(............................)$/); 
}
elsif($size eq 29){
($begin) = ($seq =~/^(.............................)/); 
($end) = ($seq =~/(.............................)$/); 
}
elsif($size eq 30){
($begin) = ($seq =~/^(..............................)/); 
($end) = ($seq =~/(..............................)$/); 
}
elsif($size eq 31){
($begin) = ($seq =~/^(...............................)/); 
($end) = ($seq =~/(...............................)$/); 
}
elsif($size eq 32){
($begin) = ($seq =~/^(................................)/); 
($end) = ($seq =~/(................................)$/); 
}
elsif($size eq 33){
($begin) = ($seq =~/^(.................................)/); 
($end) = ($seq =~/(.................................)$/); 
}
elsif($size eq 34){
($begin) = ($seq =~/^(..................................)/); 
($end) = ($seq =~/(..................................)$/); 
}
elsif($size eq 35){
($begin) = ($seq =~/^(...................................)/); 
($end) = ($seq =~/(...................................)$/); 
}
elsif($size eq 36){
($begin) = ($seq =~/^(....................................)/); 
($end) = ($seq =~/(....................................)$/); 
}
elsif($size eq 37){
($begin) = ($seq =~/^(.....................................)/); 
($end) = ($seq =~/(.....................................)$/); 
}
elsif($size eq 38){
($begin) = ($seq =~/^(......................................)/); 
($end) = ($seq =~/(......................................)$/); 
}
elsif($size eq 39){
($begin) = ($seq =~/^(.......................................)/); 
($end) = ($seq =~/(.......................................)$/); 
}
elsif($size eq 40){
($begin) = ($seq =~/^(........................................)/); 
($end) = ($seq =~/(........................................)$/); 
}
elsif($size eq 41){
($begin) = ($seq =~/^(.........................................)/); 
($end) = ($seq =~/(.........................................)$/); 
}
elsif($size eq 42){
($begin) = ($seq =~/^(..........................................)/); 
($end) = ($seq =~/(..........................................)$/); 
}
elsif($size eq 43){
($begin) = ($seq =~/^(...........................................)/); 
($end) = ($seq =~/(...........................................)$/); 
}
elsif($size eq 44){
($begin) = ($seq =~/^(............................................)/); 
($end) = ($seq =~/(............................................)$/); 
}
elsif($size eq 45){
($begin) = ($seq =~/^(.............................................)/); 
($end) = ($seq =~/(.............................................)$/); 
}
elsif($size eq 46){
($begin) = ($seq =~/^(..............................................)/); 
($end) = ($seq =~/(..............................................)$/); 
}
elsif($size eq 47){
($begin) = ($seq =~/^(...............................................)/); 
($end) = ($seq =~/(...............................................)$/); 
}
elsif($size eq 48){
($begin) = ($seq =~/^(................................................)/); 
($end) = ($seq =~/(................................................)$/); 
}
elsif($size eq 49){
($begin) = ($seq =~/^(.................................................)/); 
($end) = ($seq =~/(.................................................)$/); 
}
elsif($size eq 50){
($begin) = ($seq =~/^(..................................................)/); 
($end) = ($seq =~/(..................................................)$/); 
}
elsif($size eq 51){
($begin) = ($seq =~/^(...................................................)/); 
($end) = ($seq =~/(...................................................)$/); 
}
elsif($size eq 52){
($begin) = ($seq =~/^(....................................................)/); 
($end) = ($seq =~/(....................................................)$/); 
}
elsif($size eq 53){
($begin) = ($seq =~/^(.....................................................)/); 
($end) = ($seq =~/(.....................................................)$/); 
}
elsif($size eq 54){
($begin) = ($seq =~/^(......................................................)/); 
($end) = ($seq =~/(......................................................)$/); 
}
elsif($size eq 55){
($begin) = ($seq =~/^(.......................................................)/); 
($end) = ($seq =~/(.......................................................)$/); 
}
elsif($size eq 56){
($begin) = ($seq =~/^(........................................................)/); 
($end) = ($seq =~/(........................................................)$/); 
}
elsif($size eq 57){
($begin) = ($seq =~/^(.........................................................)/); 
($end) = ($seq =~/(.........................................................)$/); 
}
elsif($size eq 58){
($begin) = ($seq =~/^(..........................................................)/); 
($end) = ($seq =~/(..........................................................)$/); 
}
elsif($size eq 59){
($begin) = ($seq =~/^(...........................................................)/); 
($end) = ($seq =~/(...........................................................)$/); 
}
elsif($size eq 60){
($begin) = ($seq =~/^(............................................................)/); 
($end) = ($seq =~/(............................................................)$/); 
}
elsif($size eq 61){
($begin) = ($seq =~/^(.............................................................)/); 
($end) = ($seq =~/(.............................................................)$/); 
}
elsif($size eq 62){
($begin) = ($seq =~/^(..............................................................)/); 
($end) = ($seq =~/(..............................................................)$/); 
}
elsif($size eq 63){
($begin) = ($seq =~/^(...............................................................)/); 
($end) = ($seq =~/(...............................................................)$/); 
}
elsif($size eq 64){
($begin) = ($seq =~/^(................................................................)/); 
($end) = ($seq =~/(................................................................)$/); 
}
elsif($size eq 65){
($begin) = ($seq =~/^(.................................................................)/); 
($end) = ($seq =~/(.................................................................)$/); 
}
elsif($size eq 66){
($begin) = ($seq =~/^(..................................................................)/); 
($end) = ($seq =~/(..................................................................)$/); 
}
elsif($size eq 67){
($begin) = ($seq =~/^(...................................................................)/); 
($end) = ($seq =~/(...................................................................)$/); 
}
elsif($size eq 68){
($begin) = ($seq =~/^(....................................................................)/); 
($end) = ($seq =~/(....................................................................)$/); 
}
elsif($size eq 69){
($begin) = ($seq =~/^(.....................................................................)/); 
($end) = ($seq =~/(.....................................................................)$/); 
}
elsif($size eq 70){
($begin) = ($seq =~/^(......................................................................)/); 
($end) = ($seq =~/(......................................................................)$/); 
}
elsif($size eq 71){
($begin) = ($seq =~/^(.......................................................................)/); 
($end) = ($seq =~/(.......................................................................)$/); 
}
elsif($size eq 72){
($begin) = ($seq =~/^(........................................................................)/); 
($end) = ($seq =~/(........................................................................)$/); 
}
elsif($size eq 73){
($begin) = ($seq =~/^(.........................................................................)/); 
($end) = ($seq =~/(.........................................................................)$/); 
}
elsif($size eq 74){
($begin) = ($seq =~/^(..........................................................................)/); 
($end) = ($seq =~/(..........................................................................)$/); 
}
elsif($size eq 75){
($begin) = ($seq =~/^(...........................................................................)/); 
($end) = ($seq =~/(...........................................................................)$/); 
}
elsif($size eq 76){
($begin) = ($seq =~/^(............................................................................)/); 
($end) = ($seq =~/(............................................................................)$/); 
}
elsif($size eq 77){
($begin) = ($seq =~/^(.............................................................................)/); 
($end) = ($seq =~/(.............................................................................)$/); 
}
elsif($size eq 78){
($begin) = ($seq =~/^(..............................................................................)/); 
($end) = ($seq =~/(..............................................................................)$/); 
}
elsif($size eq 79){
($begin) = ($seq =~/^(...............................................................................)/); 
($end) = ($seq =~/(...............................................................................)$/); 
}
elsif($size eq 80){
($begin) = ($seq =~/^(................................................................................)/); 
($end) = ($seq =~/(................................................................................)$/); 
}
elsif($size eq 81){
($begin) = ($seq =~/^(.................................................................................)/); 
($end) = ($seq =~/(.................................................................................)$/); 
}
elsif($size eq 82){
($begin) = ($seq =~/^(..................................................................................)/); 
($end) = ($seq =~/(..................................................................................)$/); 
}
elsif($size eq 83){
($begin) = ($seq =~/^(...................................................................................)/); 
($end) = ($seq =~/(...................................................................................)$/); 
}
elsif($size eq 84){
($begin) = ($seq =~/^(....................................................................................)/); 
($end) = ($seq =~/(....................................................................................)$/); 
}
elsif($size eq 85){
($begin) = ($seq =~/^(.....................................................................................)/); 
($end) = ($seq =~/(.....................................................................................)$/); 
}
elsif($size eq 86){
($begin) = ($seq =~/^(......................................................................................)/); 
($end) = ($seq =~/(......................................................................................)$/); 
}
elsif($size eq 87){
($begin) = ($seq =~/^(.......................................................................................)/); 
($end) = ($seq =~/(.......................................................................................)$/); 
}
elsif($size eq 88){
($begin) = ($seq =~/^(........................................................................................)/); 
($end) = ($seq =~/(........................................................................................)$/); 
}
elsif($size eq 89){
($begin) = ($seq =~/^(.........................................................................................)/); 
($end) = ($seq =~/(.........................................................................................)$/); 
}
elsif($size eq 90){
($begin) = ($seq =~/^(..........................................................................................)/); 
($end) = ($seq =~/(..........................................................................................)$/); 
}
elsif($size eq 91){
($begin) = ($seq =~/^(...........................................................................................)/); 
($end) = ($seq =~/(...........................................................................................)$/); 
}
elsif($size eq 92){
($begin) = ($seq =~/^(............................................................................................)/); 
($end) = ($seq =~/(............................................................................................)$/); 
}
elsif($size eq 93){
($begin) = ($seq =~/^(.............................................................................................)/); 
($end) = ($seq =~/(.............................................................................................)$/); 
}
elsif($size eq 94){
($begin) = ($seq =~/^(..............................................................................................)/); 
($end) = ($seq =~/(..............................................................................................)$/); 
}
elsif($size eq 95){
($begin) = ($seq =~/^(...............................................................................................)/); 
($end) = ($seq =~/(...............................................................................................)$/); 
}
elsif($size eq 96){
($begin) = ($seq =~/^(................................................................................................)/); 
($end) = ($seq =~/(................................................................................................)$/); 
}
elsif($size eq 97){
($begin) = ($seq =~/^(.................................................................................................)/); 
($end) = ($seq =~/(.................................................................................................)$/); 
}
elsif($size eq 98){
($begin) = ($seq =~/^(..................................................................................................)/); 
($end) = ($seq =~/(..................................................................................................)$/); 
}
elsif($size eq 99){
($begin) = ($seq =~/^(...................................................................................................)/); 
($end) = ($seq =~/(...................................................................................................)$/); 
}
elsif($size eq 100){
($begin) = ($seq =~/^(....................................................................................................)/); 
($end) = ($seq =~/(....................................................................................................)$/); 
}
	else{
		die "size can be only from 4 to 7";
	}


	if($what eq "begin"){
		return $begin ;
	}
	elsif($what eq "end"){
		return $end ;
	}
	else{
		die ;
	}

}

sub util_ParseWebBlast{
   my ($infile) = @_ ;
   my $ifh = util_read($infile);
   my $info = {};
   my ($length,$percent);
   my $push = 0 ; 
   my $pushident = 0 ; 
   my @STRS;
   my $fulllength  ;
   my $identities  ;
   my $str ; 
   
   my $SSSS ;
   my @LLLs ;
   my @TRSnm ;
   my @SCORES ;
   my @EVALUES ;
   my @IDENTITY ;
   while(<$ifh>){
        next if(/^\s*$/);
	    chomp ;
		if(/No hits found/){
			return ($fulllength,undef) ;
		}
        if(/Length/ && ! defined $fulllength){
		   ($fulllength) = (/Length=(\d+)/);
	    }
	    elsif(/Identities/){
	 	   ($identities,$length,$percent) = (/Identities\s*=\s*(\d+)\/(\d+)\s*\((\d+)\%\)/);
		   #$str = "$identities / $length $percent";
	 	   push @IDENTITY, $percent;
	 	   #push @STRS, "$SSSS  $SCORES[@SCORES-1]";

		   my $evalue = $EVALUES[@EVALUES-1] + 0  ;
	 	   push @STRS, "\t$SSSS\t$evalue";
	 	   push @TRSnm, $SSSS ;
		   push @LLLs, $length ;
	    }
		elsif(/Score =/){
			s/,//g;
			my @lll = split;
			push @SCORES, $lll[2];
			push @EVALUES, $lll[7];
		}
	    if(/^>/){
	 	   s/>//;
	       $SSSS = $_ ;
	    }
   }
   if(!defined $fulllength){
			return (undef) ;
	}
   return ($fulllength,\@TRSnm,\@STRS,\@LLLs,\@SCORES,\@IDENTITY,\@EVALUES);
}

sub util_extractSliceFromFasta{
   my ($str,$start,$end) = @_;
   return util_extractSliceFromFastaString($str,$start,$end);
}
sub util_extractSliceFromFastaString{
   my ($str,$start,$end) = @_;
   die if(!defined ($start && $end));
	#print "util_extractSliceFromFasta $start $end\n";

   my @l = split "", $str ; 
   my $N = @l -1 ;
   my $diff = abs ($start - $end) + 1 ;
   print "$N+1 length, asked for $diff  \n" if($verbose);
   die if($start > $end);
   

   my $retstr = "";
   foreach my $idx ($start...$end){
	   #print "$idx $start $end\n";
	   $idx = $idx -1 ;
	   $retstr = $retstr . $l[$idx] ;
   }
   return $retstr ;

}

sub util_PARSEBLAST{
	my ($infile) =@_ ;
    
    my $ifh = util_read($infile);
    my $junk ;
    my $Query ;
    my $Querylength ;
    my $Subjectlength ;
	my @SAVED ;
	my $found = 0 ;
    while(<$ifh>){
    	chomp ;
		next if(/^BLAST/);
		next if(/^\s*$/);
		push @SAVED, $_ ;
        if(/Query= /){
    		($junk,$Query) = split ;
    		next ;
    	}
        if(/Length=/){
    		s/Length=//;
			if(!defined $Querylength){
    		    $Querylength =$_ ;
			}
			else{
				## for pw only
    		    $Subjectlength =$_ ;
			}
    		next ;
    	}
    	if(/Sequences producing significant alignments/){
			$found = 1 ;
			last ;
		}

    }   

    my @alllines ;
	if(!$found){
		print "Processing pairwise\n";
		@alllines = @SAVED ;
		my $CNT = 0 ;
		my @newalllines ;
		my $subject ;
		my $subjectlen ;
	    while(@alllines){
		    $_ = shift @alllines ;
    	    if(/^\s*Subject/){
			   s/\s*//g;
			   ($subject) = (/Subject=(.*)/);
		       $_ = shift @alllines ;
			   s/\s*//g;
			   ($subjectlen) = (/Length=(.*)/);
			}
    	    if(/Score/){
    		    push @newalllines, ">$CNT";
    		    push @newalllines, "Length=$subjectlen";
				$CNT++;
			}
			push @newalllines, $_ ;
    	}

		@alllines = @newalllines ;
	}
	else{
    
        while(<$ifh>){
           chomp ;
           push @alllines, $_ ;
        }
	}

	
    my $info = [];
	while(@alllines){
		$_ = shift @alllines ;
    	if(/^\s*>/){
    		ParseSingleMatch($_,\@alllines,$info);
			$found = 1 ;
    		last ;
    	}
	}


	return ($info,$Querylength);
}


sub ParseSingleMatch{
	my ($line,$alllines,$ALLINFO) = @_ ;
	my $N = @{$alllines};
	#print "$alllines $N kkkkkkkkk\n";
	$line =~ s/>//;
	my ($name,@l) = split " " ,$line ;
	my $description = join "ZZZ", @l;

	my ($junk,$subjectlength,$blastscore,$expect,$iden_percent,$ratio,$subjectmatched) ;
	my ($querystart,$queryend) ;
	my ($subjectstart,$subjectend) ;

	while(@{$alllines}){
		## recurse
		$_ = shift @{$alllines};
	    if(/^\s*>/){
		    ParseSingleMatch($_,$alllines,$ALLINFO);
	    }
		else{
           if(/Length=/){
		       s/Length=//;
		       $subjectlength =$_ ;
	       }
		   elsif(/Score/){
		   		if(defined $querystart){
				 	$description = "UUU" if($description eq "");
			       my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
			       #print "XXXXXXX $str\n";
				   push @{$ALLINFO}, $str ;
				}
		   		undef $querystart;
		   		undef $queryend;
		   		undef $subjectstart;
		   		undef $subjectend;
		   	    ($junk,$junk,$blastscore,$junk,$junk,$junk,$junk,$expect)= split ;
			    $expect =~ s/,//;
			}
		   elsif(/Identities/){
		   	      ($junk,$junk,$ratio,$iden_percent) = split ;
				   ($junk,$subjectmatched) = split "/", $ratio;
				   $iden_percent =~ s/\%//;
				   $iden_percent =~ s/,//;
				   $iden_percent =~ s/\(//;
				   $iden_percent =~ s/\)//;
		   }
		   elsif(/Query/){
		   	  my (@l) = split ;
			  my $N= @l -1;
			  $querystart = $l[1] if(!defined $querystart);;
			  $queryend = $l[$N] ;
		   }
		   elsif(/Sbjct/){
		   	  my (@l) = split ;
			  my $N= @l -1;
			  $subjectstart = $l[1] if(!defined $subjectstart);;
			  $subjectend = $l[$N] ;
			  #print "$subjectstart $subjectend \n";
		   }


		}
	}

	$description = "UUU" if($description eq "");
	my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
	#print "XXXXXXX $str\n";
    push @{$ALLINFO}, $str ;

}




sub util_R_extractIndex{
	my ($INF,$var,$idx) = @_ ;
    my $ifh = util_read($INF);
	my @retl;
    while(<$ifh>){
         next if(/^\s*$/);
	     my (@l) = split " " ,$_ ;
	      push @retl, $l[$idx];
    }
	close($ifh);
	return @retl ;
}

sub util_get_max{
	my ($l) = @_ ;
	my @l = @{$l};
	my $max = -100000;
	foreach my $i (@l){
		$max = $i if($i > $max) ;
	}
	return $max ;
}


sub util_getCodonTable{
my %DNA_code = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V',);

return \%DNA_code ;

}

sub util_getComplimentaryString{
   my ($fe) = @_ ; 

   my @l = split "", $fe ;
   
   my @final;
   foreach my $l (@l){
	   if($l eq "A"){
		   push @final, "T";
	   }
	   elsif($l eq "T"){
		   push @final, "A";
	   }
	   elsif($l eq "G"){
		   push @final, "C";
	   }
	   elsif($l eq "C"){
		   push @final, "G";
	   }
	   elsif($l eq "R"){
	   	   my $r = rand();
		   if($r > 0.5){
		       push @final, "C";
		   }
		   else{
		       push @final, "T";
		   }
	   }
   }
   
   
   my $rev = "";
   foreach my $l (@final){
	   $rev = $rev . $l ;
   }
   my $NNN = @final ;
   
   $rev = reverse $rev ;
   return $rev ;
}
sub util_ParseBlastPW{
	my ($infile) = @_ ;
	my $ifh = util_read($infile);
	while(<$ifh>){
		if(/Expect/){
			my ($val) = (/Expect = (.*),/);
			if(!defined $val){
			($val) = (/Expect = (.*)/);
			}
			return $val ;
		}
	}
	close($ifh);
	return -1 ; 
}

sub util_ReadCodonBiasFile{
   my ($codonfile) = @_ ;
   my $ifh = util_read($codonfile);
   
   my $sort = {};
   while(<$ifh>){
   	 next if(/^\s*$/);
   	 my ($nm,@l) = split ; 
   	$sort->{$nm} = $_ ;
   }
   my $info = {};
   my $mapValues = {};
   foreach my $k (sort keys %{$sort}){
   	$_ = $sort->{$k};
       my $at = 0 ;
       my $gc = 0 ;
       my $py = 0 ;
       my $pu = 0 ;
   	next if(/^\s*$/);
   	my ($nm,@l) = split ; 
   	my $sum = 0 ;
   
   	$info->{$nm} = {};
   	$mapValues->{$nm} = {};
   
   	my $start = 0 ;
   	foreach my $i (@l){
   	 	$i =~ s/=/ /;
   		my ($id,$junk) = split " ",$i ;
   	    $mapValues->{$nm}->{$id} = $junk ;
   
   
   		my $val = 10* $junk ;
   		my $end = $val+$start ;
   		#print "$id $junk $val $end \n";
   		foreach my $idx ($start..$end){
   			$info->{$nm}->{$idx} = $id ;
   		}
   		$start = $end;
   	}
   }
   return ($info,$mapValues) ;
}

sub util_ExtractOneIdxFromFile{
	my ($infile,$idx,$varname,$fix) = @_ ;
   my @list ;
   my $ifh = util_read($infile);
   while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ;
	 my $N = @l -1 ;
	 next if($idx > $N);

	 my $NNN = $l[$idx];

	 $NNN =~ s/\(//;
	 $NNN =~  s/\)//;

	 if(defined $fix){
	 	$NNN = $NNN/$fix;
	 }
	 push @list, $NNN ;
   }
   close($ifh);


   my ($mean,$sd) = util_GetMeanSD(\@list);
   my $join = join " ,", @list ;
   my $assign = " $varname <- c ($join) \n";
   return (\@list,$mean,$sd,$assign) ;
}
