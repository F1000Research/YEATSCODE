#!/bin/csh 


#### DB/all.fasta is the initial files - make sure you have the name correct for the TRS's, 
#### First get a list of trs - remove exact ones.
#### for example, no semicolons, etc

set initialfile = DB/all.fasta
if(! -e $initialfile) then
	echo "Error: Need $initialfile file";
	exit 
endif


# Make BLAST DB
if(! -e DB/all.fasta.nsd) then
    echo ---------- Made BLAST DB for $initialfile
    makeblastdbN.csh DB/all.fasta
endif

mkdir -p FASTADIR_NT/
mkdir -p ORF/
mkdir -p FASTADIR_ORF/


# 1) Get list of TRSs, removing exact ones
setenv FASTADIR $cwd/FASTADIR_NT/
if(! -e list.trs) then
	newfile.csh mapname2length
    pruneSameSequenceFromMadeFasta.pl -outf list.trs -inf DB/all.fasta -writedata 1 -length mapname2length
	replacestring.pl -wit "" -whic ";" -outf kkkk -inf mapname2length -same
    extractindexfromfile.pl -in list.trs -out ttt -idx 0 -sep ""
    \mv -f ttt list.trs

	## list.trs.mapping.sort has the exact same ones
	sort.pl -in list.trs.mapping -idx 1 > & ! /dev/null & 
endif 

set totalnumber=`wc -l list.trs`
echo "----------There are $totalnumber of transcripts"

set finalcommands=finalcommands.csh
newfile.csh $finalcommands


mkdir -p FASTADIR_NT/
mkdir -p ORF/
mkdir -p ORFTMP/
mkdir -p FASTADIR_ORF/



### Get orfs, parse them in 3 ORF per TRS

set orfsinonefile=DB/orfsinonefile.fasta
newfile.csh $orfsinonefile
foreach i (`cat list.trs`)
	if(! -e ORF/$i.orf) then 
	   getorf FASTADIR_NT/$i.ALL.1.fasta ORFTMP/$i.orf
	   replacestring.pl -wit ".ORF_" -whic "_" -inf  ORFTMP/$i.orf -outf ORF/$i.orf 
	   findrepeatedorfs.pl -trs $i -orfdir ORF/ -write 4

	   ## list.err.repeat.sort has TRS's with exactly repeating ORF's
       sort.pl -in list.err.repeat -idx 1

	   cat ORF/$i* >> $orfsinonefile
	endif 

end

## Make BLAST DB for all ORFs
if(! -e "$orfsinonefile.nsd") then
    echo ---------- Made BLAST DB for $orfsinonefile
    makeblastdbP.csh $orfsinonefile 
endif

ls FASTADIR_ORF/ > ! list.trs.orfs
replacestring.pl -whic ".ALL.1.fasta" -with "" -outf kkkk -inf list.trs.orfs -same


## BLAST each TRS to complete set
## Generate command files for BLAST of ALL TRSs to itself
newfile.csh runblasttoitself.csh 
echo "mkdir -p BLASTOUT_NT_TOITSELF" >> runblasttoitself.csh 
foreach i (`cat list.trs`)
	if(! -e BLASTOUT_NT_TOITSELF/$i.blast) then 
       echo "runoneblastn.csh DB/all.fasta FASTADIR_NT/$i.ALL.1.fasta  BLASTOUT_NT_TOITSELF//$i.blast >  & ! /dev/null & " >> runblasttoitself.csh 
	 endif
end
scheduleprocessInsertingsleep.pl -inf runblasttoitself.csh  -int 500 -sleep 10
echo source Sch.runblasttoitself.csh >> $finalcommands

    

newfile.csh mapRedundancyNT
foreach i (`cat list.trs`)
      $SRC/GNM/parseBlastLatest.pl -outf mapRedundancyNT -inf BLASTOUT_NT_TOITSELF//$i.blast -trs $i -verb 0 -percentidentity 95 -percentmatched 90 -percentlength 20
end 
groupBasedonCutoff.pl -outf mapRedundancyNT.groups -inf mapRedundancyNT -cutoff 0.001 -dir 0
sort.pl -idx 3 -in mapRedundancyNT.groups -rev
extractindexfromfile.pl -idx 5 -in mapRedundancyNT.groups.sort
# mapRedundancyNT.groups.sort.5 has the list of groups sorted with largest




echo "#You may want to run a minimizer first before you run BLAST on all the ORFs"  >> $finalcommands

## Generate command files for BLAST of ALL ORFS
newfile.csh runallblast.csh 
echo "mkdir -p BLASTOUT" >> runallblast.csh 
foreach i (`cat list.trs.orfs`)
    echo "runoneblastp.csh nr FASTADIR_ORF/$i.ALL.1.fasta  BLASTOUT/$i.blast >  & ! /dev/null & " >> runallblast.csh 
end
scheduleprocessInsertingsleep.pl -inf runallblast.csh  -int 0 -sleep 30
echo source Sch.runallblast.csh >> $finalcommands









