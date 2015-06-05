#!/bin/csh -f

if($#argv != 5  ) then 
    echo "Wrong args, required 5 " ; exit 
endif 

set trs=$1
set percentlength=$2
set percentmatched=$3
set percentidentity=$4
set expect=$5

## set up the DB
setenv BLASTDB /home/sandeepc/DATA/rafael/walnut/DB
set    TRSDB=DB.list.transcriptome.merged

# if this is normalized, keep the cutoff low
set CUTOFFTRANSCRIPTOME=0
set RAWCOUNT=/home/sandeepc/DATA/rafael/walnut/raw_bwa_counts_by_trs.txt.normalized.merged
#set RAWCOUNT=/home/sandeepc/DATA/rafael/walnut/raw_bwa_counts_by_trs.txt.normalized.merged.overlap.final

mkdir -p BLASTOUT 


echo "Info:blasting to  $TRSDB "
blastp -db $TRSDB -query $FASTADIR/$trs.ALL.1.fasta -out BLASTOUT/$trs.blast

newfile.csh blast.prot.$trs.out.reflection

echo "Info:finding related transcripts"
parseBlastLatest.pl -outf blast.prot.$trs.out.reflection -inf BLASTOUT/$trs.blast -trs $trs -percentlength $percentlength -percentmatched $percentmatched -percentidentity $percentidentity -expe $expect 

