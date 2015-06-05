#!/bin/csh 

if($#argv != 5  ) then 
  echo "Usage : db list FASTADIR BLASTOUT"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set db = $1 
set list = $2 
set FASTADIR=$3
set BLASTOUT=$4
set BLASTPROG=$5

mkdir -p BLASTOUT

foreach i (`cat $list`)
	if(! -e $FASTADIR/$i.ALL.1.fasta) then 
	   echo "Could not find  $FASTADIR/$i.ALL.1.fasta "
	   exit 
	endif 
	$BLASTPROG -db  $db -query $FASTADIR/$i.ALL.1.fasta -out $BLASTOUT/$i.blast
end 


