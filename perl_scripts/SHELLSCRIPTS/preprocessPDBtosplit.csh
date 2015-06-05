#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD=`pwd`

mkdir -p split 
cd split 

set list = $PWD/$1
set F = list.split    

touch $F
foreach i (`cat $list`)
	 if(! -e $i.split) then 
		 touch $i.split 
	     splitpdbIntChains.pl -outf $F -pr $i
	 endif 
end 

setenv PDBDIR $cwd

getPDBModel1ChainAlist.csh $F

foreach i (`cat $F`)
    if(! -e $i.ALL.1.fasta) then
        writeFasta.pl -p $i 
	endif 
end 

cat *fasta > ! F
processStringFromFasta.pl -outf oo -in F


