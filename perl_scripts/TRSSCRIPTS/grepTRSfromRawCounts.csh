#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set rawcount = $2

newfile.csh $3
#echo " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB " >> trs_counts
foreach i (`cat $list`)
     grep -i $i $rawcount >>     $3
end 

