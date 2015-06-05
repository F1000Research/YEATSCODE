#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set protein = $1
set list = $2
newfile.csh F 
foreach i (`cat $list`)
	echo -n "$i  " >> F 
    head -1 $protein.$i.pdb.out  | grep -v does >> F 
	echo "" >> F
end 

sort.pl -in F -idx 6 -out FF
extractindexfromfile.pl -idx 0 -in FF -out A 
extractindexfromfile.pl -idx 6 -in FF -out B 

concatfilePerline.pl A B -out $protein.cumu.scores
removeChainIdKeepingsameorder.csh A A1
ANNNEW A1
