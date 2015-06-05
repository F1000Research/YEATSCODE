#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set val = $2 
set lines = $3 

\rm x*
split -l $lines $list

wc -l $list
wc -l x*

foreach i (`ls x*`)
    
	echo SPLITTING $i
	checkIdentity.pl -out $i.LLL -simi $val -list $i -needleout oooooooooo -arg ~/needle.arg > & ! /dev/null
end 

cat *.LLL > ! LLL
echo Now doing for all...numebr of lines is
wc -l LLL
checkIdentity.pl -out LLL.out -simi 80 -list LLL -needleout oooooooooo -arg ~/needle.arg > & ! /dev/null

 

