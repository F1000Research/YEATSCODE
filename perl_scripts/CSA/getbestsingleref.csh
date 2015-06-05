#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set PDB = $1
set list = $2
set sortdir = $3
newfile.csh output 

foreach i (`cat $list`)
    if(-e $i/$i.$PDB.pdb.out ) then 
        echo -n "$i " >> output
        head -2 $i/$i.$PDB.pdb.out | grep RES  >> output
    endif 
end 

sort.pl -in output -out $sortdir/$PDB.output.sorted -idx 5

