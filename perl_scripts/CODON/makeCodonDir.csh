#!/bin/csh  -f

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set nm = $1
set file = $2

mkdir -p $nm/DB
\cp -f $file* $nm/DB

cd $nm/DB

if( -e $file.gz ) then 
echo unzipping
gunzip $file.gz > & ! /dev/null 
endif 

\rm all.fasta
ln -s $file all.fasta

lc.pl -out ooo -in $file -uc 
\mv -f ooo $file

cd ../
$SRC/GNM/yeatsTOP_AA.csh > & ! log & 
