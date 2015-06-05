#!/bin/csh  

if($#argv != 1  ) then 
  echo "Usage : list "
  echo "You need to set ,  list"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


set PWD=`pwd`
set list = $1 
set listuniq = $1.new.uniq

## this takes the 4 letter pdb
pdbgetlist.csh   $list 

processOnePDBforallPropertieslist.csh $listuniq

echo Doing APBS fasta ...
apbs.csh         $listuniq

echo Doing fpocket fasta ...
fpocketlist.csh     $PWD/$listuniq

