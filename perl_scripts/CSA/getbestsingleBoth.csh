#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set cutoff = $1

ln -s ../ANNOTATE .

set sortdir = SORTDIR

mkdir -p $sortdir

$SRC/CSA/getbestsinglequeryList.csh list.ref list.query $sortdir
set list=list.query
mkdir -p SINGLEQUERY
newfile.csh tmplist 
foreach i (`cat $list`)
   extractindexfromfile.pl -in $sortdir/$i.output.sorted -idx 0 -out A > & ! /dev/null
   extractindexfromfile.pl -in $sortdir/$i.output.sorted -idx 5 -out B > & ! /dev/null
   concatfilePerline.pl A B -out SINGLEQUERY/$i.list.outann 
   echo $i SINGLEQUERY/$i.list.outann >> tmplist 
end 
annotate.pl -mapping -in ~/pdb_seqres.txt -cutoffscore $cutoff -list tmplist -anndis 5 -scores




$SRC/CSA/getbestsinglerefList.csh list.ref list.query $sortdir
set list=list.ref
mkdir -p SINGLEREF
newfile.csh tmplist 
foreach i (`cat $list`)
   extractindexfromfile.pl -in $sortdir/$i.output.sorted -idx 0 -out A > & ! /dev/null
   extractindexfromfile.pl -in $sortdir/$i.output.sorted -idx 5 -out B > & ! /dev/null
   concatfilePerline.pl A B -out SINGLEREF/$i.list.outann 
   echo $i SINGLEREF/$i.list.outann >> tmplist 
end 
annotate.pl -mapping -in ~/pdb_seqres.txt -cutoffscore $cutoff -list tmplist -anndis 5 -scores

unlink tmplist 
