#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  exit 
endif 


/home/sandeepc/DATA/externaltools/probis-2.4.2/probis -f1 $PDBDIR/$1.pdb -c1 A -f2 $PDBDIR/$2.pdb -c2 A -compare -super > & ! /dev/null 
#/home/sandeepc/DATA/externaltools/probis-2.4.2/probis -f1 $PDBDIR/$1.pdb -c1 A -f2 $PDBDIR/$2.pdb -c2 A -compare -super 

parseProbis.pl -p1 $1 -p2 $2 -out $3

