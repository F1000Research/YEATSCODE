#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : protein name "
  exit 
endif 

if(-e  $PDBDIR/$1.pdb) then 
   echo "Copying from pdb"
   cp -f $PDBDIR/$1.pdb .
else
   wget 	http://www.rcsb.org/pdb/files/$1.pdb
endif 

getPDBModel1ChainA.pl -protein $1 -model



\rm $1.pdb
\rm $1*-m*.pdb $1*-c*.pdb 

mv -f A.pdb $1.pdb 
