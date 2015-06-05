#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : "
  exit 
endif 

set PWD = ` pwd`
set listref = $PWD/$1

set split=$listref.splitted
newfile.csh $split
foreach ref ( ` cat $listref` )
   $SRC/MISC/splitpdbIntChains.pl -p $ref -out $split 
end



setenv PDBDIR $PWD

orderHetAtlast.pl -li $split  
getPDBModel1ChainAlist.csh $split 
