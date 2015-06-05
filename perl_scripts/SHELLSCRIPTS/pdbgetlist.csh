#!/bin/csh -f


set list = $PWD/$1 

cd $PDBDIR 
newfile.csh $list.new
foreach i (`cat $list`)
   if(! -e ${i}A.pdb) then 
   	  if(! -e $i.pdb) then
         wget    http://www.rcsb.org/pdb/files/$i.pdb 
	  endif 
      splitpdbIntChains.pl -p $i  -outf $list.new
   endif 
   #echo ls ${i}?.pdb
   ls ${i}?.pdb >> $list.new 
end 


$SRC/SHELLSCRIPTS/removePDBfromname.csh $list.new
UNIQ $list.new

getPDBModel1ChainAlist.csh $list.new


#foreach i (`cat list.new`)
   #\mv -f $i.pdb $PDBDIR
#end 
#apbs.csh list.new
 

