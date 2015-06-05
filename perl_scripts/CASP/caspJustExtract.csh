#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : "
  exit 
endif 


set PWD = ` pwd`
set listref = $1

foreach i ( ` cat $listref` )
  mkdir -p PR.$i
  \cp -f $i PR.$i 
  cd PR.$i 
  	tar -xvzf $i
	mkdir BACK
	\cp -f * BACK
	\rm -f *
	cd T*
	mkdir APBSDIR
	setenv APBSDIR APBSDIR
	setenv PDBDIR $cwd 

	\rm *pdb 
	ls  | grep -v APBS | grep -v list > ! list
    foreach i ( ` cat list` )
		mv $i $i.pdb 
	end
	getPDBModel1ChainAlist.csh list
  cd $PWD

end

