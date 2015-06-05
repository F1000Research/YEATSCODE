#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : "
  exit 
endif 


set PWD = ` pwd`
set listref = $1

foreach i ( ` cat $listref` )
  cd PR.$i 
	cd T*
	mkdir -p APBSDIR
	setenv APBSDIR APBSDIR
	setenv PDBDIR $cwd 

	apbs.csh list 
   
    newfile.csh NEWLIST
    foreach j ( ` cat list` )
		if(! -e APBSDIR/$j/pot1.dx.atompot) then 
			#\rm -rf APBSDIR/$j/
			#echo $j >> NEWLIST
		endif 
	end
	apbs.csh NEWLIST

  cd $PWD

end

