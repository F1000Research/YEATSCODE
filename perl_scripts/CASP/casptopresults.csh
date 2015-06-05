#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : "
  exit 
endif 


set PWD = ` pwd`
set listref = $1

foreach i ( ` cat $listref` )
  cd PR.$i 
	set target=`ls | grep "^T"`

	cd T*
	pwd 

	setenv APBSDIR APBSDIR
	setenv PDBDIR $cwd 

	caspresults.csh
	caspscore.pl -in xxxx -lis results -model 1 -target $target

	\cp -f $target.out ~/Desktop

  cd $PWD

end

