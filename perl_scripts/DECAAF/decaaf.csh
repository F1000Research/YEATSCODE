#!/bin/csh

if($#argv != 3  ) then 
  echo "Usage : "
  exit 
endif 

set PWD = ` pwd`
set top = $PWD/$1
set listnumber = $PWD/$2
set listpdbs = $3

set pdb = ` cat $top` 
echo "Working for pdb $pdb"

foreach ref ( ` cat $listnumber` )

	echo "running for $ref"
	mkdir -p $ref.run
	cd $ref.run
    createCLASPinput.csh $pdb ../$ref.in 4 4 
	ln -s ANNOTATE.4 ANNOTATE 
	\cp -f ../list* . 
	cd -
end

foreach ref ( ` cat $listnumber` )
	cd $ref.run/
	$SRC/CSA/runRefExtractEasilyNamed.csh list.2OQIA list.1PTDA 0
	cd -

end
