#!/bin/tcsh -f

# This script is to be used with ANA_RMSD-Split.csh

# EzgiKaraca, 20092011, 6:19PM

# Usage: l-rmsd-calc.csh ref file.list ('ref' is name of reference PDB file and the .lzone file).

set WDIR = $PWD
set refe = $WDIR/$argv[1].pdb
set lzone = $WDIR/$argv[1].lzone
set atoms = 'CA,CB,C,O,N'

#foreach i ($argv[2])
#  pdb_segxchain $i >$i:r.tmp1

echo "" > l-rmsd_xray.disp
foreach i (` cat $argv[2] `)
  #/home/daans/pdb_related/pdb-tools/pdb_segxchain.py $i >$i:r.tmp1
  #pdb_segxchain $i >$i:r.tmp1
  echo $i
  echo $i >>l-rmsd_xray.disp
  /home/dariom/ProFitV3.3/bin/profit <<_Eod_ |grep RMS |tail -2 >>l-rmsd_xray.disp
    refe $refe
    mobi $i
    ratoms $atoms
    atom $atoms
    `cat $lzone`
    quit
_Eod_

#\rm $i:r.tmp1
end
awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' l-rmsd_xray.disp |grep pdb |awk '{print $1,$2}' > BB-CB-l-RMSD.dat
#sort -nk2 BB-CB-l-RMSD.dat > BB-CB-l-RMSD_sorted.dat
rm -rf l-rmsd_xray.disp
