#!/usr/bin/env bash
# Farzaneh Meimandi Parizi
# 10-June-2019

# This script prepares the files by adding restraints to run Modeller
set -e

if [ $# -lt 4 ]
  then
      echo "Usage: <Template_ID>  <Target_ID> <Model_directory> <_interface_cutoff_distance>"
      exit
fi

tmpl_ID=$1 # Template file
int_cut=$4 # Interface cutoff value
query_ID=$2 # Name of the target PDB ID
model_path=${3:-`pwd`}; cd $model_path # Move to the model directory if script is called automatically.


# Other variables.
model_FN='template_'$tmpl_ID'_modelled_to_'$query_ID # File name of the modelled template
run_name=$tmpl_ID'_'$query_ID

query_pept_chain_ID=(`bash /home/fmeiman/tools/find_pept_chain_id.sh $query_ID'.pdb'`) # Finds the chain ID of the peptide in model file

query_mhc_chain_ID=(`bash /home/fmeiman/tools/find_mhc_chain_id.sh $query_ID'.pdb'`) # Finds the chain ID of the MHC in query file

tmpl_pept_chain_ID=(`bash /home/fmeiman/tools/find_pept_chain_id.sh $tmpl_ID'.pdb'`) # Finds the chain ID of the peptide in model file

tmpl_mhc_chain_ID=(`bash /home/fmeiman/tools/find_mhc_chain_id.sh $tmpl_ID'.pdb'`) # Finds the chain ID of the MHC in query file

echo 'Chain ID of peptide in model is: '$query_pept_chain_ID
echo 'Chain ID of mhc in model is: '$query_mhc_chain_ID
echo 'Chain ID of peptide in template is: '$tmpl_pept_chain_ID
echo 'Chain ID of mhc in template is: '$tmpl_mhc_chain_ID

#########################################################
##### PART 1: Preparation of MODELLER python scripts#####
#########################################################
query_new_FN=$query_ID'_'$query_mhc_chain_ID$query_pept_chain_ID

/home/fmeiman/tools/PDB_related/pdb_selchain.py -$query_mhc_chain_ID $query_ID'.pdb' > $query_ID'_mhc.pdb'
egrep '^ATOM' $query_ID'_mhc.pdb' > $query_new_FN'.pdb'

/home/fmeiman/tools/PDB_related/pdb_selchain.py -$query_pept_chain_ID $query_ID'.pdb' > $query_ID'_pept.pdb'
egrep '^ATOM' $query_ID'_pept.pdb' >> $query_new_FN'.pdb'

tmpl_new_FN=$tmpl_ID'_'$tmpl_mhc_chain_ID$tmpl_pept_chain_ID
tmpl_new_FN_rechained=$tmpl_ID'_'$tmpl_mhc_chain_ID$query_pept_chain_ID

/home/fmeiman/tools/PDB_related/pdb_selchain.py -$tmpl_mhc_chain_ID $tmpl_ID'.pdb' > $tmpl_ID'_mhc.pdb'
egrep '^ATOM' $tmpl_ID'_mhc.pdb' > $tmpl_new_FN'.pdb'

/home/fmeiman/tools/PDB_related/pdb_selchain.py -$tmpl_pept_chain_ID $tmpl_ID'.pdb' > $tmpl_ID'_pept.pdb'
egrep '^ATOM' $tmpl_ID'_pept.pdb' >> $tmpl_new_FN'.pdb'

#python /home/lixue/tools/PDB_related/pdb_rechain.py -$tmpl_pept_chain_ID -$query_pept_chain_ID $tmpl_new_FN'.pdb' >> $tmpl_new_FN_rechained'.pdb'

# Retain only one conformation of the atom in the new file.
python /home/fmeiman/tools/PDB_related/pdb_selalt.py -A $tmpl_new_FN'.pdb' > query_single_conformers_no_reset.pdb

python /home/fmeiman/tools/PDB_related/pdb_reres_chain.py -1 -$tmpl_pept_chain_ID query_single_conformers_no_reset.pdb > query_single_conformers.pdb # The residue numbers of the peptide are now set to 1. NOTE: Also HOH ATOMs numbers have been changed.

## Keep water lines
#egrep '^HETATM' template_single_conformers.pdb > template_hetatm.pdb
#/home/lixue/tools/pdb_keepWat.py template_hetatm.pdb > template_hetatm.pdb_new
#mv template_hetatm.pdb_new template_hetatm.pdb
#
#egrep '^HETATM' query_single_conformers.pdb > query_hetatm.pdb
#/home/lixue/tools/pdb_keepWat.py query_hetatm.pdb > query_hetatm.pdb_new
#mv query_hetatm.pdb_new query_hetatm.pdb
#
## Add side chains
#grep ^ATOM template_single_conformers.pdb > t_pdb_ready_for_refinement.pdb; echo 'END' >> t_pdb_ready_for_refinement.pdb
#/home/lixue/tools/run_cns_generate.sh `pwd` t_pdb_ready_for_refinement.pdb $tmpl_FN
#
## Add side chains
#grep ^ATOM query_single_conformers.pdb > q_pdb_ready_for_refinement.pdb; echo 'END' >> q_pdb_ready_for_refinement.pdb
#/home/lixue/tools/run_cns_generate.sh `pwd` q_pdb_ready_for_refinement.pdb $query_FN
#
##-- put hetatm section back to the mutated PDB file.
#cat $tmpl_FN hetatm.pdb > $tmpl_FN.tmp
#mv $tmpl_FN.tmp $tmpl_FN
#
#cat $query_FN hetatm.pdb > $query_FN.tmp
#mv $query_FN.tmp $query_FN


############################################################
##### PART 2: Interface tables and distance restraints #####
############################################################

# Create distance restraints files for P2 (C alph and C beta) and for Psh9 #### RESET NUMBERS OF PEPTIDE
/home/lixue/tools/PDB_related/contact-chainID_allAtoms query_single_conformers.pdb $int_cut > contacts.list.tmp # Make table with distance between CA/CB and MHC molecules.

/home/lixue/tools/extract_contacts.py contacts.list.tmp $tmpl_pept_chain_ID 2 > contacts_P2.list # Only extract P2 contacts.
/home/lixue/tools/extract_contacts.py contacts.list.tmp $tmpl_pept_chain_ID 9 > contacts_P9.list # Only extract P9 contacts.
echo '> Contact files P2 and P9 are generated.'

# Convert the created distance to restraints to Haddock format.
cat contacts_P2.list contacts_P9.list > contacts_P2_P9.list # Get contacts of P2 and P9 in one file.

#/home/fmeiman/tools/contact2tbl.sh contacts_P2_P9.list $low_bound_flex $up_bound_flex > contacts_P2_P9.tbl # Insert an extra column with the flexibility for HADDOCK.
#echo '> Contact table for is generated.'

# Convert the distance restraints to a pymol file as well.
/home/fmeiman/tools/contact2pml.sh contacts_P2_P9.list > contacts_P2_P9.pml
awk '{print "show sticks, resi " $3 " and chain  " $2 }' contacts_P2_P9.list | sort -u >> contacts_P2_P9.pml
awk '{print "show sticks, resi " $7 " and chain  " $6 }' contacts_P2_P9.list | sort -u >> contacts_P2_P9.pml
echo '> Contact file for Pymol is generated.'


## Now create the contacts_P2_P9.tbl again, but this time for the haddock restraints instead of pymol. This file will overwrite the other file.
#/home/fmeiman/tools/contact2tbl_haddock.sh contacts_P2_P9.list $low_bound_flex $up_bound_flex > contacts_P2_P9.tbl

#Creating pdb file with selecting peptide and mhc chains each time for template and the query proteins


#query_new_FN=$query_ID'_'$query_mhc_chain_ID$query_pept_chain_ID
#
#pdb_selchain -$query_mhc_chain_ID $query_ID'.pdb' > $query_ID'_mhc.pdb'
#egrep '^ATOM' $query_ID'_mhc.pdb' > $query_new_FN'.pdb'
#
#pdb_selchain -$query_pept_chain_ID $query_ID'.pdb' > $query_ID'_pept.pdb'
#egrep '^ATOM' $query_ID'_pept.pdb' >> $query_new_FN'.pdb'
#
#tmpl_new_FN=$tmpl_ID'_'$tmpl_mhc_chain_ID$tmpl_pept_chain_ID
#
#pdb_selchain -$tmpl_mhc_chain_ID $tmpl_ID'.pdb' > $tmpl_ID'_mhc.pdb'
#egrep '^ATOM' $tmpl_ID'_mhc.pdb' > $tmpl_new_FN'.pdb'
#
#pdb_selchain -$tmpl_pept_chain_ID $tmpl_ID'.pdb' > $tmpl_ID'_pept.pdb'
#egrep '^ATOM' $tmpl_ID'_pept.pdb' >> $tmpl_new_FN'.pdb'


#Substituting the name of pdb files to be modelled as parameters in modeller files


#python /home/lixue/tools/PDB_related/pdb_rechain.py -$tmpl_pept_chain_ID -query_pept_chain_ID $tmplID.pdb >> $tmplID_$query_pept_chain_ID.pdb

########## -s -c "%s/P'/$query_pept_chain_ID'/g|x" MyLoop.py
ex -s -c "%s/P'/$tmpl_pept_chain_ID'/g|x" MyLoop.py
ex -s -c "%s/'A'/'$query_mhc_chain_ID'/g|x" MyLoop.py

aln_FN=$tmpl_ID'_'$query_ID'.ali'
ex -s -c "%s/tmplID_qryID.ali/$aln_FN/g|x" cmd_modeller.py

ex -s -c "%s/qryID/$query_new_FN/g|x" cmd_modeller.py
ex -s -c "%s/tmplID/$tmpl_new_FN/g|x" cmd_modeller.py

##########SOAP scoring############
ex -s -c "%s/tmplID_qryID.ali/$aln_FN/g|x" cmd_modeller_soap.py

ex -s -c "%s/qryID/$query_new_FN/g|x" cmd_modeller_soap.py
ex -s -c "%s/tmplID/$tmpl_new_FN/g|x" cmd_modeller_soap.py

#adding restraonts to the parameter files of modeller

##cat contacts_P2_P9.list | while read line
##do
##    ll=($line)
##    echo ${ll[3]}':'${ll[2]}:${ll[1]}"  "${ll[8]}:${ll[7]}:${ll[6]} 
##done > restraints.out
##rsr1="$(cat restraints.out | sed -n ''3p'')"
##rsr2=`(tail -4 restraints.out | head -1)`
##rsr1=($rsr1)
##rsr2=($rsr2)
##ex -s -c "%s/CA:2:C/${rsr1[0]}/|x"  MyLoop.py
##ex -s -c "%s/CZ:7:A/${rsr1[1]}/|x"  MyLoop.py
##ex -s -c "%s/CA:9:C/${rsr2[0]}/|x"  MyLoop.py
##ex -s -c "%s/CZ2:147:A/${rsr2[1]}/|x"  MyLoop.py

#rsr_pattern="$(cat P2_P9.rsr | sed -n ''2p'')"
#rsr2=`(tail -1 restraints.out)`

#echo $rsr_pattern
#echo $rsr1
#echo $rsr

#n=`(cat P2_P9.rsr | wc -l)`
#for ((i = 1; i <= N; i++));
#do
#   ex -s -c "%s/CE1:7:A CA:2:C/$rsr1/|x"  P2_P9.rsr
#   echo  $rsr_pattern >> P2_P9.rsr
#
#done

#echo  $rsr_pattern >> P2_P9.rsr
 
#sed '0,/CE1:7:A CA:2:C/s/CE1:7:A CA:2:C/'$rsr1'/' P2_P9.rsr > P2_P9_replaced.rsr
#sed '0,/CE1:7:A CA:2:C/s/CE1:7:A CA:2:C/'$rsr2'/' P2_P9_replaced.rsr > P2_P9_replaced2.rsr

#ex -s -c "%s/CE1:7:A CA:2:C/$rsr1/|x"  P2_P9.rsr
#echo  $rsr_pattern >> P2_P9.rsr
#ex -s -c "%s/CE1:7:A CA:2:C/$rsr2/|x"  P2_P9.rsr

cat contacts_P2_P9.list | while read line
do
    ll=($line)
    rsr1=(${ll[3]}':'${ll[2]}:${ll[1]} ${ll[8]}:${ll[7]}:${ll[6]})

    cat restraint_tmp.txt >> MyLoop.py
    ex -s -c "%s/atom1/${rsr1[0]}/|x"  MyLoop.py
    ex -s -c "%s/atom2/${rsr1[1]}/|x"  MyLoop.py
done

#printf %s 'P1>'$query_ID > $query_ID'.ali'
#/home/lixue/tools/PDB_related/PDB2SEQ.pl $query_ID'.pdb' $query_mhc_chain_ID  > $query_ID'_mhc.ali'
#tail -n +5 $query_ID'_mhc.ali' >> $query_ID'.ali'
#printf %s '*' >> $query_ID'.ali'
#/home/lixue/tools/PDB_related/PDB2SEQ.pl $query_ID'.pdb' $query_pept_chain_ID  >> $query_ID'_pept.ali'
#tail -n +5 $query_ID'_pept.ali' >> $query_ID'.ali'

#aln_FN=$tmpl_ID'_'$query_ID'.ali'

echo '>P1;query_'$query_new_FN > $aln_FN
echo 'sequence:::::::::' >> $aln_FN
/home/lixue/tools/PDB_related/PDB2SEQ.pl $query_ID'_mhc.pdb' > $query_ID'_mhc.ali'  
echo -n `(tail -n +5 $query_ID'_mhc.ali')` >>  $aln_FN
echo -n "/" >>  $aln_FN
/home/lixue/tools/PDB_related/PDB2SEQ.pl $query_ID'_pept.pdb' >> $query_ID'_pept.ali'
echo -n `(tail -n +5 $query_ID'_pept.ali')` >>  $aln_FN
echo  '*' >>  $aln_FN

echo '' >>  $aln_FN
echo '>P1;'$tmpl_new_FN >> $aln_FN
echo 'structure:'$tmpl_new_FN'.pdb:1:A:9:'$tmpl_pept_chain_ID'::::' >>  $aln_FN
/home/lixue/tools/PDB_related/PDB2SEQ.pl $tmpl_ID'_mhc.pdb'  > $tmpl_ID'_mhc.ali'
echo -n `(tail -n +5 $tmpl_ID'_mhc.ali')` >>  $aln_FN
echo -n "/" >>  $aln_FN
/home/lixue/tools/PDB_related/PDB2SEQ.pl $tmpl_ID'_pept.pdb'  >> $tmpl_ID'_pept.ali'
echo -n `(tail -n +5 $tmpl_ID'_pept.ali')` >>  $aln_FN
echo -n '*' >>  $aln_FN

#ex -s -c "%s/path/$model_path/|x" run_file.sh

echo  '> Your file is ready for modelling.'

# Remove junk files.
rm query_single_conformers_no_reset.pdb
rm query_single_conformers.pdb
#rm multi_mutations.list
# Remove contact files
rm contacts_P2.list
rm contacts_P9.list
#rm contacts_P2_P9.list
#rm single_conformers.atomResNum
rm contacts.list.tmp
rm restraint_tmp.txt
#rm P2_P9_replaced.rsr
rm $tmpl_ID'_mhc.ali'
rm $tmpl_ID'_pept.ali'
rm $query_ID'_mhc.ali'
rm $query_ID'_pept.ali'
rm $tmpl_ID'_mhc.pdb'
rm $tmpl_ID'_pept.pdb'
rm $query_ID'_mhc.pdb'
rm $query_ID'_pept.pdb'

echo $model_FN ## Warning: If you want to remove this line, please make a copy of this file. This line is important for the workflow.
