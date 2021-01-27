# this script maps and renumbers model.pdb according to ref.pdb
# Here is an example.

ref=$1
model=$2

egrep '^(ATOM|HETATM|END)' $model > clean_$model
for chnID in {M,P};do
../../../tools/pdb-pdbalign $ref $chnID clean_$model $chnID > common.pdb

# delete the warning line in common.pdb
sed -i '/Warning/d' common.pdb

# delete the lines with residue name of 'X'
awk '$1 == "ATOM" && substr($0,22,1) != "X"' common.pdb  >  common_$chnID.pdb

done

cat common_M.pdb common_P.pdb > common.pdb
sed '$ a END' common.pdb
