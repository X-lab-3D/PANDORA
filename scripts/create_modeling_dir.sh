#/usr/bin/env bash
#!/usr/bin/env bash
# Farzaneh Meimandi Parizi
# 05-June-2019

# Creating the required directries, downloading PDB files and also copiying related files.

set -e

if [ $# -ne 4 ]
  then
      echo "Usage: <Parent_directory_for_all_runs> <Modelling_directory_name> <Template_structure_PDB_id>  <Query_structure_PDB_id>"
      exit 1
fi


# Varibales
central_dir=$1 # Directory in which all the modelling file folders will be created.
model_dir=$2 # Folder inside the central directory that can contains multiple run_names of the same modelling. Example: "3kpp_3kpm"
template_structure_id=$3 # 4 character PDB ID. Example: "3kpm"
model_structure_id=$4 # 4 character PDB ID. Example: "3kpm"

# Move to the central directory
mkdir -p $central_dir; cd $central_dir
# Make the directory that can contain several individual runs of the same modelling
mkdir -p $model_dir; cd $model_dir

/home/fmeiman/tools/PDB_related/pdb_download.sh  $template_structure_id
/home/fmeiman/tools/PDB_related/pdb_download.sh   $model_structure_id



shopt -s extglob
#rm -r !(*.pdb)

cp /home/fmeiman/tools/modeller/MHCI_modeling/MyLoop.py  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/cmd_modeller.py  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/eval_model.py  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/run_file.sh  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/restraint_tmp.txt  ./
#cp /home/fmeiman/tools/modeller/MHCI_modeling/MyModel.py  ./
#cp /home/fmeiman/tools/modeller/MHCI_modeling/cmd_modeller_brk.py ./
#cp /home/fmeiman/tools/modeller/MHCI_modeling/P2_P9.rsr  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/cmd_modeller_soap.py  ./
cp /home/fmeiman/tools/modeller/MHCI_modeling/run_file_soap.sh  ./
# Print progress message.
echo 'The environment is ready for adding restraints for modelling!'
echo `pwd`
