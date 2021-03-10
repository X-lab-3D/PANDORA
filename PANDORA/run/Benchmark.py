from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.Database()
db.construct_database()

import dill
# db = dill.load(open("Pandora_MHCI_and_MHCII_data", 'rb'))


##------------------------------ MHCI -----------------------------------

with open('benchmark_I.csv', 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')

for k in db.MHCI_data:
    try:
        target = PMHC.Target(db.MHCI_data[k].PDB_id, db.MHCI_data[k].allele, db.MHCI_data[k].peptide, M_chain_seq= db.MHCI_data[k].M_chain_seq, MHC_class= db.MHCI_data[k].MHC_class, anchors=db.MHCI_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(benchmark=True)

        lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
        core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)

        with open('benchmark_I.csv', 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.PDB_id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.PDB_id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))

    except:
        print('Something went wrong')


##------------------------------ MHCII -----------------------------------


with open('benchmark_II.csv', 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')

for k in db.MHCII_data:
    try:
        target = PMHC.Target(db.MHCII_data[k].PDB_id, db.MHCII_data[k].allele, db.MHCII_data[k].peptide, M_chain_seq= db.MHCII_data[k].M_chain_seq, N_chain_seq= db.MHCII_data[k].N_chain_seq, MHC_class= db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(benchmark=True)

        lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
        core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)


        with open('benchmark_II.csv', 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.PDB_id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.PDB_id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))

    except:
        print('Something went wrong')
#