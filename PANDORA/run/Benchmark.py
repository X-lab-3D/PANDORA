from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
import time
import dill


db = dill.load(open("Pandora_MHCI_and_MHCII_data", 'rb'))


##------------------------------ MHCI -----------------------------------

with open('benchmark_I.csv', 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')

for k in db.MHCI_data:
    try:
        t0 = time.time()
        print('Modelling %s' %db.MHCI_data[k].PDB_id)
        target = PMHC.Target(db.MHCI_data[k].PDB_id, db.MHCI_data[k].allele, db.MHCI_data[k].peptide, chain_seq= db.MHCI_data[k].chain_seq, MHC_class= db.MHCI_data[k].MHC_class, anchors=db.MHCI_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(benchmark=True)

        lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
        core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)

        print('Mean L-RMSD: %s' %(lmrsd))
        print('Mean core L-RMSD: %s' % (core_lmrsd))
        print('Modelling took %s seconds\n' %(time.time() - t0))

        with open('benchmark_I.csv', 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.PDB_id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.PDB_id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))

    except:
        print('Something went wrong')


##------------------------------ MHCII -----------------------------------


with open('benchmark_II.csv', 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')

for k in db.MHCII_data:
    try:
        t0 = time.time()
        print('Modelling %s' %db.MHCII_data[k].PDB_id)
        target = PMHC.Target(db.MHCII_data[k].PDB_id, db.MHCII_data[k].allele, db.MHCII_data[k].peptide, chain_seq= db.MHCII_data[k].chain_seq, MHC_class= db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(benchmark=True)

        lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
        core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)

        print('Mean L-RMSD: %s' %(lmrsd))
        print('Mean core L-RMSD: %s' % (core_lmrsd))
        print('Modelling took %s seconds\n' %(time.time() - t0))

        with open('benchmark_II.csv', 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.PDB_id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.PDB_id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))

    except:
        print('Something went wrong')
#