from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.Database().load('test_db')
# db.construct_database()

# import dill
# db = dill.load(open("test_db", 'rb'))


##------------------------------ MHCI -----------------------------------
nr_models = 20
filenameI = '160321_benchmark_I.csv'

with open(filenameI, 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s\n' %(','.join(['molpdf_' + str(i+1) for i in range(nr_models)]), ','.join(['L-RMSD_' + str(i+1) for i in range(nr_models)]), ','.join(['core_L-RMSD_' + str(i+1) for i in range(nr_models)])))

for k in db.MHCI_data:
    # pass
    try:
        target = PMHC.Target(db.MHCI_data[k].id, db.MHCI_data[k].allele_type, db.MHCI_data[k].peptide, M_chain_seq= db.MHCI_data[k].M_chain_seq, MHC_class= db.MHCI_data[k].MHC_class, anchors=db.MHCI_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(benchmark=True)
        # mod.results
        moldpdf = ','.join([str(round(float(i.moldpf), 4)) for i in mod.results])
        lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
        core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)

        with open(filenameI, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide, ';'.join(mod.template.allele_type),moldpdf, lmrsd, core_lmrsd))

    except:
        print('Something went wrong')


##------------------------------ MHCII -----------------------------------

nr_models = 20
filename = '160321_benchmark_II.csv'

with open(filename, 'w') as f:
    f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s\n' %(','.join(['molpdf_' + str(i+1) for i in range(nr_models)]), ','.join(['L-RMSD_' + str(i+1) for i in range(nr_models)]), ','.join(['core_L-RMSD_' + str(i+1) for i in range(nr_models)])))

for k in db.MHCII_data:
    # pass
    try:
        target = PMHC.Target(db.MHCII_data[k].id, db.MHCII_data[k].allele_type, db.MHCII_data[k].peptide, M_chain_seq= db.MHCII_data[k].M_chain_seq, N_chain_seq= db.MHCII_data[k].N_chain_seq, MHC_class= db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
        mod = Pandora.Pandora(target, db)
        mod.model(n_models=nr_models, stdev=0.2, benchmark=True)

        moldpdf = ','.join([str(round(float(i.moldpf),4)) for i in mod.results])
        lmrsd = ','.join([str(round(i.lrmsd,4)) for i in mod.results])
        core_lmrsd = ','.join([str(round(i.core_lrmsd, 4)) for i in mod.results])

        with open(filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide, ';'.join(mod.template.allele_type),moldpdf, lmrsd, core_lmrsd))

    except:
        print('Something went wrong')
#



