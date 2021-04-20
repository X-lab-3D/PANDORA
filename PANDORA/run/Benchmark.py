import sys
sys.path.append('./')
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database
import os
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.multiprocessing import freeze_support
import time

# db = Database.Database().load('test_db')
# db.construct_database()

# import dill
# db = dill.load(open("test_db", 'rb'))


# ------ Multithreading test ------


# num_cores = 2 # multiprocessing.cpu_count()

def run_pandora(tar_temp):
    target = tar_temp[0]
    template = tar_temp[1]
    filename = '130421_20mod5_sd01_slow_benchmark_II.csv'
    nr_models = 5
    if not os.path.exists(filename):
        with open(filename, 'w') as f:
            f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s\n' % (
            ','.join(['molpdf_' + str(i + 1) for i in range(nr_models)]),
            ','.join(['L-RMSD_' + str(i + 1) for i in range(nr_models)]),#  ))#,
            ','.join(['core_L-RMSD_' + str(i + 1) for i in range(nr_models)])))

    try:
        mod = Pandora.Pandora(target, template = template)
        mod.model(n_models=nr_models, stdev=0.1, benchmark=True)

        moldpdf = ','.join([str(round(float(i.moldpf), 4)) for i in mod.results])
        lmrsd = ','.join([str(round(i.lrmsd, 4)) for i in mod.results])
        core_lmrsd = ','.join([str(round(i.core_lrmsd, 4)) for i in mod.results])

        with open(filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide,
            ';'.join(mod.template.allele_type), moldpdf, lmrsd, core_lmrsd))
    except:
        print('Something went wrong')

def run_multiprocessing(func, i, num_cores):
    with Pool(processes=num_cores) as pool:
        return pool.map(func, i)

def bench_MHCII():
    t0 = time.time()
    # print(t0)
    db = Database.Database().load('13_04_21_Pandora_db')
    num_cores = 10

    list_of_targets_templates = []
    for k in db.MHCII_data:
        try:
            tar = PMHC.Target(db.MHCII_data[k].id, db.MHCII_data[k].allele_type, db.MHCII_data[k].peptide,
                                               M_chain_seq=db.MHCII_data[k].M_chain_seq,
                                               N_chain_seq=db.MHCII_data[k].N_chain_seq,
                                               MHC_class=db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
            mod = Pandora.Pandora(tar, db)
            mod.find_template(benchmark=True)
            list_of_targets_templates.append((tar, mod.template))
        except:
            pass

        # run_pandora(list_of_targets_templates[0])

    run_multiprocessing(run_pandora, list_of_targets_templates, num_cores)

    print(time.time()-t0)

def bench_MHCI():
    t0 = time.time()
    # print(t0)
    db = Database.Database().load('13_04_21_Pandora_db')
    num_cores = 10

    list_of_targets_templates = []
    for k in db.MHCI_data:
        try:

            tar = PMHC.Target(db.MHCI_data[k].id, db.MHCI_data[k].allele_type, db.MHCI_data[k].peptide,
                                                M_chain_seq=db.MHCI_data[k].M_chain_seq,
                                                anchors=db.MHCI_data[k].anchors)

            mod = Pandora.Pandora(tar, db)
            mod.find_template(benchmark=True)
            list_of_targets_templates.append((tar, mod.template))
        except:
            pass

        # run_pandora(list_of_targets_templates[0])

    run_multiprocessing(run_pandora, list_of_targets_templates, num_cores)

    print(time.time()-t0)

if __name__ == "__main__":
    # Parallel(n_jobs=num_cores)(delayed(run_pandora)(target) for target in list_of_targets[:2])
    freeze_support()  # required to use multiprocessing
    bench_MHCII()
    # bench_MHCI()








