import sys
sys.path.append('./')
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database
import os
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.multiprocessing import freeze_support
import time
import PANDORA

# db = Database.Database().load('test_db')
# db.construct_database()

# import dill
# db = dill.load(open("test_db", 'rb'))


# ------ Multithreading test ------


# num_cores = 2 # multiprocessing.cpu_count()

# tar_temp = list_of_targets_templates[0]

def run_pandora(tar_temp):
    target = tar_temp[0]
    template = tar_temp[1]
    filename = '280721_default_benchmark_I.csv'
    nr_models = 20
    if not os.path.exists(PANDORA.PANDORA_path + '/../' + filename):
        with open(PANDORA.PANDORA_path + '/../' + filename, 'w') as f:
            f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s,%s\n' % (
            ','.join(['molpdf_' + str(i + 1) for i in range(nr_models+1)]),
            ','.join(['dope_' + str(i + 1) for i in range(nr_models + 1)]),
            ','.join(['L-RMSD_' + str(i + 1) for i in range(nr_models+1)]),#  ))#,
            ','.join(['core_L-RMSD_' + str(i + 1) for i in range(nr_models+1)])))

    try:
        mod = Pandora.Pandora(target, template=template)
        # mod.model( output_dir=PANDORA.PANDORA_data + '/outputs/' + filename.replace('.csv', ''),
        #            stdev=0.2, benchmark=True, best_n_templates=3, n_homology_models=1,
        #            n_loop_models=nr_models, loop_refinement='slow')
        mod.model( output_dir=PANDORA.PANDORA_data + '/outputs/' + filename.replace('.csv', ''),
                   benchmark=True, n_loop_models=nr_models)



        dope = ','.join([str(round(float(i.dope), 4)) for i in mod.results])
        moldpdf = ','.join([str(round(float(i.moldpf), 4)) for i in mod.results])

        try:
            lmrsd = ','.join([str(round(i.lrmsd, 4)) for i in mod.results])
        except AttributeError:
            lmrsd = ','*(nr_models-1)
        try:
            core_lmrsd = ','.join([str(round(i.core_lrmsd, 4)) for i in mod.results])
        except AttributeError:
            core_lmrsd = ','*(nr_models-1)

        if len(mod.results) == nr_models:
            dope = dope + ','
            moldpdf = moldpdf + ','
            lmrsd = lmrsd + ','
            core_lmrsd = core_lmrsd + ','

        with open(PANDORA.PANDORA_path + '/../' + filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            mod.target.id, mod.target.peptide, '|'.join(mod.target.allele_type), '|'.join([i.id for i in mod.template]),
            '|'.join([i.peptide for i in mod.template]),
            ';'.join(['|'.join([x for x in i.allele_type]) for i in mod.template]), moldpdf, dope, lmrsd, core_lmrsd))
    except:
        print('Something went wrong')

def run_multiprocessing(func, i, num_cores):
    with Pool(processes=num_cores) as pool:
        return pool.map(func, i)

def bench_MHCII():
    t0 = time.time()
    # print(t0)
    db = Database.Database().load('25_04_21_Pandora_db')
    num_cores = 8

    list_of_targets_templates = []
    for k in db.MHCII_data:
        try:
            tar = PMHC.Target(db.MHCII_data[k].id, db.MHCII_data[k].allele_type, db.MHCII_data[k].peptide,
                                               M_chain_seq=db.MHCII_data[k].M_chain_seq,
                                               N_chain_seq=db.MHCII_data[k].N_chain_seq,
                                               MHC_class=db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
            mod = Pandora.Pandora(tar, db)
            mod.find_template(benchmark=True, best_n_templates=3)
            list_of_targets_templates.append((tar, mod.template))
        except:
            pass

        # run_pandora(list_of_targets_templates[0])

    run_multiprocessing(run_pandora, list_of_targets_templates, num_cores)

    print(time.time()-t0)

def bench_MHCI():
    t0 = time.time()
    # print(t0)
    db = Database.Database().load('25_04_21_Pandora_db')
    num_cores = 4

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
    # bench_MHCII()
    bench_MHCI()








