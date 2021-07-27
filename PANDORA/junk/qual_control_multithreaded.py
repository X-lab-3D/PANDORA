import sys
sys.path.append('./')

import os
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.multiprocessing import freeze_support
import time
import PANDORA



def run_qc(pdb_path):
    # pdb_path =all_model_paths[1]

    print(pdb_path)

    #     do something

    # t0 = time.time()
    try:
        os.system('bash /Users/derek/MolProbity/MolProbity-master/test/simple_molprobity_custom.sh %s' %(pdb_path))
    except:
        print('something went wrong')
        pass
    # print(time.time()-t0)

def run_multiprocessing(func, i, num_cores):
    with Pool(processes=num_cores) as pool:
        return pool.map(func, i)

def qual_control():
    t0 = time.time()

    num_cores = 12

    pandora_out_path = PANDORA.PANDORA_data + '/outputs/250421_1mod20_sd02_slow_templ1_benchmark_II/'

    # get a list of all PANDORA output dirs
    pandora_out_dirs = os.listdir(pandora_out_path)
    pandora_out_dirs = [i for i in pandora_out_dirs if not i.startswith('.')]

    # loop through all dirs
    all_model_paths = []
    for i in pandora_out_dirs:
        all_model_paths.extend([i + '/' + p for p in os.listdir(pandora_out_path + i) if p.endswith('.pdb')])

    # filter out template pdbs, only keep models
    all_model_paths = [i for i in all_model_paths if i.count('.') >= 2]

    all_model_paths = [pandora_out_path + i for i in all_model_paths]


    # run the function
    run_multiprocessing(run_qc, all_model_paths, num_cores)

    print(time.time()-t0)


if __name__ == "__main__":
    # Parallel(n_jobs=num_cores)(delayed(run_pandora)(target) for target in list_of_targets[:2])
    freeze_support()  # required to use multiprocessing
    qual_control()
    # bench_MHCI()