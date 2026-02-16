from PANDORA.Database import Database
from PANDORA.Wrapper import Wrapper
import time
import sys

n_cores = int(sys.argv[1])

t0 = time.time()
## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load()

t1 = time.time()
## B. Create the wrapper object

t2 = time.time()
## C. Create all Target Objects based on peptides in the .tsv file
wrap = Wrapper.Wrapper('/mnt/csb/PANDORA_MHCII_paper/data/cross_val_with_anchors_MHCII.tsv', db, MHC_class='II', 
                    IDs_col=0, peptides_col=1, allele_name_col=2,
                    anchors_col=3,  benchmark=True, verbose=False, 
                    header=False, num_cores=n_cores, n_loop_models=20,
                    M_chain_col=4, N_chain_col=5, restraints_stdev=False, pickle_out=True,
                    collective_output_dir='/mnt/csb/PANDORA_MHCII_paper/data/fully_flexible_pkl/mhcii_noflex')

## C. Perform modelling
t3 = time.time()

print(f'Load db {t1-t0} seconds')
print(f'Create Wrapper {t2-t1} seconds')
print(f'Run modelling {t3-t2} seconds')