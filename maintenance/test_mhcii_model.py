from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.load()

## B. Create Target object
target = PMHC.Target(id = 'MHCII_testcase',
    MHC_class='II',
    allele_type = ['HLA-DPA1*01:03','HLA-DPB1*01:01'],
    peptide = 'GSDWRFLRGYHQYA',
    use_netmhcpan=True)

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(n_loop_models=20, fully_flexible=True, stdev=0.3)