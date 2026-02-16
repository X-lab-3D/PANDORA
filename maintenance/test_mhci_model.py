from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.load()

## B. Create Target object
target = PMHC.Target(id = 'MHCI_testcase',
    allele_type = 'HLA-A*02:01:48', mhc_class = 'I',
    peptide = 'GILGFVFTL', anchors=[0,8],
    use_netmhcpan=False, rm_netmhcpan_output=False,
    )

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(n_loop_models=20, clip_C_domain=True, benchmark=False)
