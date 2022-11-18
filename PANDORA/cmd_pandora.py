import argparse

from PANDORA import Target
from PANDORA import Pandora
from PANDORA import Database

def cmd_install_database():
    """Command line tool to run install_database.
    """    
    parser = argparse.ArgumentParser(
    prog="PANDORA",
    description="PANDORA fetch database from zotero",
    )

    parser.add_argument(
        '-d','--destination-path', default='~/PANDORA_databases/default',
        help='Full path to the database destination'
    )

    args = parser.parse_args()

    Database.install_database(db_path=args.destination_path)

def run_pandora():
    """Command line tool to run one pandora case
    """
    parser = argparse.ArgumentParser(
    prog="PANDORA",
    description="Run one pMHC modelling case",
    )

    parser.add_argument(
        '-p','--peptide', required=True, type=str,
        help='One-letter sequence of the target peptide.'
    )

    parser.add_argument(
        '-a','--allele', required=True, type=str,
        help='Name of the target MHC allele'
    )
    
    parser.add_argument(
        '-c','--mhc-class', required=True, type=str,
        help='MHC class', choices=['I','II']
    )

    parser.add_argument(
        '-i','--id', 
        help='ID of the target case.\
             If not provided, it will default to peptide_allele'
    )

    parser.add_argument(
        '-k','--anchors', 
        help='Peptide anchor positions.\
            To be provided as series of integers separated by comma only.\
            Example: 2,9 for pMHC-I or 4,7,9,12 for pMHC-II'
    )

    parser.add_argument(
        '-o','--output-path', default='./',
        help='Output folder.'
    )

    args = parser.parse_args()

    #Parse arguments
    if not args.id:
        args.id = f'{args.peptide}_{args.allele}'
    
    

    ## A. Load local Database
    db = Database.load()

    ## B. Create Target object
    target = Target(id = args.id,
        allele_type = args.allele,
        peptide = args.peptide,
        anchors = [2,9])

    ## C. Perform modelling
    case = Pandora.Pandora(target, db)
    case.model()

