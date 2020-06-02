#bash rmsd_test.sh ref.pdb query_1ogt_AC.BL00090007.pdb
ref_FN=$1
mob_FN=$2

/home/dariomarzella/ProFitV3.3/bin/profit<< END
        reference $ref_FN
        mobile $mob_FN
        atom CA,C,N,O
        `cat ref.lzone`
        fit
        quit
END

