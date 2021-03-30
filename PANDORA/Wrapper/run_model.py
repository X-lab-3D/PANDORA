#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 22:15:39 2021

@author: Dario Marzella
"""
import os
from PANDORA.Pandora import Pandora

def run_model(args):
    """
    

    Args:
        target (TYPE): DESCRIPTION.
        template (TYPE): DESCRIPTION.
        num_models (TYPE, optional): DESCRIPTION. Defaults to 20.
        benchmark (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        None.

    """
    target = args[0]
    template = args[1]
    n_models = args[2]
    n_jobs=args[3]
    benchmark = args[4]
    #print(' JOB: ', target, template, num_models, benchmark)
    '''
    if benchmark:
        filename = '220321_benchmark_I.csv'
        if not os.path.exists(filename):
            with open(filename, 'w') as f:
                f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s\n' % (
                ','.join(['molpdf_' + str(i + 1) for i in range(num_models)]),
                ','.join(['L-RMSD_' + str(i + 1) for i in range(num_models)]),
                ','.join(['core_L-RMSD_' + str(i + 1) for i in range(num_models)])))
    '''
    
    #try:
    mod = Pandora.Pandora(target, template = template)
    mod.model(n_models=n_models, n_jobs=n_jobs, stdev=0.1, benchmark=benchmark)

    '''
    moldpdf = ','.join([str(round(float(i.moldpf), 4)) for i in mod.results])
    if benchmark:
        lmrsd = ','.join([str(round(i.lrmsd, 4)) for i in mod.results])
        core_lmrsd = ','.join([str(round(i.core_lrmsd, 4)) for i in mod.results])

        with open(filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide,
            ';'.join(mod.template.allele_type), moldpdf, lmrsd, core_lmrsd))
    else:
        with open(filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s\n' % (
            mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide,
            ';'.join(mod.template.allele_type), moldpdf))
    '''
    #except:
    #    print('Something went wrong with the parallel job')