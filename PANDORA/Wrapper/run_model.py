#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 22:15:39 2021

@author: Dario Marzella
"""
import os
from PANDORA import Pandora

def run_model(target, template, benchmark=False, num_models=20):
    target = target
    template = template
    nr_models = 20
    if benchmark:
        filename = '220321_benchmark_I.csv'
        if not os.path.exists(filename):
            with open(filename, 'w') as f:
                f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,%s,%s,%s\n' % (
                ','.join(['molpdf_' + str(i + 1) for i in range(nr_models)]),
                ','.join(['L-RMSD_' + str(i + 1) for i in range(nr_models)]),
                ','.join(['core_L-RMSD_' + str(i + 1) for i in range(nr_models)])))
    
    try:
        mod = Pandora.Pandora(target, template = template)
        mod.model(n_models=nr_models, stdev=0.1, benchmark=benchmark)

        moldpdf = ','.join([str(round(float(i.moldpf), 4)) for i in mod.results])
        lmrsd = ','.join([str(round(i.lrmsd, 4)) for i in mod.results])
        core_lmrsd = ','.join([str(round(i.core_lrmsd, 4)) for i in mod.results])

        with open(filename, 'a') as f:
            f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            mod.target.id, mod.target.peptide, ';'.join(mod.target.allele_type), mod.template.id, mod.template.peptide,
            ';'.join(mod.template.allele_type), moldpdf, lmrsd, core_lmrsd))
    except:
        print('Something went wrong')