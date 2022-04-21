#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 22:15:39 2021

@author: Dario Marzella
"""
from PANDORA.Pandora import Pandora

def run_model(args):
    """Runs one modelling job. Meant to be runned from Pandora.Wrapper
    
    Args:
        args (list): List of arguments. Should be containing the following, in
            order.
        target (Pandora.PMHC.PMHC.Target): Target object.
        template (Pandora.PMHC.PMHC.Template): Template object.
        n_loop_models (int, optional): Number of loop models. Defaults to 20.
        benchmark (bool, optional): Set True if running a benchmark to retrieve
            models RMSD with reference structures. Defaults to False.

    Returns:
        None.

    """
    # Check if the output directory is provided
    # if len(args) == 6:
    #     output_dir = False
    #     start=0
    # elif len(args) == 7:
    #     output_dir = args[0]
    #     start=1
    
    target = args[0]
    template = args[1]
    output_dir = args[2]
    n_loop_models = args[3]
    n_jobs=args[4]
    benchmark = args[5]
    pickle_out = args[6]
    
       
    # Create Pandora Object
    if output_dir != '':
        mod = Pandora.Pandora(target, template=template, 
                              output_dir=output_dir)
    else:
        mod = Pandora.Pandora(target, template = template)

    # Run the modelling
    mod.model(n_loop_models=n_loop_models, n_jobs=n_jobs,
              stdev=0.1, benchmark=benchmark, pickle_out=pickle_out)