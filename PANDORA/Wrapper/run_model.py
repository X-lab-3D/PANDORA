#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 22:15:39 2021

@author: Dario Marzella
"""
from PANDORA.Pandora import Pandora

def run_model(args):
    """
    

    Args:
        target (TYPE): DESCRIPTION.
        template (TYPE): DESCRIPTION.
        n_loop_models (TYPE, optional): DESCRIPTION. Defaults to 20.
        benchmark (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        None.

    """
    target = args[0]
    template = args[1]
    n_loop_models = args[2]
    n_jobs=args[3]
    benchmark = args[4]
    
    # Check if the output directory is provided
    if len(args) == 5:
        output_dir = False
    elif len(args) == 6:
        output_dir = args[5]
       
    # Create Pandora Object
    if output_dir:
        mod = Pandora.Pandora(target, template = template, output_dir=output_dir)
    else:
        mod = Pandora.Pandora(target, template = template)

    # Run the modelling
    mod.model(n_loop_models=n_loop_models, n_jobs=n_jobs, stdev=0.1, benchmark=benchmark)