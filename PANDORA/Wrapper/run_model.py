#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

    target = args['target']
    template = args['template']
    #output_dir = args['outdir']
    n_loop_models = args['n_loop_models']
    n_jobs = args['n_jobs']
    benchmark = args['benchmark']
    pickle_out = args['pickle_out']
    clip_C_domain = args['clip_C_domain']


    # Create Pandora Object
    if args['outdir'] != '':
        output_dir = args['outdir']
    elif args['outdir'] == '' and args['collective_output_dir']:
        output_dir = args['collective_output_dir']
    else:
        output_dir = False

    mod = Pandora.Pandora(target, template=template,
                          output_dir=output_dir)

    # Run the modelling
    mod.model(n_loop_models=n_loop_models, n_jobs=n_jobs,
              stdev=0.1, benchmark=benchmark, pickle_out=pickle_out,
              clip_C_domain=clip_C_domain)
