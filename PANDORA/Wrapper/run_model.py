#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tarfile
import subprocess
import traceback
import glob
import os

from PANDORA.Pandora import Pandora

def archive_and_remove(case):       
    # create archive of the folder
    with tarfile.open(f'{case}.tar', 'w') as archive:
        case_files = glob.glob(os.path.join(case, '*'))
        for case_file in case_files:
            archive.add(case_file)
    # remove the original files from the folder
    if os.path.exists(f'{case}.tar'):
        subprocess.check_call(f"rm -r {case}", shell=True)
    else:
        print(f'Error creating archive: {case}.tar, skipping the file removal')

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
    archive_output = args['archive_output']


    # Create Pandora Object
    if args['outdir'] != '':
        output_dir = args['outdir']
    elif args['outdir'] == '' and args['collective_output_dir']:
        output_dir = args['collective_output_dir']
    else:
        output_dir = False
    
    try:
        mod = Pandora.Pandora(target, template=template,
                            output_dir=output_dir)
    except Exception as e:
        print(f"Modelling case {target.id} failed at Pandora object creation step")
        print(f"Captured error: {e}")
        print(traceback.format_exc())

    # Run the modelling
    try:
        mod.model(n_loop_models=n_loop_models, n_jobs=n_jobs,
                stdev=0.1, benchmark=benchmark, pickle_out=pickle_out,
                clip_C_domain=clip_C_domain)

    except Exception as e:
        print(f"Modelling case {target.id} failed at modelling step")
        print(f"Captured error: {e}")
        print(traceback.format_exc())

    try:
        if archive_output:
            archive_and_remove(mod.output_dir)
    except Exception as e:
        print(f"Modelling case {target.id} failed at archiving step")
        print(f"Captured error: {e}")
        print(traceback.format_exc())

