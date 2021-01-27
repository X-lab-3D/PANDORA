#!/usr/bin/env python

"""
Select alternate location of residue.
Default is A.

usage: python pdb_selalt.py -<alt> <pdb file>
example: python pdb_alt.py -A 4FIV.pdb

Author: {0} ({1})

This program is part of the PDB tools distributed with HADDOCK
or with the HADDOCK tutorial. The utilities in this package
can be used to quickly manipulate PDB files, with the benefit
of 'piping' several different commands. This is a rewrite of old
FORTRAN77 code that was taking too much effort to compile. RIP.
"""

import os
import re
import sys

__author__ = "Cunliang Geng; Joao Rodrigues"
__email__ = "gengcunliang@gmail.com; j.p.g.l.m.rodrigues@gmail.com"

USAGE = __doc__.format(__author__, __email__)

def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options."""

    if not len(args):
        # No alt, from pipe
        if not sys.stdin.isatty():
            pdbfh = sys.stdin
            alt = 'A'
        else:
            sys.stderr.write(USAGE)
            sys.exit(1)
    elif len(args) == 1:
        # Alt & Pipe _or_ file & no alt
        if re.match('\-[A-Za-z0-9]', args[0]):
            alt = args[0][1:]
            if not sys.stdin.isatty():
                pdbfh = sys.stdin
            else:
                sys.stderr.write(USAGE)
                sys.exit(1)
        else:
            if not os.path.isfile(args[0]):
                sys.stderr.write('File not found: ' + args[0] + '\n')
                sys.stderr.write(USAGE)
                sys.exit(1)
            pdbfh = open(args[0], 'r')
            alt = 'A'
    elif len(args) == 2:
        # alt & File
        if not re.match('\-[A-Za-z0-9]', args[0]):
            sys.stderr.write('Invalid alt ID: ' + args[0] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)

        if not os.path.isfile(args[1]):
            sys.stderr.write('File not found: ' + args[1] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        alt = args[0][1:]
        pdbfh = open(args[1], 'r')
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)

    return (alt, pdbfh)

def _sel_alt(fhandle, alt_id):
    """Enclosing logic in a function to speed up a bit"""

    coord_re = re.compile('^(ATOM|HETATM)')
    fhandle = fhandle
    alt_id = alt_id

    for line in fhandle:
        if coord_re.match(line):
            if line[16]==" ":
                yield line
            elif line[16]== alt_id[0]:
                yield line[:16] + " " + line[17:]
        else:
            yield line

if __name__ == '__main__':
    # Check Input
    alt, pdbfh = check_input(sys.argv[1:])

    # Do the job
    new_pdb = _sel_alt(pdbfh, alt)

    try:
        sys.stdout.write(''.join(new_pdb))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)