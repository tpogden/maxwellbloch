#!/usr/bin/env python
# coding: utf-8

""" Solves a MaxwellBloch problem given an input JSON file.
"""

import sys
import argparse

from maxwellbloch import mb_solve

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True,
        help='Path of a JSON file containing an mb_solve problem definition')
    parser.add_argument('-s', '--step', required=False, 
        help="Finite difference step to use. 'ab' or 'euler'",
            choices=['ab', 'euler'], default='ab')
    parser.add_argument('-r', '--recalc', required=False, action="store_true",
        help='Recalculate the solution even if a savefile exists?', 
        default=False)
    parser.add_argument('-p', '--pbarchunksize', required=False, type=int,
        help="How often should the progress bar be updated? Default is 10, to\
            update every 10%.", default=10)
    args = vars(parser.parse_args())

    print('Loading problem definition from file {0}'.format(args['file']))
    mbs = mb_solve.MBSolve().from_json(args['file'])
    print(mbs.to_json_str())

    if mbs.savefile_exists() and not args['recalc']:
        print(f'Savefile {mbs.savefile}.qu exists and recalc not requested.')
        print(f'Loading from file.')

    mbs.mbsolve(step=args['step'], rho0=None, recalc=args['recalc'],
                         pbar_chunk_size=args['pbarchunksize'])

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
