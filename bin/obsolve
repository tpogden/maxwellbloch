#!/usr/bin/env python
# coding: utf-8

""" Solves a MaxwellBloch Optical Bloch problem given an input JSON file.
"""

import sys
import argparse

from maxwellbloch import ob_solve

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True,
        help='Path of a JSON file containing an ob_solve problem definition')
    parser.add_argument('-r', '--recalc', required=False, action="store_true",
        help='Recalculate the solution even if a savefile exists?', 
        default=False)
    args = vars(parser.parse_args())

    print('Loading problem definition from file {0}'.format(args['file']))
    obs = ob_solve.OBSolve().from_json(args['file'])
    print(obs.to_json_str())
    obs.solve(show_pbar=True, recalc=args['recalc'])

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
