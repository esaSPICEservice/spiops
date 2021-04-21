#!/usr/bin/env python3

import os
import textwrap
import glob
import re

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from spiops.utils.utils import replace
from spiops.utils.files import update_former_versions


def main(test=False, log=False):
    execution_dir = os.getcwd()

    with open(os.path.dirname(__file__) + '/config/version',
              'r') as f:
        for line in f:
            version = line

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                                description=textwrap.dedent('''\

 SPIOPS -- Version {}, SPICE Operational Procedures for ESA Missions

   SPIOPS is a library aimed to help scientists and engineers that deal 
   with Solar System Geometry for ESA planetary science missions. More 
   information is available here:

      https://github.com/esaSPICEservice/spiops


'''.format(version)),
                                epilog='''
    __   __   __      __   __     __   ___     __   ___  __          __   ___
   /__\ /__` '__\    /__` |__) | /  ` |__     /__` |__  |__) \  / | /  ` |__
   \__, .__/ \__/    .__/ |    | \__, |___    .__/ |___ |  \  \/  | \__, |___

 esa_spice@sciops.esa.int
 http://spice.esac.esa.int

''')
    parser.add_argument('-v', '--version',
                        help='Display the version of SPIOPS',
                        action='store_true')
    parser.add_argument('-m', '--metakernel',
                        help='Generate local meta-kernels from {mission}_{type}.tm',
                        action='store_true')
    parser.add_argument('-a', '--all',
                        help='Generate local meta-kernels from all files',
                        action='store_true')
    parser.add_argument('-c', '--clean',
                        help='Remove local meta-kernels',
                        action='store_true')
    parser.add_argument('-f', '--former',
                        help='Update the meta-kernels in the former_versions directory.',
                        action='store_true')
    args = parser.parse_args()

    if args.version:
        print(version)
        return

    if args.clean:

        cwd = os.getcwd()
        mks_in_dir = glob.glob('*local*.tm')
        mks_in_dir += glob.glob('*LOCAL*.TM')

        for mk_in_dir in mks_in_dir:
            os.remove(cwd + os.sep + mk_in_dir)

    if args.metakernel:

        cwd = os.getcwd()
        local_mks = []
        mks_in_dir = glob.glob('*.tm')
        mks_in_dir += glob.glob('*.TM')


        for mk_in_dir in mks_in_dir:
            if 'local' not in mk_in_dir.lower():
                if args.all:
                    replace(mk_in_dir, "'..'",
                            "'" + cwd.rsplit('/kernels', 1)[0] + "/kernels'")
                    local_mks.append(mk_in_dir)
                else:
                    not_append = re.search(r".*_v[0-9]{3}_[0-9]{8}_[0-9]{3}.tm", mk_in_dir.lower())
                    if not_append == None:
                        replace(mk_in_dir, "'..'",
                                "'" + cwd.rsplit('/kernels', 1)[0] + "/kernels'")
                        local_mks.append(mk_in_dir)

        if local_mks:
            print('SPIOPS -- Meta-Kernel Update\nThe following meta-kernels have been generated/replaced to local:')
            for mk in local_mks:
                print(f'   {mk}')
        else:
            print(
                'SPIOPS -- Meta-Kernel Update -- No meta-kernels have been updated.')

    if args.former:
        cwd = os.getcwd()
        if 'mk/former_versions' in cwd:
            mk_dir = os.sep.join(cwd.split(os.sep)[:-1])
            kernels_dir = os.sep.join(cwd.split(os.sep)[:-2])
        else:
            mk_dir = cwd
            kernels_dir = os.sep.join(cwd.split(os.sep)[:-1])
        try:
            updated_mks = update_former_versions(mk_dir, kernels_dir)
            if updated_mks:
                print(
                    'SPIOPS -- Meta-Kernel Update\nThe following former_versions meta-kernels have been updated:')
                for mk in updated_mks:
                    print(f'   {mk}')
            else:
                print(
                        'SPIOPS -- Meta-Kernel Update -- No meta-kernels have been updated.')
        except Exception as e: print(e)


    return

