from spiops.utils.utils import get_latest_kernel, get_exe_dir, get_skd_path
from spiops.utils.utils import get_sc
import subprocess
import os

UTILITY_DIR = os.path.dirname(__file__) + '/..' + get_exe_dir()


def brief(kernel, utc=False):

    utility = UTILITY_DIR + os.sep + 'brief'
    option = '-c'

    skd_path = get_skd_path(kernel)

    try:
        lsk = get_latest_kernel('lsk', skd_path, 'naif????.tls')
    except:
        lsk = get_latest_kernel('lsk', skd_path, 'NAIF????.TLS')

    if utc:
        option += ' -utc'
        kernel += ' ' + skd_path + '/lsk/' + lsk

    command_line_process = subprocess.Popen([utility, option, kernel],
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT)

    process_output, _ = command_line_process.communicate()

    return process_output.decode("utf-8")


def ckbrief(kernel, utc=False):

    utility = UTILITY_DIR + os.sep + 'ckbrief'
    option = '-rel -n'

    skd_path = get_skd_path(kernel)
    sc = get_sc(kernel)

    try:
        lsk = get_latest_kernel('lsk', skd_path, 'naif????.tls')
    except:
        lsk = get_latest_kernel('lsk', skd_path, 'NAIF????.TLS')

    try:
        sclk = get_latest_kernel('sclk', skd_path, '{}_step_????????.tsc'.format(sc))
    except:
        sclk = get_latest_kernel('lsk', skd_path, '{}_STEP_????????.TSC'.format(sc.upper()))

    try:
        fk = get_latest_kernel('fk', skd_path, '{}_v??.tf'.format(sc))
    except:
        fk = get_latest_kernel('fk', skd_path, '{}_V??.TF'.format(sc.upper()))

    if utc:
        option += ' -utc'

    kernel += ' ' + skd_path + '/lsk/' + lsk
    kernel += ' ' + skd_path + '/sclk/' + sclk
    kernel += ' ' + skd_path + '/fk/' + fk

    command_line_process = subprocess.Popen([utility, option, kernel],
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT)

    process_output, _ = command_line_process.communicate()

    return process_output.decode("utf-8")


def optiks(mkernel, utc=False):

    if 'MEX' in mkernel:
        mission = 'MEX'
    elif 'VEX' in mkernel:
        mission = 'VEX'
    elif 'ROS' in mkernel:
        mission = 'ROS'
    elif 'em16' in mkernel:
        mission = 'TGO'
    elif 'bc' in mkernel:
        mission = 'MPO'
    elif 'JUICE' in mkernel:
        mission = 'JUICE'
    else:
        raise ValueError('OPTIKS utility could not run')

    utility = UTILITY_DIR + os.sep + 'optiks'
    option = '-half -units degrees -frame {}_SPACECRAFT -showfovframes'.format(mission)

    if utc:
        option += ' -epoch {}'.format(utc)

    print(option)

    command_line_process = subprocess.Popen([utility, option, mkernel],
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT)

    process_output, _ = command_line_process.communicate()

    return process_output.decode("utf-8")
