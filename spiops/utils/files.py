import os
import glob
import platform
import subprocess
import fnmatch
from tempfile import mkstemp
from shutil import move
from ftplib import FTP


# SOME CONSTANTS, consider move it to a config file
SERVER_HOST_NAME = "spiops.n1data.lan"
SERVER_HOST_USER = "esaspice"
SERVER_HOST_FTP_PATH = "/home/esaspice/ftp"
FTP_HOST = "spiftp.esac.esa.int"


def replace(file_path, pattern, subst):

    replaced = False

    #Create temp file
    fh, abs_path = mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:

                updated_line = line.replace(pattern, subst)
                new_file.write(updated_line)
                #flag for replacing having happened
                if updated_line != line:
                    replaced = True
    #Remove original file

    if replaced:

        os.remove(file_path)
        # Update the permissions
        os.chmod(abs_path, 0o644)
        #Move new file
        move(abs_path, file_path)

        return True

    return False


def mk2list(mk):

    path_symbol = ''
    ker_mk_list = []
    with open(mk, 'r') as f:
        for line in f:

            if path_symbol:
                if path_symbol in line:

                    kernel = line.split(path_symbol)[1]
                    kernel = kernel.strip()
                    kernel = kernel[:-1]

                    ker_mk_list.append(kernel)

            if 'PATH_SYMBOLS' in line.upper():
                path_symbol = '$' + line.split("'")[1]

    return ker_mk_list


def update_former_versions(mk_path, kernels_path, updated_mk=False):

    if not updated_mk:
        updated_mk = []
    #
    # Check the meta-kernels in the former versions directory. And
    # update them if needed.
    #
    if os.path.isdir(os.path.join(mk_path, 'former_versions')):

        os.chdir(os.path.join(mk_path, 'former_versions'))
        mks_in_dir = glob.glob('*.TM')
        mks_in_dir += glob.glob('*.tm')

        for mk_in_dir in mks_in_dir:
            updated_mk_flag = replace(mk_in_dir, "'..'", "'../..'")
            kernels_list = mk2list(mk_in_dir)

            for kernel in kernels_list:
                #
                # For the time being we only consider that non-present kernels
                # are in former versions.
                #
                ker_path = os.path.join(kernels_path, kernel[1:])
                if not os.path.exists(
                        ker_path) and 'former_versions' not in ker_path:
                    kernel_as_list = kernel.split('/')
                    kernel_former = os.path.join(
                            '/' + kernel_as_list[1] + '/former_versions',
                            kernel_as_list[2])
                    updated_mk_flag = replace(mk_in_dir, kernel,
                                              kernel_former)

            if updated_mk_flag:
                updated_mk.append(mk_in_dir)

            updated_mk_flag = False

    return updated_mk


def download_file(path, file):
    try:
        if ping(SERVER_HOST_NAME):
            path = os.path.join(SERVER_HOST_FTP_PATH, path)
            path = SERVER_HOST_USER + "@" + SERVER_HOST_NAME + ":" + path
            get_from_server(path, file)
        else:
            download_from_ftp(path, file)
    except:
        print('Warning: Error downloading file: ' + file)


def download_from_ftp(path, file):
    ftp = FTP(FTP_HOST)
    ftp.login()
    ftp.cwd(path)
    handle = open(file, 'wb')
    ftp.retrbinary('RETR %s' % file, handle.write)
    handle.close()
    return


def get_from_server(path, file):
    os.system('scp ' + os.path.join(path, file) + ' ./' + file)
    return


def ping(host):
    """
    Returns True if host (str) responds to a ping request.
    Remember that a host may not respond to a ping (ICMP) request even if the host name is valid.
    From: https://stackoverflow.com/questions/2953462/pinging-servers-in-python
    """

    # Option for the number of packets as a function of
    param = '-n' if platform.system().lower()=='windows' else '-c'

    # Building the command. Ex: "ping -c 1 google.com"
    command = ['ping', param, '1', host]

    return subprocess.call(command) == 0


def list_files_from_ftp(path, file_expression):
    ftp = FTP(FTP_HOST)
    ftp.login()
    ftp.cwd(path)
    files = []

    try:
        files = ftp.nlst()
    except FTP.error_perm as resp:
        if str(resp) == "550 No files found":
            print("No files in this directory")
        else:
            raise

    selected_files = []
    for file in files:
        if fnmatch.fnmatch(file, file_expression):
            selected_files.append(file)

    return selected_files
