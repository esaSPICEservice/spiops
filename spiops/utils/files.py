import os
import glob
import platform
import subprocess
import fnmatch
import logging
import select
from tempfile import mkstemp
from shutil import move
from ftplib import FTP
from spiceypy import spiceypy


# SOME CONSTANTS, consider move it to a config file
SERVER_HOST_NAME = "spiops.ebe.lan"
SERVER_HOST_USER = "esaspice"
SERVER_HOST_FTP_PATH = "/home/esaspice/ftp"
FTP_HOST = "spiftp.esac.esa.int"


def replace(file_path, pattern, subst):

    replaced = False

    # Create temp file
    fh, abs_path = mkstemp()
    with os.fdopen(fh, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:

                updated_line = line.replace(pattern, subst)
                new_file.write(updated_line)
                # flag for replacing having happened
                if updated_line != line:
                    replaced = True

    # Remove original file
    if replaced:
        os.remove(file_path)
        os.chmod(abs_path, 0o644)  # Update the permissions
        move(abs_path, file_path)  # Move new file

    return replaced


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
        if str(path).startswith("http"):
            get_from_url(path, file)
        else:
            if ping(SERVER_HOST_NAME):
                path = os.path.join(SERVER_HOST_FTP_PATH, path)
                path = SERVER_HOST_USER + "@" + SERVER_HOST_NAME + ":" + path
                get_from_server(path, file)
            else:
                download_from_ftp(path, file)
        print('Warning: Downloaded file ' + path + os.sep + file)
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


def get_from_url(url, file):
    os.system('wget ' + url + ' -O ' + file)
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


def get_aocs_quaternions(aocs_file):
    aocs_file_data = open(aocs_file)
    quats = []

    for line in aocs_file_data.readlines():
        data = line.replace('\n', '').replace(',', ' ').split()
        et = spiceypy.str2et(data[0].replace('Z', ''))

        quat_entry = [et, float(data[2]), float(data[3]), float(data[4]), float(data[5])]
        quats.append(quat_entry)

    return quats


def get_aem_quaternions(aem_file, def_time_system="TDB"):
    aem_file_data = open(aem_file)
    inside_data_section = False
    curr_time_system = def_time_system
    prev_et = 0
    quats = []

    for line in aem_file_data.readlines():
        if 'TIME_SYSTEM' in line:
            time_system = line.split("=")[1].strip()
            if time_system != curr_time_system:
                spiceypy.timdef('SET', 'SYSTEM', 10, time_system)
                curr_time_system = time_system

        elif 'DATA_START' in line:
            inside_data_section = True

        elif inside_data_section:

            if 'DATA_STOP' in line:
                # Data section stopped, continue with next line
                inside_data_section = False
                continue

            data = line.replace('\n', '').split()

            et = spiceypy.str2et(data[0])

            if prev_et == et:
                # In case of segment end and next segment start matches, take only the
                # next segment start, so remove last entry
                quats.pop()

            quat_entry = [et, float(data[4]), float(data[1]), float(data[2]), float(data[3])]
            quats.append(quat_entry)
            prev_et = et

    return quats


def get_tm_data(file, delimiter, columns, data_factors):
    file_data = open(file)
    data = []
    for line in file_data.readlines():
        fields = line.split(delimiter)

        et = spiceypy.str2et(fields[0].replace('Z', ''))
        line_data = [et]

        df_idx = 0
        for column in columns:
            value = float(fields[column]) * data_factors[df_idx]
            line_data.append(value)
            df_idx += 1

        data.append(line_data)

    return data


def get_csv_data(file, delimiter, columns):
    file_data = open(file)
    data = []
    for line in file_data.readlines():
        fields = line.split(delimiter)
        line_data = []
        for column in columns:
            line_data.append(fields[column])
        data.append(line_data)
    return data


# For each file, download it, add data to array, and remove it
def download_tm_data(files, path, delimiter, columns, data_factors):

    tm_data = []

    for file in files:

        download_file(path, file)

        if not os.path.isfile(file):
            print('File cannot be downloaded: ' + file)
            return None

        file_tm_data = get_tm_data(file, delimiter, columns, data_factors)
        tm_data += file_tm_data

        # Remove downloaded file
        os.remove(file)

    return tm_data


def get_kernels_from_mk(metakernel, k_type, ignore_strs):
    kernels = []
    k_type_path = '/' + k_type + '/'
    path = ""

    multiple_path_values, multiple_path_symbols = False, False
    with open(metakernel, 'r') as f:
        for line in f:
            if multiple_path_values and ')' in line:
                path.append(line.split("'")[1] + k_type_path)
                multiple_path_values = False
            if multiple_path_symbols and ')' in line:
                symbol.append(line.split("'")[1])
                multiple_path_symbols = False
            if k_type_path in line:

                ignore_strs_found = False
                for ignore_str in ignore_strs:
                    if ignore_str in line:
                        ignore_strs_found = True
                        break

                if not ignore_strs_found:
                    kernel_path = path[symbol.index(line.split('$')[1].split('/')[0])]
                    kernels.append(kernel_path + line.split(k_type_path)[-1].strip().split("'")[0])

            if 'PATH_VALUES' in line and '=' in line:
                if not ')' in line:
                    multiple_path_values = True
                path = [line.split("'")[1] + k_type_path]
            if 'PATH_SYMBOLS' in line and '=' in line:
                if not ')' in line:
                    multiple_path_symbols = True
                symbol = [line.split("'")[1]]

    kernels = list(reversed(kernels))
    return kernels

def search_pds_file(pds_root, filename):
    # Search recursive in the PDS root
    files = glob.glob(f'{pds_root}/**/{filename}', recursive=True)
    if len(files) > 0:
        return files[0]
    return None

def check_spice_path(spice_path):
    if os.path.isabs(spice_path) and \
            ((os.path.isfile(spice_path) and len(spice_path) > 100)
                    or (os.path.isdir(spice_path) and len(spice_path) > 60)):
        spice_path = os.path.relpath(spice_path, os.getcwd())

    return spice_path

def read_all_text(file_path):
    text = ""
    with open(file_path, 'r') as d:
        lines = d.readlines()
        text = "".join(lines)
    return text

def commnt_read(kernel_path, directories=None):

    exe_path = directories.executables + '/' if directories is not None else ""

    command_line_process = subprocess.Popen([exe_path + 'commnt', '-r',
                                             check_spice_path(kernel_path)],
                                             bufsize=8192, shell=False,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
    process_output = ""
    dataend = False
    while (command_line_process.returncode is None) or (not dataend):
        command_line_process.poll()
        dataend = False

        ready = select.select([command_line_process.stdout, command_line_process.stderr], [], [], 1.0)

        if command_line_process.stderr in ready[0]:
            data = command_line_process.stderr.read(1024)
            if len(data) > 0:
                raise Exception(data)

        if command_line_process.stdout in ready[0]:
            data = command_line_process.stdout.read(1024)
            if len(data) == 0:  # Read of zero bytes means EOF
                dataend = True
            else:
                process_output += data.decode('utf-8')

    return process_output

def get_kernel_comments(kernel_path, directories=None):

    # First check if is a binary kernel or text kernel
    name, extension = os.path.splitext(kernel_path)
    if str(extension).lower().startswith(".t"):
        # Is a text kernel, just read the file:
        kernel_comments = read_all_text(kernel_path)
    else:
        # Shall be a binary kernel, use commnt read:
        kernel_comments = commnt_read(kernel_path, directories)

    return kernel_comments

def get_section_text_from_kernel_comments(kernel_path, section_name, directories=None):
    try:

        kernel_comments = get_kernel_comments(kernel_path, directories)
        kernel_comments = kernel_comments.splitlines()

        section_text = ""
        inside_section = False
        previous_line = ""
        for line in kernel_comments:

            if inside_section:

                if line.startswith("----") or line.startswith("===="):
                    break

                if not (previous_line.startswith("----") or previous_line.startswith("====")):
                    section_text += previous_line + "\n"

            elif previous_line.startswith(section_name) \
                    and (line.startswith("----") or line.startswith("====")):
                inside_section = True

            previous_line = line

        return section_text

    except Exception as ex:
        logging.info('Error on get_section_text_from_kernel_comments:', ex)
        return ""
