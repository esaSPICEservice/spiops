#!/usr/bin/env python3

import math
import spiceypy
import logging
import os
import numpy as np

try:
    # Try import SpiceNOFRAMECONNECT exception from spiceypy 3.1.1
    from spiceypy.utils.exceptions import SpiceNOFRAMECONNECT as SpiceNOFRAMECONNECT
except ImportError:
    # Otherwise consider a SpiceNOFRAMECONNECT exception as an SpiceyError, for spiceypy 2.3.2
    from spiceypy.utils.support_types import SpiceyError as SpiceNOFRAMECONNECT

from spiceypy.utils.support_types import *


from spiops.utils import time
from spiops.utils.utils import plot
from spiops.utils.utils import plot_attitude_error
from spiops.utils.utils import target2frame
from spiops.utils.utils import findIntersection
from spiops.utils.utils import findNearest
from spiops.utils.files import download_file
from spiops.utils.files import list_files_from_ftp
from spiops.utils.files import get_aem_quaternions
from spiops.utils.files import get_aocs_quaternions
from spiops.utils.files import download_tm_data

from spiops.utils.naif import optiks  # Do not remove, called from spival
from spiops.utils.naif import brief  # Do not remove, called from spival

import imageio
import matplotlib as mpl
import matplotlib.pyplot as plt

from spiceypy import support_types as stypes

from bokeh.plotting import figure, output_file, output_notebook, show
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, LabelSet
from spiops.utils.time import et_to_datetime
# from spiops.utils.webmust.webmust_handler import WebmustHandler


"""
The MIT License (MIT)
Copyright (c) [2015-2017] [Andrew Annex]
Copyright (c) [2015-2017] [Andrew Annex]
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

def load(mk):
    return spiceypy.furnsh(mk)


def adcsng_fill_template(template,
                  file,
                  replacements,
                  cleanup=False):


   #
   # If the temp   late file is equal to the output file then we need to create a temporary template - which will be
   # a duplicate - in order to write in the file. A situation where we would like to have them be the same is
   # for example if we call this function several times in a row, replacing keywords in the template in steps
   #
   if template == file:
       with open(file, "r") as f:
           with open('fill_template.temp', "w+") as t:
               for line in f:
                   t.write(line)

       template = 'fill_template.temp'

   with open(file, "w+") as f:
       #
       # Items are replaced as per correspondance in between the replacements dictionary
       #
       with open(template, "r+") as t:
           for line in t:
               if '{' in line:
                   for k, v in replacements.items():
                       if '{' + k + '}' in line: line = line.replace('{' + k + '}', v)
               f.write(line)

               #
               # If the option cleanup is set as true, we remove the keyword assignments in the filled templated which are
               # unfilled (they should be optional)
               #
   if cleanup:

       with open(file, "r") as f:
           with open('fill_template.temp', "w+") as t:
               for line in f:
                   t.write(line)

       template = 'fill_template.temp'

       with open(file, "w+") as f:
           with open('fill_template.temp', "r") as t:
               for line in t:
                   if '{' not in line:
                       f.write(line)

                       #
                       # The temporary files are removed
                       #
   if os.path.isfile('fill_template.temp'):
       os.remove('fill_template.temp')


# Originally an adcsng funtion, needs to be re-arranged in adcsng to be made
# more generic
def adcsng_hk_quaternions2ck_reader(tm_file,
                                    input_time_format='UTC',
                                    input_time_field_number='1',
                                    delimiter=',',
                                    input_processing=False,
                                    qs_col=1, qx_col=2, qy_col=3, qz_col=4):

    #
    # We obtain the number of data fields and its correspondance
    #
    input_data_field_numbers = [qx_col, qy_col, qz_col, qs_col]

    tm_list = []
    previous_row_time = ''

    sclk_partition = '1'
    sclk_delimiter = '.'


    filter_flag = False
    index = 0
    row_prev = []
    sclk_fraction_prev = ''
    with open(tm_file, 'r') as t:

        for line in t:

            #
            # TODO: Main difference from fucntion from adcsng
            #
            if '#' not in line and 'Date' not in line and input_time_format not in line:
                index += 1

                row_data = []

                # We need to remove the end of line character:
                line = line.split('\n')[0]

                try:
                    if ',' in delimiter:

                        if input_time_format == 'SCLK':
                            if ',' in input_time_field_number:
                                row_time = sclk_partition + '/' + str(line.split(delimiter)[
                                                   int(input_time_field_number[0]) - 1]) + \
                                           sclk_delimiter + str(line.split(delimiter)[
                                                   int(input_time_field_number[2]) - 1])

                            else:
                                input_time_field_number = int(input_time_field_number)
                                row_time = str(line.split(delimiter)[
                                                input_time_field_number - 1])

                        else:
                            row_time = str(line.split(delimiter)[input_time_field_number-1])

                        if (' ' in row_time):
                            if input_time_format == 'SCLK':
                                row_time = row_time.replace(' ','')
                            else:
                                row_time = row_time.replace(' ','T')

                        for data_element_field_number in input_data_field_numbers:
                            row_data.append(float(line.split(',')[data_element_field_number-1]))

                    else:

                        proc_line = line.strip()

                        row_time = str(proc_line.split(delimiter)[input_time_field_number - 1])

                        for data_element_field_number in input_data_field_numbers:
                            #
                            # We need to check that
                            #
                            row_data.append(float(line.split()[data_element_field_number-1]))
                except:
                    logging.info('   HM TM Processing: Found incomplete data line in line {}:'.format(index))
                    logging.info('   {}'.format(line))
                    continue

                row = row_time + ' '

                # As indicated by Boris Semenov in an e-mail "ROS and MEX "measured" CKs"
                # sometimes the scalar value is negative and the sign of the rest of the
                # components of the quaternions needs to be changed!
                if row_data[-1] < 0:
                    neg_data = [-x for x in row_data]

                    logging.info('   HM TM Processing: Found negative QS on input line {}:'.format(row_data))
                    logging.info('   ' + neg_data)
                    row_data = neg_data

                for element in row_data:

                    row += str(element) + ' '

                # We filter out "bad quaternions"

                row += '\n'

                # We remove the latest entry if a time is duplicated
                if row_time == previous_row_time:
                    logging.info(
                        '   HM TM Processing: Found duplicate time at {}'.format(
                                row_time))
                else:
                    # We do not include the entry if one element equals 1 or gt 1
                    append_bool = True
                    for quaternion in row_data:
                        if quaternion >= 1.0:
                            append_bool = False
                            logging.info(
                                '   HM TM Processing: Found quaternion GT 1 on input line {}:'.format(
                                    row_data))
                            logging.info('   ' + str(row))

                    # This is a special filter that has been set for ExoMars2016
                    # More explanations in [1]
                    if input_processing:
                        sclk_fraction = line.split(':')[-1].split(' ')[0]

                        if filter_flag:
                            if sclk_fraction == sclk_fraction_prev:
                                row_prev.append(row)
                            elif len(row_prev) <= 5 and sclk_fraction == sclk_initial:

                                logging.info(
                                    '   HM TM Processing: Coarse quaternion: Spurious SCLK fractions before input line {}:'.format(
                                            index))

                                for element in row_prev:

                                    logging.info('   ' + str(element).split('\n')[0])
                                    tm_list.remove(element)

                                filter_flag = False
                                tm_list = []
                                row_prev = []
                                sclk_fraction_prev = sclk_fraction
                            else:
                                row_prev = []
                                filter_flag = False

                        if sclk_fraction_prev and sclk_fraction != sclk_fraction_prev and not filter_flag:
                                filter_flag = True
                                row_prev.append(row)
                                sclk_initial = sclk_fraction_prev

                        sclk_fraction_prev = sclk_fraction

                    if append_bool:
                        tm_list.append(row)

                previous_row_time = row_time

    # We remove the carriage return from the last line
    last_line = tm_list[-1].split('\n')[0]
    tm_list = tm_list[:-1]
    tm_list.append(last_line)

    return(tm_list)


def fov_illum(mk, sensor, time=None, angle='DEGREES', abcorr='LT+S',
              report=False, unload=False):
    """
    Determine the Illumination of a given FoV (for light scattering computations
    for example). This function is based on  the following SPICE APIs:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/getfov_c.html
    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/spkezp_c.html


    :param mk: Meta-kernel to load the computation scenario
    :type mk: str
    :param sensor: Sensor ID code or name
    :type sensor: Union[str, int]
    :param time: Time to compute the quantity
    :type time: Union[str, float]
    :param angle: Angular unit; it can be 'DEGREES' or 'RADIANS'. Default is 'DEGREES'
    :type angle: str
    :param abcorr: Aberration correction. Default and recommended is 'LT+S'
    :type abcorr: str
    :param report: If True prints the resulting illumination angle on the screen
    :type report: bool
    :param unload: If True it will unload the input meta-kernel
    :type unload: bool
    :return: Angle in between a sensor's boresight and the sun-sc direction
    :rtype: float
    """
    room = 99
    shapelen = 1000
    framelen = 1000
    angle = angle.upper()

    spiceypy.furnsh(mk)

    if time:
        time = spiceypy.utc2et(time)
    else:
        time = spiceypy.utc2et('2016-08-10T00:00:00')

    if angle != 'DEGREES' and angle != 'RADIANS':
        print('angle should be either degrees or radians')

    if isinstance(sensor, str):
        instid = spiceypy.bodn2c(sensor)
    else:
        instid = sensor

    shape, frame, bsight, n, bounds = spiceypy.getfov(instid, room, shapelen, framelen)

    rotation = spiceypy.pxform(frame, 'J2000', time)

    bsight = spiceypy.mxv(rotation, bsight)

    # The following assumes that the IDs of the given S/C FK have been defined
    # according to the NAIF/ESS standards:
    #
    #    -NXXX
    #
    #       where:
    #          N is the SC id and can consist on a given number of digits
    #          XXX are three digits that identify the sensor
    sc_id = int(str(instid)[:-3])

    ptarg, lt = spiceypy.spkezp(10, time, 'J2000', abcorr, sc_id)

    fov_illumination = spiceypy.vsep(bsight, ptarg)

    if unload:
        spiceypy.unload(mk)

    if angle == 'DEGREES':
        fov_illumination = math.degrees(fov_illumination)

    if report:
        print('Illumination angle of {} is {} [{}]'.format(sensor,
                                                           fov_illumination,
                                                           angle))

    return fov_illumination


def cov_spk_obj(mk, object, time_format='TDB', global_boundary=False,
                report=False, unload=False):
    """
    Provides time coverage summary for a given object for a list of
    binary SPK files provided in a meta-kernel. Several options are
    available. This function is based on the following SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/spkcov_c.html

    The NAIF utility BRIEF can be used for the same purpose.

    :param mk: Meta-kernel to load the computation scenario
    :type mk: str
    :param object: Ephemeris Object to obtain the coverage from
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB in calendar format) or 'TDB'. Default is 'TDB'
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the coverage windows or only the absolute start and finish coverage times
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen
    :type report: bool
    :param unload: If True it will unload the input meta-kernel
    :type unload: bool
    :return: Returns a list with the coverage intervals
    :rtype: list
    """
    spiceypy.furnsh(mk)
    boundaries_list = []
    et_boundaries_list = []

    object_id = spiceypy.bodn2c(object)
    maxwin = 2000
    spk_count = spiceypy.ktotal('SPK') - 1

    while spk_count >= 0:

        spk_kernel = spiceypy.kdata(spk_count, 'SPK', 155, 155, 155)

        spk_ids = spiceypy.spkobj(spk_kernel[0])

        for id in spk_ids:

            if id == object_id:

                object_cov = SPICEDOUBLE_CELL(maxwin)
                spiceypy.spkcov(spk_kernel[0], object_id, object_cov)

                boundaries = time.cov_int(object_cov=object_cov,
                                          object_id=object_id,
                                          kernel=spk_kernel[0],
                                          global_boundary=global_boundary,
                                          time_format=time_format,
                                          report=report)

                boundaries_list.append(boundaries)

                #
                # We need to have the boundaries in TDB in order to sort out the
                # min and max to obtain the global ones for multiple kernels
                #
                if global_boundary:
                    et_boundaries_list.append(time.cov_int(
                                                      object_cov=object_cov,
                                                      object_id=object_id,
                                                      kernel=spk_kernel[0],
                                                      global_boundary=True,
                                                      time_format='TDB',
                                                      report=False))

        spk_count -= 1

    if global_boundary:
        start_time = min(et_boundaries_list)[0]
        finish_time = max(et_boundaries_list)[1]

        boundaries_list = time.et2cal([start_time, finish_time],
                                      format=time_format)

        if report:
            print("Global Coverage for {} [{}]: {} - {}".format(
                str(spiceypy.bodc2n(object_id)), time_format, boundaries_list[0],
                boundaries_list[1]))


    if unload:
        spiceypy.unload(mk)

    return boundaries_list


def cov_spk_ker(spk, object=False, time_format='TDB', support_ker ='',
                report=False, unload=True):
    """
    Provides time coverage summary for a given object for a given SPK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/spkcov_c.html

    The NAIF utility BRIEF can be used for the same purpose.

    :param spk: SPK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least it should be a leapseconds kernel (LSK) and optionally a meta-kernel (MK)
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object or list of objects to obtain the coverage from
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' or 'SPICE' (for TDB in calendar format) or 'TDB'. Default is 'TDB'
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the coverage windows or only the absolute start and finish coverage times
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen
    :type report: bool
    :param unload: If True it will unload the input meta-kernel
    :type unload: bool
    :return: Returns a list with the coverage intervals
    :rtype: list
    """
    spiceypy.furnsh(spk)
    object_id = []
    boundaries = []

    if object and not isinstance(object, list):
        object = [object]

    if support_ker:

        if isinstance(support_ker, str):
            support_ker = [support_ker]

        for ker in support_ker:
            spiceypy.furnsh(ker)

    maxwin = 2000

    spk_ids = spiceypy.spkobj(spk)

    if not object:
        object_id = spk_ids
        object = []
        for id in spk_ids:
            object.append(spiceypy.bodc2n(id))
    else:
        for element in object:
            object_id.append(spiceypy.bodn2c(element))

    for id in object_id:

        if id in spk_ids:

            object_cov = SPICEDOUBLE_CELL(maxwin)
            spiceypy.spkcov(spk, id, object_cov)

            cov = time.cov_int(object_cov=object_cov,
                                      object_id=id,
                                      kernel=spk,
                                      time_format=time_format,
                                      report=report)

        else:
            if report:
                print('{} with ID {} is not present in {}.'.format(object,
                                                             id, spk))
            if unload:
                spiceypy.unload(spk)
                if support_ker:

                    if isinstance(support_ker, str):
                        support_ker = [support_ker]

                    for ker in support_ker:
                        spiceypy.unload(ker)
            return False

        if time_format == 'SPICE':
            boundaries.append(object_cov)
        else:
            boundaries.append(cov)

    if unload:
        spiceypy.unload(spk)
        if support_ker:

            if isinstance(support_ker, str):
                support_ker = [support_ker]

            for ker in support_ker:
                spiceypy.unload(ker)

    return boundaries


def spkVsOem(sc, spk, plot_style='line', notebook=True):

    spiceypy.timdef('SET', 'SYSTEM', 10, 'TDB')
    spiceypy.furnsh(spk)

    if sc == 'MPO':
        file = spk.split('/')[-1].replace('\n', '').replace('bc_mpo_fcp_', '').split('_')[0]
        file = 'BCCruiseOrbit__' + file + '.bc'
        download_file("data/ANCDR/BEPICOLOMBO/fdy", file)
    else:
        print('Unsupported spacecraft: ' + sc)
        return None, None

    print('OEM file: ' + file)
    if not os.path.isfile(file):
        print('OEM file cannot be downloaded!')
        return None, None

    oemfile = open(file)
    error = []
    pos_norm_error = []
    vel_norm_error = []
    data_list = []
    for line in oemfile.readlines():
        if 'CENTER_NAME' in line:
            center = line.split('= ')[1].replace('\n', '')
        if line[:2] == '20':
            data = line.replace('\n', '').split()
            data_list.append(data)
    for i in range(0, len(data_list)-1, 1):
        #
        # skip OEM lines with repeated time tags (typically at the end of a
        # segment) as are superseeded by the latest line with that time tag
        #
        if data_list[i][0] != data_list[i+1][0]:
            data = data_list[i]
            et = spiceypy.str2et(data[0])
            state = spiceypy.spkezr(sc, et, 'J2000', 'NONE', center)[0]
            curr_error = [et,
                          abs(state[0] - float(data[1])),
                          abs(state[1] - float(data[2])),
                          abs(state[2] - float(data[3])),
                          abs(state[3] - float(data[4])),
                          abs(state[4] - float(data[5])),
                          abs(state[5] - float(data[6]))]
            error.append(curr_error)
            pos = np.asarray(curr_error[1:4])
            vel = np.asarray(curr_error[4:7])
            pos_norm_error.append(spiceypy.vnorm(pos))
            vel_norm_error.append(spiceypy.vnorm(vel))

    max_pos_norm_error = max(pos_norm_error)
    max_vel_norm_error = max(vel_norm_error)

    error = np.asarray(error)
    print('Avg position error [km]: ' + str(np.mean(pos_norm_error)) +
          ' , max position error [km]: ' + str(max_pos_norm_error))
    print('Avg velocity error [km/s]: ' + str(np.mean(vel_norm_error)) +
          ' , max velocity error [km/s]: ' + str(max_vel_norm_error))

    plot(error[:, 0],
         [error[:, 1], error[:, 2], error[:, 3]],
         yaxis_name=['X', 'Y', 'Z'],
         title='Source OEM to generated SPK position difference',
         format=plot_style,
         yaxis_units='Position error Km',
         notebook=notebook)
    plot(error[:, 0],
         [error[:, 4], error[:, 5], error[:, 6]],
         yaxis_name=['VX', 'VY', 'VZ'],
         title='Source OEM to generated SPK velocity difference',
         format=plot_style,
         yaxis_units='Km/s',
         notebook=notebook)

    os.remove(file)
    spiceypy.timdef('SET', 'SYSTEM', 10, 'UTC')
    spiceypy.unload(spk)

    return max_pos_norm_error, max_vel_norm_error


def ckVsAEM(sc, ck, plot_style='line', notebook=True):

    spiceypy.timdef('SET', 'SYSTEM', 10, 'TDB')
    spiceypy.furnsh(ck)

    if sc == 'MPO':
        file = ck.split('/')[-1].replace('\n', '').split('_')[4]
        file = 'AttitudePredictionST__' + file + '.bc'
        download_file("data/ANCDR/BEPICOLOMBO/fdy", file)
    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    print('AEM file: ' + file)
    if not os.path.isfile(file):
        print('AEM file cannot be downloaded!')
        return None

    aem_guats = get_aem_quaternions(file)

    if len(aem_guats):
        # If any quaternion inserted, remove the first element to create a
        # margin with the start of the CK
        aem_guats.pop(0)

    error, max_ang_error = get_quats_ang_error(aem_guats, sc)

    plot_attitude_error(np.asarray(error),
                        max_ang_error,
                        'Source AEM Quaternions to generated CK orientation difference',
                        plot_style,
                        notebook)

    os.remove(file)
    spiceypy.unload(ck)

    return max_ang_error


def ckVsAocs(sc, ck, plot_style='line', notebook=True):

    spiceypy.timdef('SET', 'SYSTEM', 10, 'UTC')
    spiceypy.furnsh(ck)

    if sc == 'MPO':
        file = ck.split('/')[-1].replace('\n', '').split('_')[5]
        file = 'mpo_raw_hk_aocs_measured_attitude_' + file + '.tab'
        download_file("data/ANCDR/BEPICOLOMBO/hkt", file)
    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    print('AOCS tab file: ' + file)
    if not os.path.isfile(file):
        print('AOCS tab file cannot be downloaded!')
        return None

    aocs_quats = get_aocs_quaternions(file)

    error, max_ang_error = get_quats_ang_error(aocs_quats, sc)

    plot_attitude_error(error,
                        max_ang_error,
                        'Source AOCS Measured Quaternions to generated CK orientation difference',
                        plot_style,
                        notebook)

    os.remove(file)
    spiceypy.unload(ck)

    return max_ang_error


def get_quats_ang_error(quats, sc):
    error = []
    max_ang_error = 0

    for quat in quats:
        et = quat[0]
        q_spice = spiceypy.m2q(spiceypy.pxform('J2000', sc + '_SPACECRAFT', et))

        if quat[1] < 0:
            q_spice[0] *= -1
            q_spice[1] *= -1
            q_spice[2] *= -1
            q_spice[3] *= -1

        quat[2] *= -1
        quat[3] *= -1
        quat[4] *= -1

        q_error = [abs(q_spice[0] - quat[1]),
                   abs(q_spice[1] - quat[2]),
                   abs(q_spice[2] - quat[3]),
                   abs(q_spice[3] - quat[4])]

        mrot_spice = spiceypy.q2m(q_spice)
        mrot_quats = spiceypy.q2m(quat[1:5])

        vz_spice = spiceypy.mxv(mrot_spice, [0, 0, 1])
        vz_quats = spiceypy.mxv(mrot_quats, [0, 0, 1])

        ang_error = spiceypy.vsep(vz_spice, vz_quats)
        max_ang_error = max(max_ang_error, abs(ang_error))

        curr_error = [et]
        curr_error.extend(q_error)
        error.append(curr_error)

    max_ang_error = np.rad2deg(max_ang_error) * 1000
    return np.asarray(error), max_ang_error


def saa_vs_hk_sa_position(sc, plot_style='line', notebook=True):

    spiceypy.timdef('SET', 'SYSTEM', 10, 'UTC')

    sa_angles = []  # Read angles from TM, List of items as [et, angle_deg]

    if sc == 'MPO':

        # Set some mission specific constants
        sadm_frame = 'MPO_SA'               # SA Rotating Frame
        sadm_ref_frame = 'MPO_SA_SADM'      # SA Fixed Frame
        ref_vector = np.asarray([0, 0, 1])  # SA Rotating Plane normal
        ref_cross_vector = np.asarray([0, 1, 0])  # Common rotation vector btw SA Rotating Frm and Fixed Frm

        # For MPO SA the TM files are given in daily basis, so we need to
        # concatenate N of them to obtain a greater period coverage.
        hkt_path = "data/ANCDR/BEPICOLOMBO/hkt/"
        hkt_expression = 'mpo_raw_hk_sa_position_????????.tab'
        num_sa_files = 7  # Compare last week

        # Determine files to use to fetch TM data
        sa_files = list_files_from_ftp(hkt_path, hkt_expression)
        sa_files = sa_files[-num_sa_files:]

        # For each file, download it, add data to array, and remove it
        sa_angles = download_tm_data(sa_files, hkt_path, ",", [2], [180.0 / math.pi])
        if not len(sa_angles):
            print("Cannot obtain required TM data, aborting.")
            return None

    elif sc == 'MTM':

        # Set some mission specific constants
        sadm_frame = 'MTM_SA+X'               # SA Rotating Frame
        sadm_ref_frame = 'MTM_SA+X_ZERO'      # SA Fixed Frame
        ref_vector = np.asarray([0, 1, 0])    # SA Rotating Plane normal
        ref_cross_vector = np.asarray([1, 0, 0])  # Common rotation vector btw SA Rotating Frm and Fixed Frm

        # For MTM SA the TM files are given in daily basis, so we need to
        # concatenate N of them to obtain a greater period coverage.
        hkt_path = "data/ANCDR/BEPICOLOMBO/hkt/"
        hkt_expression = 'mtm_raw_hk_sa_position_????????.tab'
        num_sa_files = 7  # Compare last week

        # Determine files to use to fetch TM data
        sa_files = list_files_from_ftp(hkt_path, hkt_expression)
        sa_files = sa_files[-num_sa_files:]

        # For each file, download it, add data to array, and remove it
        sa_angles = download_tm_data(sa_files, hkt_path, ",", [2], [180.0 / math.pi])
        if not len(sa_angles):
            print("Cannot obtain required TM data, aborting.")
            return None

    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    # Compare with SPICE SA Angles
    error = []
    max_ang_error = 0
    num_gaps = 0

    for sa_angle in sa_angles:

        try:
            et = sa_angle[0]
            hk_sa_angle = sa_angle[1]

            # Determine the rotation matrix to pass from the SA Rotating Frame
            # to the SA Fixed frame
            sadm_rot = spiceypy.pxform(sadm_frame, sadm_ref_frame, et)

            # Covert the SA reference vector in the rotating frame into
            # the SA fixed frame
            sadm_vector = spiceypy.mxv(sadm_rot, ref_vector)
            sadm_angle = np.rad2deg(spiceypy.vsep(ref_vector, sadm_vector))

            # Because vsep always is positive, we are going to get the cross product to
            # determine if is a positive or negative rotation
            sadm_cross_vector = np.cross(ref_vector, sadm_vector)

            # The dot product of the normalised vectors shall be or 1 or -1.
            sadm_angle = np.dot(spiceypy.unorm(ref_cross_vector)[0], spiceypy.unorm(sadm_cross_vector)[0]) * sadm_angle

            ang_error = abs(sadm_angle - hk_sa_angle) * 1000  # mdeg
            max_ang_error = max(max_ang_error, ang_error)

            error.append([et, ang_error])

        except SpiceNOFRAMECONNECT:
            # There is a gap in the CK file, ignore this SA sample.
            num_gaps += 1
            continue

    # Plot error
    if len(error):
        error = np.asarray(error)
        print('Max angular error [mdeg]: ' + str(max_ang_error))

        plot(error[:, 0],
             [error[:, 1]],
             yaxis_name=['Angular error'],
             title=sc + " SA angular error between TM and SPICE",
             format=plot_style,
             yaxis_units='mdeg',
             notebook=notebook)

        return max_ang_error

    else:
        print('Angular error cannot be computed. Found ' + str(num_gaps) + '  of ' + str(len(sa_angles)) + ' samples without data.')
        return None




"""
def sadmCkVsMust(sc, start_time, end_time, plot_style='line', notebook=True):

    # TODO: METHOD NOT FINISHED!!!

    if sc == 'MPO':
        must_sadm_param = 'NCADAF41'  # TODO: Set correct parameter for SADM
        sadm_frame = 'MPO_SA'
        sadm_ref_frame = 'MPO_SA_SADM'
        ref_vector = np.asarray([0, 0, 1])
        mission_phase = 'BEPICRUISE'
    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    et_start = spiceypy.utc2et(start_time)
    et_end = spiceypy.utc2et(end_time)

    start_time = et_to_datetime(et_start).strftime('%Y-%b-%d %H:%M:%S')
    end_time = et_to_datetime(et_end).strftime('%Y-%b-%d %H:%M:%S')

    error = []
    max_ang_error = 0

    tm = WebmustHandler(mission_phase=mission_phase)
    df_0 = tm.get_tm([must_sadm_param], start_time, end_time)

    for row in df_0:
        utc = row[0]
        must_angle = row[1]  # TODO: Must be in degrees
        et = spiceypy.utc2et(utc)

        sadm_rot = spiceypy.pxform(sadm_ref_frame, sadm_frame, et)
        sadm_vector = spiceypy.mxv(sadm_rot, ref_vector)
        sadm_angle = np.rad2deg(spiceypy.vsep(ref_vector, sadm_vector))

        ang_error = abs(sadm_angle - must_angle) * 1000  # mdeg
        max_ang_error = max(max_ang_error, ang_error)

        error.append([et, ang_error])

    # Plot error
    error = np.asarray(error)
    print('Max angular error [mdeg]: ' + str(max_ang_error))

    plot(error[:, 0],
         [error[:, 1]],
         yaxis_name=['Angular error'],
         title=sc + " SA angular error between WebMUST and SPICE",
         format=plot_style,
         yaxis_units='mdeg',
         notebook=notebook)

    return max_ang_error
"""


def cov_ck_obj(mk, object, time_format= 'UTC', global_boundary=False,
               report=False, unload=False):
    """
    Provides time coverage summary for a given object for a list of
    binary CK files provided in a meta-kernel. Several options are
    available. This function is based on the following SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param mk: Meta-kernel to load the computation scenario.
    :type mk: str
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """
    spiceypy.furnsh(mk)
    boundaries_list = []
    et_boundaries_list = []

    object_id = spiceypy.namfrm(object)
    MAXIV = 2000
    ck_count = spiceypy.ktotal('CK') - 1
    WINSIZ = 2 * MAXIV
    MAXOBJ = 10000

    while ck_count >= 0:

        ck_ids = spiceypy.support_types.SPICEINT_CELL(MAXOBJ)
        ck_kernel = spiceypy.kdata(ck_count, 'CK', 155, 155, 155)
        try:
            ck_ids = spiceypy.ckobj(ck=ck_kernel[0], outCell=ck_ids)
        except:
            ck_ids = spiceypy.ckobj(ck=ck_kernel[0])

        for id in ck_ids:
            if id == object_id:
                object_cov = spiceypy.support_types.SPICEDOUBLE_CELL(WINSIZ)
                object_cov = spiceypy.ckcov(ck=ck_kernel[0], idcode=object_id,
                                          needav=False, level='SEGMENT',
                                          tol=0.0, timsys='TDB',
                                          cover=object_cov)

                boundaries = time.cov_int(object_cov=object_cov,
                                          object_id=object_id,
                                          kernel=ck_kernel[0],
                                          global_boundary=global_boundary,
                                          time_format=time_format,
                                          report=report)

                boundaries_list.append(boundaries)

                #
                # We need to have the boundaries in TDB in order to sort out the
                # min and max to obtain the global ones for multiple kernels
                #
                if global_boundary:
                    et_boundaries_list.append(time.cov_int(
                                                      object_cov=object_cov,
                                                      object_id=object_id,
                                                      kernel=ck_kernel[0],
                                                      global_boundary=True,
                                                      time_format='TDB',
                                                      report=False))

        ck_count -= 1

    if global_boundary:
        start_time = min(et_boundaries_list)[0]
        finish_time = max(et_boundaries_list)[1]

        boundaries_list = time.et2cal([start_time, finish_time],
                                      format=time_format)

        if report:

            try:
                body_name = spiceypy.bodc2n(object_id)
            except:
                body_name = spiceypy.frmnam(object_id, 60)

            print("Global Coverage for {} [{}]: {} - {}".format(
                   body_name, time_format, boundaries_list[0],
                boundaries_list[1]))

    if unload:
        spiceypy.unload(mk)

    return boundaries_list


def cov_ck_ker(ck, object, support_ker=list(), time_format='UTC',
               report=False, unload=True):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended. It
       is optional since the kernels could have been already loaded.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format), 'TDB' or 'SPICE'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """
    spiceypy.furnsh(ck)

    if support_ker:

        if isinstance(support_ker, str):
            support_ker = [support_ker]

        for ker in support_ker:
            spiceypy.furnsh(ker)

    object_id = spiceypy.namfrm(object)
    MAXIV = 200000
    WINSIZ = 2 * MAXIV
    MAXOBJ = 100000

    ck_ids = spiceypy.support_types.SPICEINT_CELL(MAXOBJ)
    try:
        ck_ids = spiceypy.ckobj(ck, outCell=ck_ids)
    except:
        ck_ids = spiceypy.ckobj(ck)

    if object_id in ck_ids:

        object_cov = spiceypy.support_types.SPICEDOUBLE_CELL(WINSIZ)
        spiceypy.scard, 0, object_cov
        object_cov = spiceypy.ckcov(ck=ck, idcode=object_id,
                                      needav=False, level='INTERVAL',
                                      tol=0.0, timsys='TDB',
                                      cover=object_cov)

    else:
        #print('{} with ID {} is not present in {}.'.format(object,
         #                                                  object_id, ck))
        if unload:
            spiceypy.unload(ck)
            if support_ker:

                if isinstance(support_ker, str):
                    support_ker = [support_ker]

                for ker in support_ker:
                    spiceypy.unload(ker)

        return False

    if time_format == 'SPICE':
        boundaries = object_cov

    else:
        boundaries = time.cov_int(object_cov=object_cov,
                                  object_id=object_id,
                                  kernel=ck,
                                  time_format=time_format, report=report)

    if unload:
        spiceypy.unload(ck)
        if support_ker:

            if isinstance(support_ker, str):
                support_ker = [support_ker]

            for ker in support_ker:
                spiceypy.unload(ker)

    return boundaries


def time_correlation(sc, ck, plot_style='line', notebook=True):

    # Downloads a telemetry file of a given CK and computes
    # the time difference between the UTC time (1st column)
    # and the clock string (2nd column) in milliseconds

    spiceypy.timdef('SET', 'SYSTEM', 10, 'UTC')

    if sc == 'MPO':
        file = ck.split('/')[-1].replace('\n', '').split('_')[5]
        file = 'mpo_raw_hk_aocs_measured_attitude_' + file + '.tab'
        download_file("data/ANCDR/BEPICOLOMBO/hkt", file)
    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    print('AOCS tab file: ' + file)
    if not os.path.isfile(file):
        print('AOCS tab file cannot be downloaded!')
        return None

    tabfile = open(file)
    times = []
    time_diff = []

    try:
        sc_id = spiceypy.bodn2c(sc)
    except Exception as e:
        print('Spacecraft not found: ' + sc + ", err: " + str(e))
        return None

    for line in tabfile.readlines():
        data = line.replace('\n', '').replace(',', ' ').split()
        utc_et = spiceypy.str2et(data[0].replace('Z', ''))
        scs_et = spiceypy.scs2e(sc_id, data[1])

        times.append(utc_et)
        time_diff.append((utc_et - scs_et) * 1000)

    time_diff = np.abs(np.asarray(time_diff))
    max_time_diff = np.max(time_diff)
    print('Avg time difference [ms]: ' + str(np.mean(time_diff)))
    print('Max time difference [ms]: ' + str(max_time_diff))

    plot(times, time_diff,
         yaxis_name='Time diff (UTC - SCS)',
         title='Time difference between UTC and Clock String in milliseconds',
         format=plot_style,
         yaxis_units='milliseconds',
         notebook=notebook)

    os.remove(file)

    return max_time_diff


def flyby_ca_altitudes(sc, target, spk_expression, num_spk_files, from_date, to_date,
                       distance_flyby, num_samples, plot_style='line', notebook=True, plot_prefix=""):

    spiceypy.timdef('SET', 'SYSTEM', 10, 'TDB')

    target = target.upper()
    target_frame = "IAU_" + target
    start_time = spiceypy.utc2et(from_date)
    stop_time = spiceypy.utc2et(to_date)
    times = np.linspace(start_time, stop_time, num_samples)
    maxwin = 200000

    if sc == 'MPO':
        spk_path = "data/SPICE/BEPICOLOMBO/kernels/spk/"
    else:
        print('Unsupported spacecraft: ' + sc)
        return None

    # Download num_spk_files and find the flyby for each one
    spk_files = list_files_from_ftp(spk_path, spk_expression)
    spk_files = spk_files[-num_spk_files:]

    flybys_spks = []
    flybys_ets = []
    flybys_alts = []
    flybys_alt_list = []

    for spk_file in spk_files:

        # Download spk file
        download_file(spk_path, spk_file)

        if not os.path.isfile(spk_file):
            print('OEM file cannot be downloaded!')
            return None

        # Obtain the flyby data
        spiceypy.furnsh(spk_file)

        cnfine = SPICEDOUBLE_CELL(maxwin)
        result = SPICEDOUBLE_CELL(maxwin)

        spiceypy.scard(0, cnfine)
        spiceypy.wninsd(start_time, stop_time, cnfine)

        spiceypy.gfdist(target=target,
                      abcorr='NONE',
                      obsrvr=sc,
                      relate='<',
                      refval=distance_flyby,
                      step=spiceypy.spd() / 4.,
                      nintvls=maxwin,
                      cnfine=cnfine,
                      adjust=0.0,
                      result=result)

        final = SPICEDOUBLE_CELL(maxwin)

        spiceypy.gfdist(target=target,
                      abcorr='NONE',
                      obsrvr=sc,
                      relate='LOCMIN',
                      refval=distance_flyby,
                      step=spiceypy.spd() / 4.,
                      nintvls=maxwin,
                      cnfine=result,
                      adjust=0.0,
                      result=final)

        number_of_results = spiceypy.wncard(final)

        if number_of_results == 0:
            print('No ' + target + ' flyby found for that period at SPK: ' + spk_file)
            return None

        if number_of_results > 1:
            print('Error: Multiple ' + target + ' flybys found for that period at SPK: ' + spk_file)
            return None

        # Get the flyby closest approach time and distance
        flyby_et = spiceypy.wnfetd(final, 0)[0]
        (state, lt) = spiceypy.spkezr(target, flyby_et, target_frame, 'NONE', sc)
        spoint = spiceypy.sincpt('ELLIPSOID', target, flyby_et, target_frame,
                                 'NONE', sc, target_frame, state[:3])[0]
        flyby_altitude = spiceypy.vnorm(spoint + state[:3])

        flybys_spks.append(spk_file)
        flybys_ets.append(flyby_et)
        flybys_alts.append(flyby_altitude)

        # Get the flyby closest approach distance evolution
        altitudes = []
        for et in times:
            (state, lt) = spiceypy.spkezr(target, et, target_frame, 'NONE', sc)
            spoint = spiceypy.sincpt('ELLIPSOID', target, et, target_frame,
                                     'NONE', sc, target_frame, state[:3])[0]
            altitudes.append(spiceypy.vnorm(spoint + state[:3]))
        flybys_alt_list.append(altitudes)

        # Unload and clean spk file
        spiceypy.unload(spk_file)
        os.remove(spk_file)

    # Reduce the SPK names to only the SPK number to reduce the legend size
    spk_numbers = []
    if sc == 'MPO':
        for spk in flybys_spks:
            spk_number = int(spk.split("_")[3])
            spk_numbers.append(spk_number)

    if len(plot_prefix):
        plot_prefix = plot_prefix + " "

    # Plot Flyby CA altitude vs spk number
    plot(spk_numbers,
         flybys_alts,
         title=plot_prefix + target + ' Flyby CA altitude vs spk number',
         format="scatter",
         xaxis_name='SPK Number',
         yaxis_name=['Altitude'],
         yaxis_units='Km',
         notebook=notebook)

    # Plot Flyby CA time vs spk number
    plot(flybys_ets,
         spk_numbers,
         title=plot_prefix + target + ' Flyby CA spk number vs time',
         format="scatter",
         yaxis_name=['SPK Number'],
         yaxis_units='SPK Number',
         notebook=notebook)

    # Plot Flyby altitude evolution vs time per SPK
    plot(times,
         flybys_alt_list,
         yaxis_name=spk_numbers,
         title=plot_prefix + target + ' Flyby altitude evolution vs time',
         format=plot_style,
         yaxis_units='Km',
         plot_height=400,
         notebook=notebook)

    return np.max(flybys_alts)


def fk_body_ifj2000(mission, body, pck, body_spk, frame_id, report=False,
                    unload=False, file=True):
    """
    Generates a given Solar System Natural Body Inertial frame at J2000. This
    function is based on a FORTRAN subroutine provided by Boris Semenov
    (NAIF/JPL)

    The frame definition would be as follows:

    {Body} Inertial Frame at J2000 ({MISSION}_{BODY}_IF_J2000)

    Definition:

    The {body} Inertial Frame at J2000 is defined as follows:

       -  +Z axis is parallel to {body} rotation axis
          at J2000, pointing toward the North side of the
          invariable plane;

       -  +X axis is aligned with the ascending node of the {Body}
          orbital plane with the {Body} equator plane at J2000;

       -  +Y axis completes the right-handed system;

       -  the origin of this frame is the center of mass of {Body}.

    All vectors are geometric: no aberration corrections are used.


    Remarks:

    This frame is defined as a fixed offset frame using constant vectors
    as the specification method. The fixed offset for these vectors were
    based on the following directions (that also define a two-vector
    frame):

      - +Z axis along Right Ascension (RA) and Declination (DEC) of {Body}
        pole at J2000 epoch in J2000 inertial frame;

      - +X axis along the RA/DEC of {Body} instantaneous orbital plane
        ascending node on {Body} equator at J2000 epoch in J2000
        inertial frame;

    This frame has been defined based on the IAU_{BODY} frame, whose
    evaluation was based on the data included in the loaded PCK file.

    In addition {body_spk} ephemeris have been used to compute the {Body}
    instantaneous orbital plane ascending node on {Body} equator at
    J2000 epoch in J2000 inertial frame.

    :param mission: Name of the mission to use the frame
    :type mission: str
    :param body: Natural body for which the frame is defined
    :type body: str§
    :param pck: Planetary Constants Kernel to be used to extract the Pole information from
    :type pck: str
    :param body_spk: SPK kernels that contain the ephemeris of the Natural body
    :type body_spk: Union[str, list]
    :param frame_id: ID for the new frame. It is recommended to follow the convention recommended by NAIF: -XYYY where X is the ID of the mission S/C and YYY is a number between 900 and 999.
    :type frame_id: str
    :param report: If True prints some intermediate results.
    :type report: bool
    :param unload: If True it will unload the input PCK and SPK.
    :type unload: bool
    :param file: If True it generates the frame definition in a file with the following name: {MISSION}_{BODY}_IF_J2000.tf
    :type file: bool
    :return: Returns the Euler angles to transform the computed frame with J2000. Only if parameter file is False
    :rtype: str
    """
    body = body.upper()
    mission = mission.upper()

    spiceypy.furnsh(pck)

    #
    # This can actually be a list of bodies.
    #
    spiceypy.furnsh(body_spk)

    #
    # Get instantaneous Body state at J2000 and compute instantaneous
    # orbital normal.
    #
    state, lt = spiceypy.spkezr(body, 0.0, 'J2000', 'NONE', 'SUN')
    normal = spiceypy.ucrss(state[0:3:1], state[3:6:1])

    #
    # Get J2000 -> IAU_{BODY} rotation at J2000 and compute Body pole
    # direction in J2000 at J2000.
    #
    mat = spiceypy.pxform('IAU_{}'.format(body), 'J2000', 0.0)
    z = spiceypy.vpack(0.0, 0.0, 1.0)
    pole = spiceypy.mxv(mat, z)

    #
    # Compute direction Body orbit's ascending node on Body equator at
    # J2000 in J2000 and print it and Body pole as RA/DEC in J2000 in
    # degrees
    #
    ascnod = spiceypy.ucrss(pole, normal)
    r, ra, dec = spiceypy.recrad(pole)

    if report:
        print('POLE RA/DEC = {}/{}'.format(ra*spiceypy.dpr(), dec*spiceypy.dpr()))

    r, ra, dec = spiceypy.recrad(ascnod)

    if report:
        print('ASCNOD RA/DEC = {}/{}'.format(ra * spiceypy.dpr(), dec * spiceypy.dpr()))

    #
    # Build two vector from a with POLE as Z and ASNOD as X and print rotation
    # from that frame to J200 as Euler angles.
    #
    mat = spiceypy.twovec(pole, 3, ascnod, 1)
    matxp = spiceypy.xpose(mat)
    r3, r2, r1 = spiceypy.m2eul(matxp, 3, 2, 3)

    if file:
      body_id = spiceypy.bodn2c(body)
      with open('{}_{}_IF_J2000.tf'.format(mission, body), 'w+') as f:

         f.write(r"\begindata")
         f.write('\n \n')
         f.write('    FRAME_{}_{}_IF_J2000   = {}\n'.format(mission, body,
                                                            frame_id))
         f.write("    FRAME_{}_NAME              = '{}_{}_IF_J2000'\n".format(
                 frame_id, mission, body))
         f.write('    FRAME_{}_CLASS             =  4\n'.format(frame_id))
         f.write('    FRAME_{}_CLASS_ID          = {}\n'.format(frame_id,
                                                              frame_id))
         f.write('    FRAME_{}_CENTER            =  {}\n'.format(frame_id,
                                                                 body_id))
         f.write('\n')
         f.write("    TKFRAME_{}_SPEC            = 'ANGLES'\n".format(frame_id))
         f.write("    TKFRAME_{}_RELATIVE        = 'J2000'\n".format(frame_id))
         f.write('    TKFRAME_{}_ANGLES          = (\n'.format(frame_id))
         f.write('                                        {}\n'.format(r3 *
                                                                     spiceypy.dpr()))
         f.write('                                        {}\n'.format(r2 *
                                                                     spiceypy.dpr()))
         f.write('                                        {}\n'.format(r1 *
                                                                     spiceypy.dpr()))
         f.write('                                     )\n')
         f.write('    TKFRAME_{}_AXES            = (\n'.format(frame_id))
         f.write('                                        3,\n')
         f.write('                                        2,\n')
         f.write('                                        3\n')
         f.write('                                     )\n')
         f.write("    TKFRAME_{}_UNITS           = 'DEGREES'\n".format(frame_id))
         f.write('\n')
         f.write(r"\begintext")

    else:
        return '{}_IF->J2000 (3-2-3): {} - {} - {}'.format(body,
            r3 * spiceypy.dpr(),
            r2 * spiceypy.dpr(),
            r1 * spiceypy.dpr())

    if unload:
        spiceypy.unload(pck)
        spiceypy.unload(body_spk)

    return


def eul_angle_report(et_list, eul_ck1, eul_ck2, eul_num, tolerance, name=''):

    eul_error = list(numpy.degrees(abs(numpy.array(eul_ck1) - numpy.array(eul_ck2))))

    count = 0
    interval_bool = False
    eul_tol_list = []

    with open('euler_angle_{}_{}_report.txt'.format(eul_num, name), 'w+') as f:
        f.write('EULER ANGLE {} REPORT \n'.format(eul_num))
        f.write('==================== \n')

        for element in eul_error:

            if element >= tolerance:
                if interval_bool:
                    eul_tol_list.append(element)
                else:
                    interval_bool = True
                    eul_tol_list.append(element)
                    utc_start = spiceypy.et2utc(et_list[count], 'ISOC', 2)

            else:
                if interval_bool:
                    utc_finish = spiceypy.et2utc(et_list[count], 'ISOC', 2)

                    f.write('TOLERANCE of ' + str(tolerance) + ' DEG exceeded from ' + utc_start + ' until ' +
                          utc_finish + ' with an average angle of ' + str(numpy.mean(eul_tol_list)) + ' DEG \n')

                interval_bool = False

            count += 1

        f.write('\nMAX Error:  {} DEG\n'.format(str(max(eul_error))))
        f.write('MIN Error:   {} DEG\n'.format(str(min(eul_error))))
        f.write('MEAN Error: {} DEG\n'.format(str(numpy.mean(eul_error))))

    return


def attitude_error_report(et_list, ang_ck1, ang_ck2, tolerance, name=''):

    ang_error = list(abs(numpy.array(ang_ck1) - numpy.array(ang_ck2)))

    count = 0
    interval_bool = False
    ang_tol_list = []

    with open('attitude_error_{}_report.txt'.format(name), 'w+') as f:
        f.write('ATTITUDE ERROR REPORT \n')
        f.write('==================== \n')

        for element in ang_error:

            if element >= tolerance:
                if interval_bool:
                    ang_tol_list.append(element)
                else:
                    interval_bool = True
                    ang_tol_list.append(element)
                    utc_start = spiceypy.et2utc(et_list[count], 'ISOC', 2)

            else:
                if interval_bool:
                    utc_finish = spiceypy.et2utc(et_list[count], 'ISOC', 2)

                    f.write('TOLERANCE of ' + str(tolerance) + ' DEG exceeded from ' + utc_start + ' until ' +
                          utc_finish + ' with an average angle of ' + str(numpy.mean(ang_tol_list)) + ' DEG \n')

                interval_bool = False

            count += 1

        f.write('\nMAX Error:  {} ARCSECONDS\n'.format(str(max(ang_error))))
        f.write('MIN Error:   {} ARCSECONDS\n'.format(str(min(ang_error))))
        f.write('MEAN Error: {} ARCSECONDS\n'.format(str(numpy.mean(ang_error))))

    return


def state_report(et_list, pos_spk1, pos_spk2, vel_spk1, vel_spk2, pos_tolerance,
                 vel_tolerance, name=''):

    pos_error = list(abs(numpy.array(pos_spk1) - numpy.array(pos_spk2)))
    vel_error = list(abs(numpy.array(vel_spk1) - numpy.array(vel_spk2)))

    count = 0
    interval_bool = False
    pos_tol_list = []

    with open('state_{}_report.txt'.format(name), 'w+') as f:
        f.write('STATE REPORT \n')
        f.write('============ \n')

        for element in pos_error:

            if element >= pos_tolerance:
                if interval_bool:
                    pos_tol_list.append(element)
                else:
                    interval_bool = True
                    pos_tol_list.append(element)
                    utc_start = spiceypy.et2utc(et_list[count], 'ISOC', 2)

            else:
                if interval_bool:
                    utc_finish = spiceypy.et2utc(et_list[count], 'ISOC', 2)

                    f.write('TOLERANCE of ' + str(pos_tolerance) + ' KM exceeded from ' + utc_start + ' until ' +
                          utc_finish + ' with an average distance of ' + str(numpy.mean(pos_tol_list)) + ' KM \n')

                interval_bool = False

            count += 1

        count = 0
        interval_bool = False
        vel_tol_list = []

        for element in vel_error:

            if element >= vel_tolerance:
                if interval_bool:
                    vel_tol_list.append(element)
                else:
                    interval_bool = True
                    vel_tol_list.append(element)
                    utc_start = spiceypy.et2utc(et_list[count], 'ISOC', 2)

            else:
                if interval_bool:
                    utc_finish = spiceypy.et2utc(et_list[count], 'ISOC', 2)

                    f.write('TOLERANCE of ' + str(vel_tolerance) + ' KM/S exceeded from ' + utc_start + ' until ' +
                          utc_finish + ' with an average velocity of ' + str(numpy.mean(vel_tol_list)) + ' KM/S \n')

            count += 1

        f.write('\nMAX Error:  {} KM\n'.format(str(max(pos_error))))
        f.write('MIN Error:   {} KM\n'.format(str(min(pos_error))))
        f.write('MEAN Error: {} KM\n'.format(str(numpy.mean(pos_error))))

        f.write('\nMAX Error:  {} KM/S\n'.format(str(max(vel_error))))
        f.write('MIN Error:   {} KM/S\n'.format(str(min(vel_error))))
        f.write('MEAN Error: {} KM/S\n'.format(str(numpy.mean(vel_error))))

    return


def ckdiff_euler(mk, ck1, ck2, spacecraft_frame, target_frame, resolution, tolerance,
           utc_start='', utc_finish='', plot_style='line', report=True,
           notebook=False):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """

    # Compute time windows
    spiceypy.furnsh(mk)

    windows_ck1 = cov_ck_ker(ck1, object=spacecraft_frame, time_format='SPICE')
    spiceypy.unload(ck1)

    windows_ck2 = cov_ck_ker(ck2, object=spacecraft_frame, time_format='SPICE')
    spiceypy.unload(ck2)

    windows_intersected = spiceypy.wnintd(windows_ck1, windows_ck2)

    number_of_intervals = list(range(spiceypy.wncard(windows_intersected)))

    et_boundaries_list = []
    for element in number_of_intervals:
        et_boundaries = spiceypy.wnfetd(windows_intersected, element)
        et_boundaries_list.append(et_boundaries[0])
        et_boundaries_list.append(et_boundaries[1])

    start = True
    for et_start, et_finish in zip(et_boundaries_list[0::2], et_boundaries_list[1::2]):

        if start:
            et_list = numpy.arange(et_start, et_finish, resolution)
            start = False

        et_list = numpy.append(et_list, numpy.arange(et_start, et_finish, resolution))

    if utc_start:
        et_start = spiceypy.utc2et(utc_start)

    if utc_finish:
        et_finish = spiceypy.utc2et(utc_finish)
        et_list = numpy.arange(et_start, et_finish, resolution)

    # Process CK1
    spiceypy.furnsh(ck1)

    eul1_ck1 = []
    eul2_ck1 = []
    eul3_ck1 = []
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame,  target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck1.append(math.degrees(euler[0]))
        eul2_ck1.append(math.degrees(euler[1]))
        eul3_ck1.append(math.degrees(euler[2]))

    spiceypy.unload(ck1)

    # Process CK2
    spiceypy.furnsh(ck2)

    eul1_ck2 = []
    eul2_ck2 = []
    eul3_ck2 = []
    for et in et_list:
        rot_mat = spiceypy.pxform(spacecraft_frame, target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck2.append(math.degrees(euler[0]))
        eul2_ck2.append(math.degrees(euler[1]))
        eul3_ck2.append(math.degrees(euler[2]))

    spiceypy.unload(ck2)

    # Plot angles
    ck1_filename = ck1.split('/')[-1].split('.')[0]
    ck2_filename = ck2.split('/')[-1].split('.')[0]

    eul1_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul2_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul3_name = '{}_{}'.format(ck1_filename, ck2_filename)

    plot(et_list, [eul1_ck1, eul1_ck2], yaxis_name=['Euler Angle 1 CK1',
                                                   'Euler Angle 1 CK2'],
                                                    title='Euler Angle 1 {}'.format(eul1_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    plot(et_list, [eul2_ck1, eul2_ck2], yaxis_name=['Euler Angle 2 CK1',
                                                   'Euler Angle 2 CK2'],
                                                    title='Euler Angle 2 {}'.format(eul2_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    plot(et_list, [eul3_ck1, eul3_ck2], yaxis_name=['Euler Angle 3 CK1',
                                                   'Euler Angle 3 CK2'],
                                                    title='Euler Angle 3 {}'.format(eul3_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    # Generate reports
    if report:
        eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=eul1_name)
        eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=eul2_name)
        eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=eul3_name)

    return


def ckdiff(ck1, ck2, spacecraft_frame, target_frame, resolution, tolerance,
           utc_start='', utc_finish='', mk='', output='boresight', boresight = [0,0,1],
           plot_style='line', report=False, notebook=False):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """

    if mk:
        spiceypy.furnsh(mk)

    windows_ck1 = cov_ck_ker(ck1, object=spacecraft_frame, time_format='SPICE')
    spiceypy.unload(ck1)

    windows_ck2 = cov_ck_ker(ck2, object=spacecraft_frame, time_format='SPICE')
    spiceypy.unload(ck2)

    windows_intersected = spiceypy.wnintd(windows_ck1, windows_ck2)

    number_of_intervals = list(range(spiceypy.wncard(windows_intersected)))

    et_boundaries_list = []
    for element in number_of_intervals:
        et_boundaries = spiceypy.wnfetd(windows_intersected, element)
        et_boundaries_list.append(et_boundaries[0])
        et_boundaries_list.append(et_boundaries[1])

    start = True
    for et_start, et_finish in zip(et_boundaries_list[0::2], et_boundaries_list[1::2]):

        if start:
            et_list = numpy.arange(et_start, et_finish, resolution)
            start = False

        et_list = numpy.append(et_list, numpy.arange(et_start, et_finish, resolution))

    if utc_start:
        et_start = spiceypy.utc2et(utc_start)

    if utc_finish:
        et_finish = spiceypy.utc2et(utc_finish)
        et_list = numpy.arange(et_start, et_finish, resolution)

    spiceypy.furnsh(ck1)
    eul1_ck1 = []
    eul2_ck1 = []
    eul3_ck1 = []
    bsight_ck1 = []
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame, target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck1.append(math.degrees(euler[0]))
        eul2_ck1.append(math.degrees(euler[1]))
        eul3_ck1.append(math.degrees(euler[2]))

        bsight = spiceypy.mxv(rot_mat, boresight)
        bsight_ang = spiceypy.vsep(bsight, boresight)
        bsight_ck1.append(bsight_ang*spiceypy.dpr())
    spiceypy.unload(ck1)

    spiceypy.furnsh(ck2)
    eul1_ck2 = []
    eul2_ck2 = []
    eul3_ck2 = []
    bsight_ck2 = []
    for et in et_list:
        rot_mat = spiceypy.pxform(spacecraft_frame, target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck2.append(math.degrees(euler[0]))
        eul2_ck2.append(math.degrees(euler[1]))
        eul3_ck2.append(math.degrees(euler[2]))

        bsight = spiceypy.mxv(rot_mat, boresight)
        bsight_ang = spiceypy.vsep(bsight, boresight)
        bsight_ck2.append(bsight_ang*spiceypy.dpr())
    spiceypy.unload(ck2)

    ck1_filename = ck1.split('/')[-1].split('.')[0]
    ck2_filename = ck2.split('/')[-1].split('.')[0]

    title_name = '{}_{}'.format(ck1_filename, ck2_filename)

    if output == 'euler_angles':

        plot(et_list, [eul1_ck1,eul1_ck2], yaxis_name=['Degrees', 'Degrees'],
                                                        title='Euler Angle 1 {}'.format(title_name),
                                                        format=plot_style,
                                                        notebook=notebook)

        plot(et_list, [eul2_ck1,eul2_ck2], yaxis_name=['Degrees', 'Degrees'],
                                                        title='Euler Angle 2 {}'.format(title_name),
                                                        format=plot_style,
                                                        notebook=notebook)

        plot(et_list, [eul3_ck1,eul3_ck2], yaxis_name=['Degrees', 'Degrees'],
                                                        title='Euler Angle 3 {}'.format(title_name),
                                                        format=plot_style,
                                                        notebook=notebook)
    else:
        plot(et_list, [bsight_ck1], yaxis_name=['Degrees', 'Degrees'],
             title='+Z Axis Angle Difference {}'.format(title_name),
             format=plot_style,
             notebook=notebook)

    if report:
        eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=title_name)
        eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=title_name)
        eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=title_name)

    return


def get_euler_boresights_angles(ck, et_list, spacecraft_frame,
                               target_frame, boresight):
    spiceypy.furnsh(ck)

    eul1 = []
    eul2 = []
    eul3 = []
    bsights = []
    angles = []

    for et in et_list:
        rot_mat = spiceypy.pxform(spacecraft_frame, target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1.append(math.degrees(euler[0]))
        eul2.append(math.degrees(euler[1]))
        eul3.append(math.degrees(euler[2]))

        bsight = spiceypy.mxv(rot_mat, boresight)
        bsight_ang = spiceypy.vsep(bsight, boresight)
        bsights.append(spiceypy.convrt(bsight_ang, 'RADIANS', 'ARCSECONDS'))

        (rot_axis, rot_angle) = spiceypy.raxisa(rot_mat)
        angles.append(spiceypy.convrt(rot_angle, 'RADIANS', 'ARCSECONDS'))

    spiceypy.unload(ck)

    return eul1, eul2, eul3, bsights, angles


def ckdiff_error(ck1, ck2, spacecraft_frame, target_frame, resolution, tolerance,
                 mk='', utc_start='', utc_finish='', output='',
                 boresight=[0,0,1], plot_style='line', report=False,
                 notebook=False):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """
    if mk:
        spiceypy.furnsh(mk)

    try:
        windows_ck1 = cov_ck_ker(ck1, object=spacecraft_frame, time_format='SPICE')
        spiceypy.unload(ck1)

        windows_ck2 = cov_ck_ker(ck2, object=spacecraft_frame, time_format='SPICE')
        spiceypy.unload(ck2)
    except:
        print('WARNING: No Time Window could be determined')
        return None

    windows_intersected = spiceypy.wnintd(windows_ck1, windows_ck2)

    number_of_intervals = list(range(spiceypy.wncard(windows_intersected)))

    if not len(number_of_intervals):
        print('WARNING: No Time Windows intersected')
        return None

    et_boundaries_list = []
    for element in number_of_intervals:
        et_boundaries = spiceypy.wnfetd(windows_intersected, element)
        et_boundaries_list.append(et_boundaries[0])
        et_boundaries_list.append(et_boundaries[1])

    start = True
    for et_start, et_finish in zip(et_boundaries_list[0::2], et_boundaries_list[1::2]):

        if start:
            et_list = numpy.arange(et_start, et_finish, resolution)
            start = False

        et_list = numpy.append(et_list, numpy.arange(et_start, et_finish, resolution))

    if utc_start:
        et_start = spiceypy.utc2et(utc_start)

    if utc_finish:
        et_finish = spiceypy.utc2et(utc_finish)
        et_list = numpy.arange(et_start, et_finish, resolution)

    if not len(et_list):
        print('WARNING: No valid time period')
        return None

    eul1_ck1, eul2_ck1, eul3_ck1, bsight_ck1, angle_ck1 = get_euler_boresights_angles(ck1, et_list, spacecraft_frame,
                                                                                      target_frame, boresight)

    eul1_ck2, eul2_ck2, eul3_ck2, bsight_ck2, angle_ck2 = get_euler_boresights_angles(ck2, et_list, spacecraft_frame,
                                                                                      target_frame, boresight)

    angle_diff = [abs(i - j) for i, j in zip(angle_ck1, angle_ck2)]

    if output == 'euler_angles':

        eul1_diff = [i - j for i, j in zip(eul1_ck1, eul1_ck2)]
        eul2_diff = [i - j for i, j in zip(eul2_ck1, eul2_ck2)]
        eul3_diff = [i - j for i, j in zip(eul3_ck1, eul3_ck2)]

        plot(et_list, [eul1_diff, eul2_diff, eul3_diff],
             yaxis_name=['Degrees', 'Degrees', 'Degrees'],
             title='Euler Angle Differences',
             format=plot_style, yaxis_units='deg',
             notebook=notebook)

    elif output == 'boresight':

        bsight_diff = [np.abs(i - j) for i, j in zip(bsight_ck1, bsight_ck2)]

        plot(et_list, bsight_diff,
             yaxis_name='',
             title='Boresight Angle Difference',
             format=plot_style, yaxis_units='arcsec',
             notebook=notebook)

    # Attitude Error
    else:

        plot(et_list, angle_diff,
             yaxis_name='ang_diff',
             title='Attitude Error',
             format=plot_style, yaxis_units='arcsec',
             notebook=notebook)

    if report:

        ck1_filename = ck1.split('/')[-1].split('.')[0]
        ck2_filename = ck2.split('/')[-1].split('.')[0]

        bsight_name = '{}_{}'.format(ck1_filename, ck2_filename)

        attitude_error_report(et_list, bsight_ck1, bsight_ck2, tolerance, name=bsight_name)

        if output == 'euler_angles':

            eul1_name = '{}_{}'.format(ck1_filename, ck2_filename)
            eul2_name = '{}_{}'.format(ck1_filename, ck2_filename)
            eul3_name = '{}_{}'.format(ck1_filename, ck2_filename)

            eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=eul1_name)
            eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=eul2_name)
            eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=eul3_name)

        elif output == 'rotaxis':

            rotaxis_name = '{}_{}'.format(ck1_filename, ck2_filename)

            attitude_error_report(et_list, angle_ck1, angle_ck2, tolerance, name=rotaxis_name)

    if mk:
        spiceypy.unload(mk)

    return np.max(angle_diff)


def ckplot(ck1, spacecraft_frame, target_frame, resolution,
           mk = '', utc_start='', utc_finish='', notebook=False,
           plot_style='circle'):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """

    if mk:
        spiceypy.furnsh(mk)
    spiceypy.furnsh(ck1)

    et_boundaries_list = cov_ck_ker(ck1, support_ker=mk, object=spacecraft_frame,
                                        time_format='TDB')


    start = True
    for et_start, et_finish in zip(et_boundaries_list[0::2], et_boundaries_list[1::2]):

        if start:
            et_list = numpy.arange(et_start, et_finish, resolution)
            start = False

        et_list = numpy.append(et_list, numpy.arange(et_start, et_finish, resolution))


    # TODO: if we want to really use start and end times and intersect it with the available intervals we need to develop this
    if utc_start:
        et_start = spiceypy.utc2et(utc_start)

    if utc_finish:
        et_finish = spiceypy.utc2et(utc_finish)
        et_list = numpy.arange(et_start, et_finish, resolution)


    eul1 = []
    eul2 = []
    eul3 = []
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame,  target_frame,et)
        euler = spiceypy.m2eul(rot_mat, 1, 2, 3)
        eul1.append(math.degrees(euler[0]))
        eul2.append(math.degrees(euler[1]))
        eul3.append(math.degrees(euler[2]))

    spiceypy.unload(ck1)

    plot(et_list, [eul1,eul2,eul3],
         yaxis_name=['Euler Angle 1', 'Euler Angle 2', 'Euler Angle 3'],
         title='Euler Angles for {}'.format(ck1.split('/')[-1]), notebook=notebook, format=plot_style)

    return


def spkdiff(mk, spk1, spk2, spacecraft, target, resolution, pos_tolerance,
            vel_tolerance, target_frame='', utc_start='', utc_finish='',
            plot_style='line', report=True):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/spiceypy/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least
       it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel
       (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
    :param object: Ephemeris Object to obtain the coverage from.
    :type object: str
    :param time_format: Output time format; it can be 'UTC', 'CAL' (for TDB
       in calendar format) or 'TDB'. Default is 'TDB'.
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the
       coverage windows or only the absolute start and finish coverage times.
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen.
    :type report: bool
    :param unload: If True it will unload the input meta-kernel.
    :type unload: bool
    :return: Returns a list with the coverage intervals.
    :rtype: list
    """

    if not target_frame:
        target_frame = 'IAU_{}'.format(target.upper())

    spiceypy.furnsh(mk)

    windows_spk1 = cov_spk_ker(spk1, object=spacecraft, time_format='SPICE')
    spiceypy.unload(spk1)

    windows_spk2 = cov_spk_ker(spk2, object=spacecraft, time_format='SPICE')
    spiceypy.unload(spk2)

    windows_intersected = spiceypy.wnintd(windows_spk1, windows_spk2)

    number_of_intervals = list(range(spiceypy.wncard(windows_intersected)))

    et_boundaries_list = []
    for element in number_of_intervals:
        et_boundaries = spiceypy.wnfetd(windows_intersected, element)
        et_boundaries_list.append(et_boundaries[0])
        et_boundaries_list.append(et_boundaries[1])

    start = True
    for et_start, et_finish in zip(et_boundaries_list[0::2], et_boundaries_list[1::2]):

        if start:
            et_list = np.arange(et_start, et_finish, resolution)
            start = False

        et_list = numpy.append(et_list, numpy.arange(et_start, et_finish, resolution))

    if utc_start:
        et_start = spiceypy.utc2et(utc_start)

    if utc_finish:
        et_finish = spiceypy.utc2et(utc_finish)
        et_list = numpy.arange(et_start, et_finish, resolution)

    spiceypy.furnsh(spk1)

    state_spk1 = []
    state_spk2 = []
    pos_spk1 = []
    pos_spk2 = []
    vel_spk1 = []
    vel_spk2 = []

    for et in et_list:
        state = spiceypy.spkezr(target, et, target_frame, 'NONE', spacecraft)[0]

        state_spk1.append(state)
        pos_spk1.append(np.sqrt(state[0]*state[0] +
                                state[1]*state[1] +
                                state[2]*state[2]))
        vel_spk1.append(np.sqrt(state[3]*state[3] +
                                state[4]*state[4] +
                                state[5]*state[5]))

    spiceypy.unload(spk1)
    spiceypy.furnsh(spk2)
    for et in et_list:
        state = spiceypy.spkezr(target, et, target_frame, 'NONE', spacecraft)[0]

        state_spk2.append(state)
        pos_spk2.append(np.sqrt(state[0]*state[0] +
                                state[1]*state[1] +
                                state[2]*state[2]))
        vel_spk2.append(np.sqrt(state[3]*state[3] +
                                state[4]*state[4] +
                                state[5]*state[5]))


    plot(et_list, [pos_spk1, pos_spk2], yaxis_name=['Position SPK1',
                                                    'Position SPK2'],
         title='Position of {} w.r.t {} ({})'.format(spacecraft, target, target_frame),
         format=plot_style)

    spiceypy.unload(spk2)
    if report:

        spk1_filename = spk1.split('/')[-1].split('.')[0]
        spk2_filename = spk2.split('/')[-1].split('.')[0]

        state_report(et_list, pos_spk1, pos_spk2, vel_spk1, vel_spk2, pos_tolerance, vel_tolerance,
                     name='{}_{}'.format(spk1_filename, spk2_filename))


    return


def pck_body_placeholder(bodies):
    """

    :param bodies:
    :type bodies:
    :return:
    :rtype:
    """

    with open('update_to_pck.tpc', 'w+') as f:

        pl_id = 517
        for body in bodies:

            #
            # Get body NAIF ID.
            #
            try:
                id =  spiceypy.bodn2c(str(body.upper))
            except:
                id = pl_id
                pl_id += 1

            f.write('       {0}  {1}        1       1       1       -    Placeholder radii\n'.format(id, body[:1].upper() + body[1:].lower()))

        f.write('\n\n')
        pl_id = 517
        for body in bodies:

            #
            # Get body NAIF ID.
            #
            try:
                id =  spiceypy.bodn2c(str(body.upper))
            except:
                id = pl_id
                pl_id += 1

            f.write('BODY{}_RADII = (1        1       1   )\n'.format(id))

        f.write('\n\n')
        pl_id = 517
        for body in bodies:

            #
            # Get body NAIF ID.
            #
            try:
                id =  spiceypy.bodn2c(str(body.upper))
            except:
                id = pl_id
                pl_id += 1

            f.write("        FRAME_IAU_{0} = {1}\n".format(body.upper(), id))
            f.write("        FRAME_{0}_NAME = 'IAU_{1}'\n".format(id, body.upper()))
            f.write("        FRAME_{}_CLASS = 2\n".format(id))
            f.write("        FRAME_{0}_CLASS_ID = {0}\n".format(id))
            f.write("        FRAME_{0}_CENTER = {0}\n".format(id))
            f.write("        BODY{}_POLE_RA       = (    0.        0.         0.  )\n".format(id))
            f.write("        BODY{}_POLE_DEC      = (   90.        0.         0.  )\n".format(id))
            f.write("        BODY{}_PM            = (  -90.        0.         0.  )\n".format(id))
            f.write("        BODY{}_LONG_AXIS     = (    0.                       )\n\n".format(id))

    return


def read_ik_with_sectors(sensor_name):

    #
    # Since all IK variable names contain NAIF ID of the instrument,
    # the input sensor acronym, NNN, needs to be expanded into its
    # full name, ROS_RPC_NNN, which then can be used to find the
    # sensor's NAIF ID code.
    #
    sensnm = sensor_name

    secsiz = 0
    secsis = 0

    try:
        sensid = spiceypy.bodn2c(sensnm)
    except:
        print('Cannot determine NAIF ID for {}'.format(sensnm))
        return sensnm, 0, 0, secsiz, secsis, '', []

    #
    # No IK routines can be used to retrieve loaded data.  First,
    # retrieve the number of sectors provided in the
    # INS-NNNNNN_NUMBER_OF_SECTORS keyword (here -NNNNNN is the NAIF ID
    # of the sensor.)
    #
    ikkwd = 'INS#_NUMBER_OF_SECTORS'
    ikkwd = spiceypy.repmi(ikkwd, "#", sensid)

    try:
        secnum = spiceypy.gipool(ikkwd, 0, 2)
    except:
        print('Loaded IK does not contain {}.'.format(ikkwd))
        return sensnm, sensid, 0, secsiz, secsis, '', []

    #
    # Second, retrieve the sector size provided in the
    # INS-NNNNNN_SECTOR_SIZE or INS-NNNNNN_SECTOR_SIZES keyword.
    #
    ikkwd = 'INS#_SECTOR_SIZES'
    ikkwd = spiceypy.repmi(ikkwd, '#', sensid)

    try:
        secsis = spiceypy.gdpool(ikkwd, 0, 2)

        #
        # We need to search for INS-NNNNNN_SECTOR_SIZE in the second place
        # for it would also be found by INS-NNNNNN_SECTOR_SIZES
        #
    except:

        ikkwd = 'INS#_SECTOR_SIZE'
        ikkwd = spiceypy.repmi(ikkwd, '#', sensid)

        try:
            room = int(secnum[0]*secnum[1]*2)
            secsiz = spiceypy.gdpool(ikkwd, 0, room)
        except:
            print('Loaded IK does not contain {}.'.format(ikkwd))
            return sensnm, sensid, secnum, secsiz, secsis, '', []

    #
    # Third, retrieve the frame in which sector view direction are
    # defined. It is provided in the INS-NNNNNN_FRAME keyword.
    #
    ikkwd = 'INS#_FRAME'
    ikkwd = spiceypy.repmi(ikkwd, '#', sensid)

    try:
        secfrm = spiceypy.gcpool(ikkwd, 0, 1)
    except:
        print('Loaded IK does not contain {}.'.format(ikkwd))
        return sensnm, sensid, secnum, secsiz, secsis, secfrm, []

    #
    # Last, retrieve the sector view directions provided in the
    # INS-NNNNNN_SECTOR_DIRECTIONS keyword.
    #
    ikkwd = 'INS#_SECTOR_DIRECTIONS'
    ikkwd = spiceypy.repmi(ikkwd, '#', sensid)

    try:
        room = int(secnum[0]*secnum[1]*3)
        secdir = spiceypy.gdpool(ikkwd, 0, room)


        #
        # Re-arrange the secdir list into a list of lists in which each
        # individual list is a sector direction vector
        #
        secdir_list = []
        secdir_line = []
        count = 0
        for element in secdir:  # Start counting from 1
            secdir_line.append(element)
            count += 1
            if count % 3 == 0:
                secdir_list.append(secdir_line)
                secdir_line = []
                count = 0
        secdir = secdir_list

    except:
        print('Loaded IK does not contain {}.'.format(ikkwd))
        return sensnm, sensid, secnum, secsiz, secsis, secfrm, []

    return sensnm, sensid, secnum, secsiz, secsis, secfrm, secdir


#def mex_tgo_occultations(interval, refval):
#    (out, radii) = spiceypy.bodvrd('MARS', 'RADII', 3)
#
#    # Compute flattening coefficient.
#    re = radii[0]
#    rp = radii[2]
#    f = (re - rp) / re
#
#    a = re
#    b = radii[1]
#    c = rp
#
#    MAXIVL = 10000
#    MAXWIN = 2 * MAXIVL
#    TDBFMT = 'YYYY MON DD HR:MN:SC.### (TDB) ::TDB'
#
#    # Initialize the "confinement" window with the interval
#    # over which we'll conduct the search.
#    cnfine = stypes.SPICEDOUBLE_CELL(2)
#    spiceypy.wninsd(interval.start, interval.finish, cnfine)
#
#    #
#    # In the call below, the maximum number of window
#    # intervals gfposc can store internally is set to MAXIVL.
#    # We set the cell size to MAXWIN to achieve this.
#    #
#    riswin = stypes.SPICEDOUBLE_CELL(MAXWIN)
#
#    #
#    # Now search for the time period, within our confinement
#    # window, during which the apparent target has elevation
#    # at least equal to the elevation limit.
#    #
#    #   VARIABLE        I/O  DESCRIPTION
#    #   --------------- ---  -------------------------------------------------
#    #   SPICE_GF_CNVTOL  P   Convergence tolerance.
#    #   occtyp           I   Type of occultation.
#    #   front            I   Name of body occulting the other.
#    #   fshape           I   Type of shape model used for front body.
#    #   fframe           I   Body-fixed, body-centered frame for front body.
#    #   back             I   Name of body occulted by the other.
#    #   bshape           I   Type of shape model used for back body.
#    #   bframe           I   Body-fixed, body-centered frame for back body.
#    #   abcorr           I   Aberration correction flag.
#    #   obsrvr           I   Name of the observing body.
#    #   step             I   Step size in seconds for finding occultation
#    #                        events.
#    #   cnfine          I-O  SPICE window to which the search is restricted.
#    #   result           O   SPICE window containing results.
#    #
#    spiceypy.gfoclt('ANY', 'MARS', 'ELLIPSOID', 'IAU_MARS', 'MEX',
#                    'POINT', '', 'NONE', 'TGO', 60, cnfine, riswin)
#
#    #
#    # Now we perform another search to constrain the number of occultations by a
#    # distance criteria
#    #
#    cnfine = riswin
#
#    riswin = stypes.SPICEDOUBLE_CELL(MAXWIN)
#
#    #
#    # We're not using the adjustment feature, so
#    # we set `adjust' to zero.
#    #
#    adjust = 0.0
#
#    #
#    # We use a step size of 1 hour
#    #
#    step = 60 * 60
#
#    # nintvls  =  2*n  +  ( m / step )
#    #
#    # where
#    #
#    #    n     is the number of intervals in the confinement
#    #          window
#    #
#    #    m     is the measure of the confinement window, in
#    #          units of seconds
#    #
#    #    step  is the search step size in seconds
#    #
#    ndays = 100
#    nintvls = int(2 * 1 + (ndays * 24 * 60 * 60 / step))
#
#    #
#    # Now search for the time period, within our confinement
#    # window, during which the apparent target has elevation
#    # at least equal to the elevation limit.
#    #
#    #   VARIABLE         I/O  DESCRIPTION
#    #   ---------------  ---  ------------------------------------------------
#    #   SPICE_GF_CNVTOL   P   Convergence tolerance
#    #   target            I   Name of the target body.
#    #   abcorr            I   Aberration correction flag.
#    #   obsrvr            I   Name of the observing body.
#    #   relate            I   Relational operator.
#    #   refval            I   Reference value.
#    #   adjust            I   Adjustment value for absolute extrema searches.
#    #   step              I   Step size used for locating extrema and roots.
#    #   nintvls           I   Workspace window interval count.
#    #
#    #   cnfine           I-O  SPICE window to which the search is confined.
#    #   result            O   SPICE window containing results.
#    #
#    spiceypy.gfdist('MEX', 'NONE', 'TGO', '<', refval, adjust, step, nintvls,
#                    cnfine, riswin)
#
#    #
#    # The function wncard returns the number of intervals
#    # in a SPICE window.
#    #
#    winsiz = spiceypy.wncard(riswin)
#
#    lat_mid_list = []
#    lon_mid_list = []
#    dist_mid_list = []
#
#    lat_list = []
#    lon_list = []
#    dist_list = []
#
#    x, y, z = [], [], []
#
#    if winsiz == 0:
#        print('No events were found.')
#
#    else:
#
#        #
#        # Display the visibility time periods.
#        #
#        print(
#                'Occultation times of {0:s} as seen from {1:s} when the distance is '
#                'less than {2:f} km:\n'.format('MEX', 'TGO', refval))
#
#        for i in range(winsiz):
#            #
#            # Fetch the start and stop times of
#            # the ith interval from the search result
#            # window riswin.
#            #
#            [intbeg, intend] = spiceypy.wnfetd(riswin, i)
#
#            #
#            # Convert the rise time to a TDB calendar string.
#            #
#            timstr = spiceypy.timout(intbeg, TDBFMT)
#            et_rise = intbeg
#
#            #
#            # Write the string to standard output.
#            #
#            # if i == 0:
#            #
#            #    print('Occultation start time:'
#            #          '  {:s}'.format(timstr))
#            # else:
#            #
#            #    print('Occultation start time:'
#            #          '  {:s}'.format(timstr))
#            #
#            #
#            # Convert the set time to a TDB calendar string.
#            #
#            timstr = spiceypy.timout(intend, TDBFMT)
#            et_set = intend
#
#            #
#            # Write the string to standard output.
#            #
#            # if i == (winsiz - 1):
#            #
#            #    print('Occultation or window stop time: '
#            #          '  {:s}'.format(timstr))
#            # else:
#            #
#            #    print('Occultation stop time: '
#            #          '  {:s}'.format(timstr))
#            #
#            # print(' ')
#
#            #
#            # Generate a Time Window with the rise and set times
#            #
#            utc_rise = spiceypy.et2utc(et_rise, 'ISOC', 3)
#            utc_set = spiceypy.et2utc(et_set, 'ISOC', 3)
#
#            time_window = spiops.TimeWindow(utc_rise, utc_set, resolution=1)
#
#            interval = time_window.window
#            num = 0
#            for et in interval:
#
#                num += 1
#
#                (linept, lt) = spiceypy.spkpos('MARS', et, 'IAU_MARS', 'NONE',
#                                               'TGO')
#                (linedr, lt) = spiceypy.spkpos('MEX', et, 'IAU_MARS', 'NONE',
#                                               'TGO')
#
#                #
#                #  Variable  I/O  Description
#                #  --------  ---  --------------------------------------------------
#                #  a          I   Length of ellipsoid's semi-axis in the x direction
#                #  b          I   Length of ellipsoid's semi-axis in the y direction
#                #  c          I   Length of ellipsoid's semi-axis in the z direction
#                #  linept     I   Point on line
#                #  linedr     I   Direction vector of line
#                #  pnear      O   Nearest point on ellipsoid to line
#                #  dist       O   Distance of ellipsoid from line
#                #
#                (pnear, dist) = spiceypy.npedln(a, b, c, linept, linedr)
#
#                (lon, lat, alt) = spiceypy.recpgr('MARS', pnear, re, f)
#
#                lon = spiceypy.dpr() * lon
#                lat = spiceypy.dpr() * lat
#
#                lon_list.append(lon)
#                lat_list.append(lat)
#                dist_list.append(spiceypy.vnorm(linedr))
#
#                if num == int(len(interval) / 2):
#                    lon_mid_list.append(lon)
#                    lat_mid_list.append(lat)
#                    dist_mid_list.append(spiceypy.vnorm(linedr))
#
#    spiops.plot(lon_mid_list, [lat_mid_list],
#                xaxis_name='Longitude [deg]',
#                yaxis_name=['Latitude [deg]'],
#                title='TGO-MEX Occultation Groundtrack for MEX-TGO Distance < {}km'.format(
#                    refval),
#                plot_height=500,
#                plot_width=900,
#                format='circle',
#                background_image=True,
#                line_width=6)
#
#    spiops.plot(lon_list, [lat_list],
#                xaxis_name='Longitude [deg]',
#                yaxis_name=['Latitude [deg]'],
#                title='TGO-MEX Occultation Groundtrack for MEX-TGO Distance < {}km'.format(
#                    refval),
#                plot_height=500,
#                plot_width=900,
#                format='circle',
#                background_image=True,
#                line_width=1)
#
#    return


def sensor_with_sectors(sensor, mk, fk=''):

    #
    # Load ROS FK and RPC IK files.
    #
    spiceypy.furnsh(mk)
    if fk:
        spiceypy.furnsh(fk)

    #
    # Get ELS IK data.
    #
    sensnm, sensid, secnum, secsiz, secsis, secfrm, secdir = read_ik_with_sectors(sensor)

    #
    # Report ELS IK data.
    #
    print('SENSOR NAIF NAME:  {}'.format(sensnm))
    print('SENSOR NAIF ID:    {}'.format(sensid))
    print('NUMBER OF SECTORS: {}'.format(secnum))

    #if secsiz != 0:
    #    print('SECTOR SIZE:       {}'.format(secsiz))
    #else:
    #    print('SECTOR SIZES:      {}'.format(secsis))


    print('REFERENCE FRAME:   {}'.format(secfrm))
    print('SECTOR DIRECTIONS: {}'.format(secdir))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    for element in secdir:

        x = element[0]
        y = element[1]
        z = element[2]

        ax.scatter(x, y, z, c='r', marker='o')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

    ax.autoscale(tight=True)

    plt.show()

    return


def get_angle(frame1, frame2, et):

    angle = 0
    angle_bool = False

    try:
        # Get the rotation matrix between two frames
        cmat = spiceypy.pxform(frame1, frame2, et)

        (angle3, angle2, angle1) = spiceypy.m2eul(cmat, 3, 2, 1)
        for tmp_angle in [angle3, angle2, angle1]:
            if np.around(tmp_angle, 2) != 0:
                angle = np.rad2deg(tmp_angle)
                angle_bool = True

    except ValueError as e:
        print(e)

    return angle, angle_bool


def get_earth_angle(frame, et, obs):
    try:
        (earth_vec, lt) = spiceypy.spkezr('EARTH', et, frame, 'LT+S', obs)
        return np.rad2deg(spiceypy.vsep([0, 0, 1], earth_vec[:3]))

    except:
        return 0


def hga_angles(sc, et):

    hga_el_frame = sc + '_HGA_EL'
    hga_az_frame = sc + '_HGA_AZ'
    hga_frame = sc + '_HGA'

    if sc == 'MPO':

        # First azimuth and then the elevation
        hga_az, hga_az_bool = get_angle('MPO_HGA_APM', hga_az_frame, et)
        if hga_az_bool:
            hga_az = -hga_az + 180  # Invert azimuth and add half revolution
        hga_el, hga_el_bool = get_angle(hga_az_frame, hga_el_frame, et)

    elif sc == 'MTM':
        return []

    else:
        hga_zero_frame = sc + '_SPACECRAFT'

        # First elevation and then the azimuth
        hga_el, hga_el_bool = get_angle(hga_zero_frame, hga_el_frame, et)
        hga_az, hga_az_bool = get_angle(hga_el_frame, hga_az_frame, et)

    hga_earth = get_earth_angle(hga_frame, et, sc)

    return [hga_az, hga_el], hga_earth


def mga_angles(sc, et):

    if sc == 'MPO':

        # First azimuth and then the elevation
        mga_az, mga_az_bool = get_angle('MPO_MGA_BOOM-H', 'MPO_MGA_BOOM', et)
        mga_el, mga_el_bool = get_angle('MPO_MGA_ZERO', 'MPO_MGA', et)
        mga_earth = get_earth_angle('MPO_MGA', et, 'MPO')

        return [mga_az, mga_el], mga_earth

    return [0, 0], 0


def solar_aspect_angles(sc, time):

    sa_frame = ''

    if sc == 'TGO':

        sa_p_frame = sc+'_SA+Z'
        sa_n_frame = sc+'_SA-Z'

    elif sc == 'MPO':

        sa_frame = sc+'_SA'

    elif sc == 'MTM':
        sa_p_frame = sc + '_SA+X'
        sa_n_frame = sc + '_SA-X'

    else:
        sa_p_frame = sc+'_SA+Y'
        sa_n_frame = sc+'_SA-Y'


    sc_id = spiceypy.bodn2c(sc)

    try:

        # If there is only one Solar Array e.g.: BEPICOLOMBO MPO
        if sa_frame:

            (sun_vec, lt) = spiceypy.spkezp(10, time, sa_frame, 'NONE', sc_id)
            saa_sa = np.rad2deg(spiceypy.vsep([1, 0, 0], sun_vec))

        else:

            (sun_vec, lt) = spiceypy.spkezp(10, time, sa_p_frame, 'NONE', sc_id)
            saa_sa_p = np.rad2deg(spiceypy.vsep([1, 0, 0], sun_vec))

            (sun_vec, lt) = spiceypy.spkezp(10, time, sa_n_frame, 'NONE', sc_id)
            saa_sa_n = np.rad2deg(spiceypy.vsep([1, 0, 0], sun_vec))

        (sun_vec, lt) = spiceypy.spkezp(10, time, sc+'_SPACECRAFT', 'NONE', sc_id)
        saa_sc_x = np.rad2deg(spiceypy.vsep([1, 0, 0], sun_vec))
        saa_sc_y = np.rad2deg(spiceypy.vsep([0, 1, 0], sun_vec))
        saa_sc_z = np.rad2deg(spiceypy.vsep([0, 0, 1], sun_vec))

    except:

        #print('No CK information for {}'.format(time))
        saa_sa, saa_sa_p, saa_sa_n, saa_sc_x, saa_sc_y, saa_sc_z = 0,0,0,0,0,0

    if sa_frame:

        return([saa_sa], [saa_sc_x, saa_sc_y, saa_sc_z])

    else:

        return ([saa_sa_p, saa_sa_n], [saa_sc_x, saa_sc_y, saa_sc_z])


def solar_array_angles(sa_frame, time):

    # Rotation axis must be angle 3 to have a range of [-pi, pi], the
    # rotation axis is derived from the FK.

    if 'MPO' in sa_frame:
        sa_zero_frame = 'MPO_SA_SADM'
    elif 'MEX' in sa_frame:
        sa_zero_frame = sa_frame + '_GIMBAL'
    else:
        sa_zero_frame = sa_frame + '_ZERO'
    try:

        #TODO This works for  JUICE only in principle.

        # Get the rotation matrix between two frames
        cmat = spiceypy.pxform(sa_frame, sa_zero_frame, time)

        (angle3, angle2, angle1) = spiceypy.m2eul(cmat, 2, 3, 1)

    except:

        # print('No CK information for {}'.format(time))
        angle3 = 0
        angle2 = 0
        angle1 = 0

    return(np.round(angle3*spiceypy.dpr(),3),
           np.round(angle2*spiceypy.dpr(),3),
           np.round(angle1*spiceypy.dpr(),3))


def structures_position(sc_frame, kernel, time):

    return


def body_distance_to_plane(body_distance, body_plane, time):

    body_1 = body_plane
    body_2 = body_distance

    if isinstance(time, str):
        time = spiceypy.utc2et(time)

    id_1 = spiceypy.bodn2c(body_1)
    id_2 = spiceypy.bodn2c(body_2)

    mat = spiceypy.pxform('MEX_SIDING_SPRING_PLANE','IAU_MARS', time)
    vec1_1 = spiceypy.mxv(mat, [1,0,0])
    vec2_1 = spiceypy.mxv(mat, [0,1,0])

    state_1 = spiceypy.spkgeo(id_2, time, 'IAU_MARS', id_1)[0]

    pos_1 = state_1[0:3]
    vel_1 = state_1[2:5]

    pos_2 = [0,0,0]

    norm_1 = np.cross(vec1_1,vec2_1)
    norm_1 = norm_1/np.linalg.norm(norm_1)

    # https://mathinsight.org/distance_point_plane

    a1, b1, c1 = norm_1[0], norm_1[1], norm_1[2]
    d1 = -1*norm_1[0]*pos_1[0] - norm_1[1]*pos_1[1] - norm_1[2]*pos_1[2]

    dist_1 = abs(a1 * pos_2[0] + b1 * pos_2[1] + c1 * pos_2[2] + d1) / np.sqrt(
        np.square(a1) + np.square(b1) + np.square(c1))

    dist_real = np.linalg.norm(pos_1)

    return dist_1, dist_real


def angle_between_planes(body_1, body_2, time):

    if isinstance(time, str):
        time = spiceypy.utc2et(time)

    mat = spiceypy.pxform('MEX_SIDING_SPRING_PLANE', 'HEE', time)
    norm_1 = spiceypy.mxv(mat, [0,0,1])

    norm_1 = norm_1 / np.linalg.norm(norm_1)

    angle = 180 - spiceypy.dpr()*spiceypy.vsep(norm_1,[0,0,1])

    return angle


def plane_ellipsoid(body_1, body_2, time):

    id_1 = spiceypy.bodn2c(body_1)
    id_2 = spiceypy.bodn2c(body_2)

    mat = spiceypy.pxform('MEX_SIDING_SPRING_PLANE','IAU_MARS', time)
    vec1 = spiceypy.mxv(mat, [1,0,0])
    vec2 = spiceypy.mxv(mat, [0,1,0])

    state1 = spiceypy.spkgeo(id_2, time, 'IAU_'+body_2, id_1)[0]
    pos1 = state1[0:3]

    plane = spiceypy.psv2pl(pos1, vec1, vec2)

    # Get the body semi-axis lenght
    (num, semi_axis) = spiceypy.bodvcd(id_2, "RADII", 3)

    a = semi_axis[0]
    b = semi_axis[1]
    c = semi_axis[2]

    try:
        ellipse = spiceypy.inedpl(a, b, c, plane)
    except:
        ellipse = 0

    return ellipse


def beta_angle(observer, target, time):

    # Provided by Bernhard Geiger

    if not isinstance(time, float):
        et = spiceypy.utc2et(time)
    else:
        et = time

    #
    # compute the Sun position relative to Mars; vector from Mars to Sun
    #
    vec_tar_sun, lt = spiceypy.spkpos( 'SUN', et, 'J2000', 'None', target)

    #
    # beta angle
    #
    sta_tar_obs, lt = spiceypy.spkezr(observer, et, 'J2000', 'None', target)

    #
    # orbital plane is defined by the cross-product of position and velocity
    # vector
    #
    vec_orbit = spiceypy.vcrss(sta_tar_obs[:3], sta_tar_obs[3:])

    #
    # the beta angle can be computed from the orbital plane and Sun vectors
    #
    beta = abs(90.-spiceypy.vsep(vec_orbit, vec_tar_sun)*spiceypy.dpr())

    return beta


def ck_coverage_timeline(metakernel, frame_list, notebook=True, html_file_name='test',
                         plot_width=975, plot_height=700):

    if notebook:
        output_notebook()
    else:
        output_file(html_file_name + '.html')

    cov_start = []
    cov_finsh = []
    kernels = []

    with open(metakernel, 'r') as f:
        for line in f:
            if '/ck/' in line and 'prelaunch' not in line:
                kernels.append(line.split('/ck/')[-1].strip().split("'")[0])
            if 'PATH_VALUES' in line and '=' in line:
                path = line.split("'")[1] + '/ck/'

    kernels = list(reversed(kernels))
    ck_kernels = []
    colors = []

    for kernel in kernels:
        for frame in frame_list:
            cov = cov_ck_ker(path + kernel, frame, support_ker=metakernel, time_format='TDB')

            if cov:
                color = "lawngreen"
                if 'MPO' in frame or 'MMO' in frame or 'MTM' in frame or 'TGO' in frame:
                    type = kernel.split('_')[3]
                    if type[2] == 'p': color = 'orange'
                    elif type[2] == 'r': color = 'green'
                    elif type[2] == 't': color = 'red'
                    elif type[2] == 'c': color = 'purple'
                    elif type[2] == 'm': color = 'blue'
                cov_start.append(cov[0])
                cov_finsh.append(cov[-1])
                ck_kernels.append(kernel)
                colors.append(color)

    spiceypy.furnsh(metakernel)
    date_format = 'UTC'
    start_dt =[]
    finsh_dt =[]
    for element in cov_start:
        start_dt.append(et_to_datetime(element, date_format))
    for element in cov_finsh:
        finsh_dt.append(et_to_datetime(element, date_format))



    source = ColumnDataSource(data=dict(start_dt=start_dt,
                                        finsh_dt=finsh_dt,
                                        ck_kernels=ck_kernels))

    title = "CK Kernels Coverage"
    if 'ops' in metakernel.lower():
        title += ' - OPS Metakernel'
    elif 'plan' in metakernel.lower():
        title += ' - PLAN Metakernel'
    p = figure(y_range=ck_kernels, plot_height=plot_height, plot_width=plot_width, title=title, )
    p.hbar(y=ck_kernels, height=0.2, left=start_dt, right=finsh_dt, color=colors)

    labels = LabelSet(x='start_dt', y='ck_kernels', text='ck_kernels', level='glyph',
                  x_offset=-2, y_offset=5, source=source, render_mode='canvas')

    p.xaxis.formatter = DatetimeTickFormatter(seconds=["%Y-%m-%d %H:%M:%S"],
                                              minsec=["%Y-%m-%d %H:%M:%S"],
                                              minutes=["%Y-%m-%d %H:%M:%S"],
                                              hourmin=["%Y-%m-%d %H:%M:%S"],
                                              hours=["%Y-%m-%d %H:%M:%S"],
                                              days=["%Y-%m-%d %H:%M:%S"],
                                              months=["%Y-%m-%d %H:%M:%S"],
                                              years=["%Y-%m-%d %H:%M:%S"])

    p.xaxis.major_label_orientation = 0#pi/4
    p.yaxis.visible = False
    p.xaxis.axis_label_text_font_size = "5pt"


    p.add_layout(labels)

    show(p)


def spk_coverage_timeline(metakernel, sc_list, notebook=True, html_file_name='test',
                         plot_width=975, plot_height=500):

    if notebook:
        output_notebook()
    else:
        output_file(html_file_name + '.html')

    cov_start = []
    cov_finsh = []
    kernels = []

    with open(metakernel, 'r') as f:
        for line in f:
            if '/spk/' in line and 'prelaunch' not in line:
                kernels.append(line.split('/spk/')[-1].strip().split("'")[0])
            if 'PATH_VALUES' in line and '=' in line:
                path = line.split("'")[1] + '/spk/'

    kernels = list(reversed(kernels))
    spk_kernels = []
    colors = []

    for kernel in kernels:
        for sc in sc_list:
            cov = cov_spk_ker(path+kernel, sc.upper(), support_ker=metakernel,
                                time_format='TDB')
            if cov:
                color = "lawngreen"
                if 'MPO' in sc or 'MMO' in sc or 'MTM' in sc:
                    type = kernel.split('_')[2]
                    if type[2] == 'p':
                        color = 'orange'
                    elif type[2] == 'r':
                        color = 'green'
                    elif type[2] == 't':
                        color = 'red'
                    elif type[2] == 'c':
                        color = 'purple'
                    elif type[2] == 'm':
                        color = 'blue'
                cov_start.append(cov[0][0])
                cov_finsh.append(cov[0][-1])
                spk_kernels.append(kernel)
                colors.append(color)

    spiceypy.furnsh(metakernel)
    date_format = 'UTC'
    start_dt =[]
    finsh_dt =[]
    for element in cov_start:
        start_dt.append(et_to_datetime(element, date_format))
    for element in cov_finsh:
        finsh_dt.append(et_to_datetime(element, date_format))



    source = ColumnDataSource(data=dict(start_dt=start_dt,
                                        finsh_dt=finsh_dt,
                                        spk_kernels=spk_kernels))

    title = "SPK Kernels Coverage"
    if 'ops' in metakernel.lower():
        title += ' - OPS Metakernel'
    elif 'plan' in metakernel.lower():
        title += ' - PLAN Metakernel'
    p = figure(y_range=spk_kernels,  plot_height=plot_height, plot_width=plot_width, title=title, )
    p.hbar(y=spk_kernels, height=0.2, left=start_dt, right=finsh_dt, color=colors)

    labels = LabelSet(x='start_dt', y='spk_kernels', text='spk_kernels', level='glyph',
                  x_offset=-2, y_offset=5, source=source, render_mode='canvas')

    p.xaxis.formatter = DatetimeTickFormatter(seconds=["%Y-%m-%d %H:%M:%S"],
                                              minsec=["%Y-%m-%d %H:%M:%S"],
                                              minutes=["%Y-%m-%d %H:%M:%S"],
                                              hourmin=["%Y-%m-%d %H:%M:%S"],
                                              hours=["%Y-%m-%d %H:%M:%S"],
                                              days=["%Y-%m-%d %H:%M:%S"],
                                              months=["%Y-%m-%d %H:%M:%S"],
                                              years=["%Y-%m-%d %H:%M:%S"])

    p.xaxis.major_label_orientation = 0#pi/4
    p.yaxis.visible = False
    p.xaxis.axis_label_text_font_size = "5pt"

    p.add_layout(labels)

    show(p)

    spiceypy.unload(metakernel)


#
#   The camera has no distortion;  the image of a point
#   is determined by the intersection of the focal plane
#   and the line determined by the point and the camera's
#   focal point.
#
#   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/inrypl_c.html
#
def target_center_pixel(et, camera, target, target_frame):

    camera_id = spiceypy.bodn2c(camera)

    focal_lenght = spiceypy.gdpool('INS{}_FOCAL_LENGTH'.format(camera_id), 0, 80)[0]
    focal_lenght /= 1000000 # in milimiters (original routine for OSIRIS was 717.322)


    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id, 100)
    visible = spiceypy.fovtrg(camera, target, 'POINT', target_frame, 'LT+S', camera, et)


    if visible:

        (ptarg, lt) = spiceypy.spkpos(camera, et, frame, 'LT+S', target)
        (ptarg, norm) = spiceypy.unorm(ptarg)

        focus = spiceypy.vscl(focal_lenght, [0,0,1])


        #
        #   The camera's focal plane contains the origin in
        #   camera coordinates, and the z-vector is orthogonal
        #   to the plane.  Make a CSPICE plane representing
        #   the focal plane.
        #
        focal_plane = spiceypy.nvc2pl([0,0,1], 0.)

        #
        #    The image of the target body's center in the focal
        #    plane is defined by the intersection with the focal
        #    plane of the ray whose vertex is the focal point and
        #    whose direction is dir.
        #
        (nxpts, image_focal) = spiceypy.inrypl(focus, ptarg, focal_plane)

        if nxpts != 1:
            print('Something went wrong')
            return 0

    else:
        print('{}: {} Center is not in the image'.format(spiceypy.et2utc(et, 'ISOC', 3, 20),target))
        return

    pixel_size = spiceypy.gdpool('INS{}_PIXEL_SIZE'.format(camera_id), 0, 80)[0]
    ccd_center = spiceypy.gdpool('INS{}_CCD_CENTER'.format(camera_id), 0, 80)

    image  = image_focal * (1000000000/pixel_size)  # Pixel Size in Microns

    return (image_focal, (ccd_center[1]+image[1], ccd_center[0]+image[0]))


def pixel_center_distance(et, camera, pixel_x, pixel_y):


    pix_x = 1024 - pixel_y
    pix_y = 1024 - pixel_x

    camera_name = camera
    camera_id = spiceypy.bodn2c(camera_name)


    focal_lenght = 717.322/1000000

    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id, 100)


    (image_focal, target_pixel) = target_center_pixel(et, camera, report=False)

    center_x = image_focal[0]
    center_y = image_focal[1]

    tar_cent_vec = [center_x, center_y, bsight[2]*focal_lenght]
    pix_vec = [pix_x*13.5/1000000000, pix_y*13.5/1000000000, bsight[2]*focal_lenght]

    (ptarg, lt) = spiceypy.spkpos(camera, et, frame, 'LT+S', '67P/C-G')
    (ptarg, norm) = spiceypy.unorm(ptarg)


    pixel_vector = [13.5/1000000000, 0, bsight[2]*focal_lenght]

    pix_cent_dist = np.tan(spiceypy.vsep(pix_vec,tar_cent_vec))*norm

    pixel_size = np.tan(spiceypy.vsep([0,0,bsight[2]*focal_lenght],pixel_vector))*norm

    #print('Comet offset', spiceypy.dpr()*spiceypy.vsep(bsight*focal_lenght,tar_cent_vec))
    #print('Pixe-Comet offset', spiceypy.dpr()*spiceypy.vsep(pix_vec,tar_cent_vec))


    return pix_cent_dist, target_pixel, pixel_size


def simulate_image(utc, camera, mission_targets, camera_spk=False,
                    pixel_lines=False, pixel_samples=False, dsks=False,
                    generate_image=False, report=False, name=False,
                    illumination=True,  metakernel='', unload_kernels=True, log=False):
    '''

    :param utc: Image acquisition time in UTC format e.g.: 2016-01-01T00:00:00
    :type utc: str
    :param metakernel: SPICE Kernel Dataset     Meta-Kernel
    :type metakernel: str
    :param camera: Name of the camera to be used. Usually found in the
    instrument kernel (IK) e.g.: 'ROS_NAVCAM-A'
    :type camera: str
    :param mission_targets: Targets of the observation, e.g.:'67P/C-G'
    :type mission_targets: list
    :param pixel_lines: Number of pixel lines usually provided by the IK.
    :type pixel_lines: int
    :param pixel_samples: Number of pixel samples per line usually provided by
    the IK.
    :type pixel_samples: int
    :param dsk: Digital Shape Model to be used for the computation. Not required
    of included in the Meta-Kernel.
    :type dsk: str
    :param generate_image: Flag to determine whether if the image is saved or
    plotted.
    :type generate_image: bool
    :param plot_image: Flag to determine whether if the image is to be plotted or
    plotted.
    :type generate_image: bool
    :param report: Flag for processing report.
    :type generate_image: bool
    :param name: Name to be provided to the image
    :type generate_image: str
    :return: Name of the output image
    :rtype: str
    '''
    if metakernel:
        spiceypy.furnsh(metakernel)
    if dsks:
        for dsk in dsks:
            spiceypy.furnsh(dsk)
    et = spiceypy.utc2et(utc)

    #
    # We retrieve the camera information using GETFOV. More info available:
    #
    #   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html
    #
    camera_name = camera
    camera_id = spiceypy.bodn2c(camera_name)
    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id, 100)


    #
    # TODO: In the future all the sensors should be epehmeris objects, see
    # https://issues.cosmos.esa.int/socci/browse/SPICEMNGT-77
    #
    if not camera_spk:
        if camera.split('_')[0] == 'ROS':
            observer = 'ROSETTA'
        elif camera.split('_')[0] == 'MEX':
            observer = 'MEX'
        elif camera.split('_')[0] == 'VEX':
            observer = 'VEX'
        elif camera.split('_')[0] == 'JUICE':
            observer = 'JUICE'
        elif camera.split('_')[0] == 'HERA':
            observer = 'HERA'
        elif camera.split('_')[0] == 'MTM':
            observer = 'MTM'
    else:
        if isinstance(camera_spk, str):
                observer = camera_spk
        else:
            observer = camera

    #
    # We check if the resolution of the camera has been provided as an input
    # if not we try to obtain the resolution of the camera from the IK
    #
    if not pixel_lines or not pixel_samples:
        try:
            pixel_samples = int(spiceypy.gdpool('INS'+str(camera_id) + '_PIXEL_SAMPLES',0,1))
            pixel_lines = int(spiceypy.gdpool('INS' + str(camera_id) + '_PIXEL_LINES',0,1))
        except:
            pass
            print("PIXEL_SAMPLES and/or PIXEL_LINES not defined for "
                  "{}".format(camera))
            return

    #
    # We generate a matrix using the resolution of the framing camera as the
    # dimensions of the matrix
    #
    nx, ny = (pixel_lines, pixel_samples)
    x = np.linspace(bounds[0][0], bounds[2][0], nx)
    y = np.linspace(bounds[0][1], bounds[2][1], ny)

    #
    # We define the matrices that will be used as outputs and the
    #
    phase_matrix = np.zeros((nx, ny))
    emissn_matrix = np.zeros((nx, ny))
    solar_matrix = np.zeros((nx, ny))
    target_matrix = np.zeros((nx, ny))

    #
    # Now we look for additional targets.
    #
    targets_frames = []
    methods = []
    targets = []

    for target in mission_targets:

        try:
            target_frame = target2frame(target)

            if dsks:
                for dsk in dsks:
                    ids = spiceypy.dskobj(dsk)
                    if spiceypy.bodn2c(target) in ids:
                        method = 'DSK/UNPRIORITIZED'
                        break
                    else:
                        method = 'ELLIPSOID'
            else:
                method = 'ELLIPSOID'

            visible = spiceypy.fovtrg(camera, target, 'POINT', target_frame, 'NONE',
                                      observer, et)

            if visible:
                print('{} center is visible'.format(target))

            targets.append(target)
            targets_frames.append(target_frame)
            methods.append(method)
        except Exception as e:
            print(e)
            pass

    r = []
    for target in targets:

        r.append(np.linalg.norm(
                spiceypy.spkpos(target, et, 'J2000', 'NONE', observer)[
                    0]))
    #
    # If we have targets we order the list in target proximity order
    #
    if len(r) > 1:
        targets = np.asarray(targets)
        targets = targets[np.argsort(r)]
        targets_frames = np.asarray(targets_frames)
        targets_frames = targets_frames[np.argsort(r)]
        methods = np.asarray(methods)
        methods = methods[np.argsort(r)]

    #
    # For each pixel we compute the possible intersection with the target, if
    # the target is intersected we then compute the illumination angles. We
    # use the following SPICE APIs: SINCPT and ILLUMF
    #
    #   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sincpt_c.html
    #   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/illumf_c.html
    #
    isvisible, isiluminated = [], []
    for i in range(0, len(x), 1):
        for j in range(0, len(y), 1):



            #
            # List of pixel's boresight
            #
            ibsight = [x[i], y[j], bsight[2]]

            #
            # We do another loop in order to determine if we have other
            # 'targets' in addition to the 'main' target
            #
            for k in range(0, len(targets), 1):
                try:
                    (spoint, trgepc, srfvec ) = spiceypy.sincpt(methods[k], targets[k], et,
                                                targets_frames[k], 'NONE', observer, frame, ibsight)

                    target_matrix[i, j] = 255


                    (trgenpc, srfvec, phase, solar,
                     emissn, visiblef, iluminatedf) = spiceypy.illumf(methods[k], targets[k], 'SUN', et,
                                                      targets_frames[k], 'LT+S', observer, spoint)

                    emissn_matrix[i, j] = emissn
                    phase_matrix[i, j] = phase

                    #
                    # Add to list if the point is visible to the camera
                    #
                    if visiblef == True:
                        isvisible.append(visiblef)

                    #
                    # Add to list if the point is illuminated and seen by the camera
                    #
                    if iluminatedf == True:

                        isiluminated.append(iluminatedf)
                        solar_matrix[i, j] = solar

                    else:
                        #
                        # And we set the not illuminated pixels with np.pi/2
                        #
                        solar_matrix[i, j] = np.pi/2

                    break

                except:
                    pass

                    #
                    # If SINCPT raises an error, we set that we see nothing in
                    # the pixel.
                    #
                    emissn_matrix[i,j] = 0
                    phase_matrix[i,j] = np.pi
                    solar_matrix[i,j] = np.pi/2

    #
    # We transform the matrix from illumination angles to greyscale [0-255]
    #
    if solar_matrix.max() == solar_matrix.min():
        rescaled = solar_matrix
    else:
        rescaled = (255 / (solar_matrix.max()-solar_matrix.min()) * (solar_matrix - solar_matrix.min())).astype(np.uint8)

    rescaled = - np.flip(rescaled, 0) + 255

    #
    # We generate the plot
    #
    if generate_image:
        if not name:

            name_illum = '{}_{}.PNG'.format(camera.upper(),
                                            utc.upper())
            name_tar = '{}_{}_TAR.PNG'.format(camera.upper(),
                                            utc.upper())
        else:

            name_illum = '{}_{}_{}.PNG'.format(name.upper(),
                                         camera.upper(),
                                         utc.upper())
            name_tar = '{}_{}_{}_TAR.PNG'.format(name.upper(),
                                         camera.upper(),
                                         utc.upper())


        if np.count_nonzero(target_matrix) >= 1.0:
            if illumination:
                imageio.imwrite(name_illum, rescaled)
            else:
                imageio.imwrite(name_tar, target_matrix)

    if not generate_image:
        plt.imshow(rescaled, cmap='gray')
        plt.axis('off')
        plt.show()

    if report:
        print('{} {} {} {} {}'.format(utc, pixel_samples*pixel_lines, np.count_nonzero(target_matrix), len(isvisible), camera))

    #if report:
    #    print('Pixel report for {} w.r.t {} @ {}'.format(camera,target,utc))
    #    print('   Total number of pixels: ', pixel_samples*pixel_lines)
    #    print('   Illuminated pixels:     ', len(isvisible))
    #    print('   Hidden pixels:          ', pixel_samples*pixel_lines - len(isvisible))
    #    print('   Shadowed points:        ', pixel_samples*pixel_lines - len(isiluminated))

    if log:
        with open(log, "a") as f:
            f.write('{} {} {} {} {}\n'.format(utc, pixel_samples*pixel_lines, np.count_nonzero(target_matrix), len(isvisible), camera))

    if unload_kernels:
        spiceypy.kclear()

    return name


def sc_dsk_view(utc,mk, dsks, observer, sc_targets, sc_frames=False,
                pixels=150, name=False, generate_image=True, illumination=True,
                show3Dplot=False, unload_kernels=False):

    if not sc_frames:
        sc_frames = sc_targets

    mpl.rcParams['figure.figsize'] = (26.0, 26.0)

    spiceypy.furnsh(mk)

    for dsk in dsks:
        spiceypy.furnsh(dsk)

    utcstr = utc.replace(':','')
    et = spiceypy.utc2et(utc)

    nx, ny = (pixels, pixels)   # resolution of the image
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xv, yv = np.meshgrid(x, y)


    solar_matrix = np.zeros((nx, ny))
    target_matrix = np.zeros((nx, ny))

    isvisible, isiluminated = [], []
    r, lt = spiceypy.spkpos(observer, et, 'J2000', 'NONE', sc_targets[0])

    #
    # We define a 'Nadir frame' w.r.t. J000 to make it general regardless of
    #
    #
    zN = r
    zN = zN/np.linalg.norm(zN)
    xN = np.array([1,0,0]) - np.dot(np.dot([1,0,0], zN), zN)/np.linalg.norm(zN)**2
    xN = xN/np.linalg.norm(xN)
    yN = np.cross(zN, xN)
    yN = yN/np.linalg.norm(yN)
    RotM = np.linalg.inv(np.array([xN, yN, zN]))

    spoints = []
    flag = False
    f = 0.5 # Factor for the FOV

    for i, x in enumerate(xv):
        for j, y in enumerate(yv):
            dpxy = [x[i], y[i], -np.linalg.norm(r)*1000*f]
            ibsight = spiceypy.mxv(RotM, dpxy)

            #
            # We do another loop in order to determine if we have other
            # 'targets' in addition to the 'main' target
            #
            spoint_per_target = []
            distance_per_target = []
            for k in range(0, len(sc_targets), 1):
                try:
                    (spoint, trgepc, srfvec) = spiceypy.sincpt('DSK/UNPRIORITIZED', sc_targets[k], et, sc_frames[k], 'NONE', observer, 'J2000', ibsight)
                    spoint_per_target.append(spoint)
                    distance_per_target.append(spiceypy.vnorm(srfvec))

                except:
                    spoint_per_target.append([])
                    distance_per_target.append([])
                    pass

            for spoint in spoint_per_target:
                if not spoint.__class__ == list:
                    flag = True

            if flag:
                #
                # we get the minimum distance and the range.
                #
                distance_per_target_floats = []
                for element in distance_per_target:
                    if not isinstance(element, list):
                        distance_per_target_floats.append(element)
                min_dist = min(distance_per_target_floats)

                for k in range(0, len(sc_targets), 1):
                    if min_dist == distance_per_target[k]:
                        target = sc_targets[k]
                        tar_frame = sc_frames[k]
                        spoints.append(spoint_per_target[k])
                        spoint = spoint_per_target[k]

                        try:
                            (trgenpc, srfvec, phase, solar,
                            emissn, visiblef, iluminatedf) = spiceypy.illumf('DSK/UNPRIORITIZED', target, 'SUN', et, tar_frame, 'NONE', observer, spoint)

                            if visiblef == True:
                                isvisible.append(visiblef)

                                target_matrix[i, j] = 255
                            if iluminatedf == True:
                                isiluminated.append(iluminatedf)

                                if solar > np.pi / 2:
                                    solar_matrix[i, j] = np.pi - solar
                                else:
                                    solar_matrix[i, j] = solar
                            else:
                                solar_matrix[i, j] = np.pi / 2  # not illuminated
                        except:
                            solar_matrix[i, j] = np.pi / 2
            else:
                solar_matrix[i, j] = np.pi / 2


            flag = False

    if show3Dplot:
        spoints = np.asarray(spoints)
        fig = plt.figure(figsize=(9, 9))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(spoints[:, 0], spoints[:, 1], spoints[:, 2], marker='.')
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title('')
        plt.axis('equal')
        plt.show()

    print('total number of points: ', pixels*pixels)
    print('occulted points: ', pixels*pixels - len(isvisible))
    print('not iluminated points: ', pixels*pixels - len(isiluminated))

    #
    # We transform the matrix from illumination angles to greyscale [0-255]
    #
    if illumination:
        rescaled = (255 / (solar_matrix.max()-solar_matrix.min()) * (solar_matrix - solar_matrix.min())).astype(np.uint8)
        rescaled = - np.flip(rescaled, 0) + 255
    else:
        rescaled = target_matrix


    #
    # We generate the plot
    #
    if generate_image:
        if not name:

            name = '{}_{}.PNG'.format(sc_targets[0].upper(),
                                         utcstr.upper())
        else:

            name = '{}_{}_{}.PNG'.format(name.upper(),
                                         sc_targets[0].upper(),
                                         utcstr.upper())


        imageio.imwrite(name, rescaled)

    else:
        plt.imshow(rescaled, cmap='gray')
        plt.axis('off')
        plt.show()


    if unload_kernels:
        spiceypy.kclear()

    return


def getXYforPlanet(time_et, planet, camera, observer=''):
    """
    compute for all time instances in this class the xy position of a planet
    within a camera_name field-of-view. If not visible, return (-1,-1).
    If planet is visible, also return the size of the planet in pixels.
    Routine is tested for planet_name=Earth
    """
    #
    camera_id = spiceypy.bodn2c(camera)

    #
    # In case we do not have an SPK for the camera
    #
    if not observer:
        observer = camera

    r_planet = (spiceypy.bodvrd(planet, 'RADII', 3))[1][0]
    #
    # get instrument related info
    #

    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id,
                                                               100)
    mat = spiceypy.pxform(frame, 'J2000', time_et)
    for bs in range(0, 4):
        bounds[bs, :] = spiceypy.mxv(mat, bounds[bs, :])

    [pos, ltime] = spiceypy.spkpos(planet, time_et, 'J2000',
                                    'LT+S', observer)
    visible = spiceypy.fovray(camera, pos, 'J2000', 'S', 'MPO',
                               time_et)
    #
    # only get detector position, if target is visible
    #
    x = 0.0
    y = 0.0
    s = 0.0
    if visible:
        hit = []
        for p in range(0, 4):
            #
            # compute the plane that is build up by the coordinate origin and two FOV corner vectors
            #
            plane = spiceypy.psv2pl([0, 0, 0], bounds[p, :],
                                 bounds[(p + 1) % 4, :])
            #
            # compute the projection of the target vector onto that plane
            #
            vout = (spiceypy.unorm(spiceypy.vprjp(pos, plane)))[0]
            #
            # calculate the angle between this vector and the original corner vectors
            #
            alpha = spiceypy.vsep(bounds[p, :], vout)
            beta = spiceypy.vsep(bounds[(p + 1) % 4, :], vout)
            #
            # the ratio of these angles also give the ratio of the detector on the edge
            # of the field of view, in a first approximation, average of the two opposite
            # FOV corner values: these are the x,y coordinates on the detector
            hit.append(1024 * alpha / (alpha + beta))

        # get intersection of the points
        (x, y) = findIntersection(hit[0], 0, hit[1], 1023, 0,
                                       hit[1], 1023, hit[2])
        size = 2 * r_planet * 500 / (
                        np.tan(35. * 2 * np.pi / 360.) * spiceypy.vnorm(pos))

    else:
        print('Planet {} not visible by {} at {}'.format(planet, camera, spiceypy.et2utc(time_et,'ISOC',1,25)))
        return (False, False, False, False)

    return (time, x, y, size)


def pixel_radec(time, camera, pixel, units='radians'):
    """
    Obtain the Right Ascension and Declination in J2000 of a pixel of a given
    sensor at a given time.

    @param time: Input time in UTC
    @type time:  str
    @param camera: SPICE name for the camera sensor (requires IK kernel)
    @type camera: str
    @param pixel: Pixel location in the camera sensor CCD. Provided as two
    values [x y]
    @type pixel: list
    @param units: Angular units for Right Ascension and Declination: radians or
    degrees
    @type units: str
    @return: Right Ascension and Declination in the indicated units
    @rtype: tuple
    """

    et = spiceypy.utc2et(time)

    #
    # We retrieve the camera information using GETFOV. More info available:
    #
    #   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html
    #
    camera_name = camera
    camera_id = spiceypy.bodn2c(camera_name)
    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id, 100)

    pixel_samples = \
        int(spiceypy.gdpool(''.join(('INS', str(camera_id), '_PIXEL_SAMPLES')),
                            0, 1))
    pixel_lines = \
        int(spiceypy.gdpool(''.join(('INS', str(camera_id), '_PIXEL_LINES')), 0,
                            1))

    #
    # We generate a matrix using the resolution of the framing camera as the
    # dimensions of the matrix
    #
    nx, ny = (pixel_samples, pixel_lines)
    x = np.linspace(bounds[0][0], bounds[2][0], nx)
    y = np.linspace(bounds[0][1], bounds[2][1], ny)

    (i, j) = pixel

    #
    # List of pixel's boresight
    #
    ibsight = [x[i], y[j], bsight[2]]

    mat = spiceypy.pxform(frame, 'J2000', et)
    ibsight = spiceypy.mxv(mat, ibsight)
    (r, ra, dec) = spiceypy.recrad(ibsight)

    if units == 'degrees':
        ra = np.degrees(ra)
        dec = np.degrees(dec)

    return ra, dec


def radec_in_fov(time, ra, dec, camera, observer=False, units='degrees'):
    """
    Determine whether if a given Right Ascension and Declination coordinate in
    J2000 (ultimately a given Star) is present at a given sensor Field-of-View
    at a given time.

    @param time: Input time in UTC
    @type time:  str
    @param ra: Right Ascension in the indicated units (w.r.t J2000)
    @type ra: float
    @param dec: Declination in the indicated units (w.r.t J2000)
    @type dec: float
    @param camera: SPICE name for the camera sensor (requires IK kernel)
    @type camera: str
    @param observer: SPICE name for the camera sensor position
    @type observer: str
    @param units: Angular units for Right Ascension and Declination: radians or
    degrees
    @type units: str
    @return: True/False if the RA, DEC is in the Field-of-View
    @rtype: bool
    """

    #
    # If an observer is not provided then it is assumed that the observer is the
    # camera itself. Please note that this then requires the strctures SPK to
    # be present in the meta-kernel (or the loaded kernels)
    #
    if not observer:
        observer = camera

    if units == 'degrees':
        ra = np.radians(ra)
        dec = np.radians(dec)

    et = spiceypy.utc2et(time)

    #
    # Create a unit direction vector pointing from the given S/C
    # to the specified star. For details on corrections such
    # as parallax, please see the example in GFRFOV.
    #
    raydir = spiceypy.radrec(1.0, ra, dec)

    #
    # Is the star in the field-of-view of the given sensor?
    #
    visible = spiceypy.fovray(camera, raydir, 'J2000', 'S', observer, et)

    if visible:
        return True
    else:
        return False


def gf_radec_in_fov(start_time, finish_time, ra, dec, camera, step=60,
                    units='degrees', observer=False):
    """
    This functions provides the time windows for which a given Right Ascension
    and Declination in J2000 (ultimately a given Star) is present in the given
    camera FOV for a given start and finish UTC times.

    @param start_time: Search time window start in UTC
    @type start_time: str
    @param finish_time: Search time window finish in UTC
    @type finish_time: str
    @param ra: Right Ascension in the indicated units (w.r.t J2000)
    @type ra: float
    @param dec: Declination in the indicated units (w.r.t J2000)
    @type dec: float
    @param camera: SPICE name for the camera sensor (requires IK kernel)
    @type camera: str
    @param step: Step with which the search will be performed in seconds
    @type step: float
    @param units: Angular units for Right Ascension and Declination: radians or
    degrees
    @type units: str
    @param observer: SPICE name for the camera sensor position
    @type observer: str
    @return: List of Time Windows
    @rtype: list
    """

    #
    # If an observer is not provided then it is assumed that the observer is the
    # camera itself. Please note that this then requires the strctures SPK to
    # be present in the meta-kernel (or the loaded kernels)
    #
    if not observer:
        observer = camera

    if units == 'degrees':
        ra = np.radians(ra)
        dec = np.radians(dec)

    et_start = spiceypy.utc2et(start_time)
    et_finish = spiceypy.utc2et(finish_time)

    #
    # Create a unit direction vector pointing from the given S/C
    # to the specified star. For details on corrections such
    # as parallax, please see the example in GFRFOV.
    #
    raydir = spiceypy.radrec(1.0, ra, dec)

    MAXIVL = 10000
    MAXWIN = 2 * MAXIVL
    TDBFMT = 'YYYY MON DD HR:MN:SC.### (TDB) ::TDB'

    # Initialize the "confinement" window with the interval
    # over which we'll conduct the search.
    cnfine = stypes.SPICEDOUBLE_CELL(2)
    spiceypy.wninsd(et_start, et_finish, cnfine)

    #
    # In the call below, the maximum number of window
    # intervals gfposc can store internally is set to MAXIVL.
    # We set the cell size to MAXWIN to achieve this.
    #
    reswin = stypes.SPICEDOUBLE_CELL(MAXWIN)

    #
    # Now search for the time period, within our confinement
    # window, during which the RA, DEC ray is in the camera FOV.
    #
    # VARIABLE  I/O  DESCRIPTION
    # --------  ---  --------------------------------------------------
    # INST       I   Name of the instrument.
    # RAYDIR     I   Ray's direction vector.
    # RFRAME     I   Reference frame of ray's direction vector.
    # ABCORR     I   Aberration correction flag.
    # OBSRVR     I   Name of the observing body.
    # STEP       I   Step size in seconds for finding FOV events.
    # CNFINE     I   SPICE window to which the search is restricted.
    # RESULT     O   SPICE window containing results.
    spiceypy.gfrfov(camera, raydir, 'J2000', 'S', camera, step, cnfine, reswin)

    #
    # The function wncard returns the number of intervals
    # in a SPICE window.
    #
    winsiz = spiceypy.wncard(reswin)

    if winsiz == 0:
        print('No events were found.')

    else:

        #
        # Display the event time periods.
        #
        print('Time Windows for RA={0:f}, DEC={1:f} [DEG] in '
              '{2:s} FOV:'.format(np.degrees(ra), np.degrees(dec), camera))

        #
        # Store the values in a list
        #
        intervals = []

        for i in range(winsiz):
            #
            # Fetch the start and stop times of the ith interval from the search
            # result window reswin.
            #
            [intbeg, intend] = spiceypy.wnfetd(reswin, i)
            intervals.append([intbeg, intend])

            #
            # Convert the start and finish times to a TDB calendar string.
            #
            print(spiceypy.timout(intbeg, TDBFMT), ',',
                  spiceypy.timout(intend, TDBFMT))

    return intervals


def radec2pixel(time, ra, dec, camera, observer=False, units='degrees'):
    """
    This function determines the pixel location for a given camera of a
    Right Ascension and Declination coordinate in J2000 (ultimately a star
    position) at a given time.

    @param time: Input time in UTC
    @type time:  str
    @param ra: Right Ascension in the indicated units (w.r.t J2000)
    @type ra: float
    @param dec: Declination in the indicated units (w.r.t J2000)
    @type dec: float
    @param camera: SPICE name for the camera sensor (requires IK kernel)
    @type camera: str
    @param observer: SPICE name for the camera sensor position
    @type observer: str
    @param units: Angular units for Right Ascension and Declination: radians or
    degrees
    @type units: str
    @return: Pixel location of a given Right Ascension and Declination
    (if present in the FOV)
    @rtype: tuple
    """

    #
    # We first check whether if the RA, DEC is in the FOV
    #
    if radec_in_fov(time, ra, dec, camera, observer=observer, units='degrees'):

        (ra_matrix, dec_matrix) = camera_radec(time, camera, units=units)

        #
        # Now we look for the pixel location of the RA, DEC
        #
        (ra_idx, ra_value) = findNearest(ra_matrix, ra)
        (dec_idx, dec_value) = findNearest(dec_matrix, dec)

        #
        # If the indexes are not the same, indicate and provide both
        #
        if ra_idx != dec_idx:
            print('RA and DEC indexes are not the same: {}, {}'.format(ra_idx,
                                                                       dec_idx))

        return ra_idx

    else:
        return 'RA, DEC are not in the FOV'


def camera_radec(time, camera, units='radians', plot=False):
    """
    This function provides two mesh grids with Right Ascensions and Declinations
    for a given camera at a given time.

    @param time: Input time in UTC
    @type time:  str
    @param camera: SPICE name for the camera sensor (requires IK kernel)
    @type camera: str
    @param units: Angular units for Right Ascension and Declination: radians or
    degrees
    @type units: str
    @param plot: Indicate whether if the mesh grids should be ploted or not
    @type plot: bool
    @return:
    @rtype:
    """

    et = spiceypy.utc2et(time)

    #
    # We retrieve the camera information using GETFOV. More info available:
    #
    #   https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html
    #
    camera_name = camera
    camera_id = spiceypy.bodn2c(camera_name)
    (shape, frame, bsight, vectors, bounds) = spiceypy.getfov(camera_id, 100)

    pixel_samples = \
        int(spiceypy.gdpool(''.join(('INS', str(camera_id), '_PIXEL_SAMPLES')),
                            0, 1))
    pixel_lines = \
        int(spiceypy.gdpool(''.join(('INS', str(camera_id), '_PIXEL_LINES')), 0,
                            1))

    #
    # We generate a matrix using the resolution of the framing camera as the
    # dimensions of the matrix
    #
    nx, ny = (pixel_samples, pixel_lines)
    x = np.linspace(bounds[0][0], bounds[2][0], nx)
    y = np.linspace(bounds[0][1], bounds[2][1], ny)
    xv, yv = np.meshgrid(x, y)

    #
    # We define the matrices that will be used as outputs and the
    #
    ra_matrix = np.zeros((nx, ny))
    dec_matrix = np.zeros((nx, ny))

    #
    # For each pixel we compute the RA/DEC in J2000.
    #
    for i, x in enumerate(xv):
        for j, y in enumerate(yv):
            #
            # List of pixel's boresight
            #
            ibsight = [x[i], y[j], bsight[2]]

            mat = spiceypy.pxform(frame, 'J2000', et)
            ibsight = spiceypy.mxv(mat, ibsight)
            (r, ra, dec) = spiceypy.recrad(ibsight)

            if units == 'degrees':
                ra_matrix[i, j] = np.degrees(ra)
                dec_matrix[i, j] = np.degrees(dec)

    if plot:
        plt.imshow(ra_matrix)
        plt.axis('off')
        plt.show()
        plt.imshow(dec_matrix)
        plt.axis('off')
        plt.show()

    return ra_matrix, dec_matrix


def groundtrack_velocity(time, observer, target, target_frame):

    radii = spiceypy.bodvrd(target, 'RADII', 3)
    re = radii[0]  # equatorial
    rp = radii[2]  # polar
    f = (re-rp)/re # target flattening factor

    #
    # compute the state vector (position(1: 3), speed(4: 6)) for each second
    #
    (state, lt) = spiceypy.spkezr(observer, time, target_frame, 'NONE', target)
    (lon, lat, alt) = spiceypy.recgeo(state[:3], re, f)

    #
    #  for each second compute the Jacobian Matrix to convert the speed to
    #  a body-fixed reference frame (in this case in Geodetic coordinates)
    #
    jacobi = spiceypy.dgeodr(state[0], state[1], state[2],re, f)
    geodetic_speed = spiceypy.mxv(jacobi, state[4:])

    #
    #  from the geodetic speed extract the radial component
    #  radial_speed[i] = geodetic_speed[i,2]
    #
    radial_component = geodetic_speed[-1]

    #
    #  the tangential speed is obtained as product of the local radius of the
    #  observed body with the tangential angular speed:
    #
    #  latitudinal component
    #  ^  x
    #  | / tangential component
    #  |/
    #  o---> longitudinal component (the cos is to compensate the "shrinking"
    #        of longitude incerasing the latitude)
    local_radius = (re * rp) / (
        np.power(re ^ 2 * math.sin(lat) ^ 2) + (rp ^ 2 * math.cos(lat) ^ 2),0.5)

    #tangential_speed = local_radius * np.power(
    #    (REFORM(geodetic_speed[i, 0]) * cos(lat[i])) ^ 2 + (
    #        REFORM(geodetic_speed[i, 1])) ^ 2)


    return

def roll(time):

    # Rotation axis must be angle 3 to have a range of [-pi, pi], the
    # rotation axis is derived from the FK.
    try:
    # Get the rotation matrix between two frames
        cmat = spiceypy.pxform('SOLO_SRF', 'SOLO_ORBIT_NORM', time)

        (angle3, angle2, angle1) = spiceypy.m2eul(cmat, 2, 3, 1)

    except:

        #print('No CK information for {}'.format(time))
        angle3 = 0
        angle2 = 0
        angle1 = 0

    return(np.round(angle3*spiceypy.dpr(),3),
           np.round(angle2*spiceypy.dpr(),3),
           np.round(angle1*spiceypy.dpr(),3))


def check_rotation_matrices():

    # Tests that all defined frames of class 4 and spec=Matrix contains
    # a proper rotation matrix
    frame_class = 4
    frame_ids = spiceypy.kplfrm(frame_class)
    all_matrices_ok = True
    isrot_ntol = 1e-3
    isrot_dtol = 1e-3

    for frame_id in frame_ids:
        frame_name = spiceypy.frmnam(frame_id)

        # Find the TKFRAME SPEC variable, first using frame Id
        frame_spec = ""
        frame_id_spec = "TKFRAME_{}_SPEC".format(frame_id)
        frame_name_spec = "TKFRAME_{}_SPEC".format(frame_name)
        try:
            if spiceypy.dtpool(frame_id_spec):
                frame_spec = frame_id_spec
            else:
                all_matrices_ok = False  # Supposed to be not reachable

        except Exception as e:
            frame_spec = ""

            # Find the TKFRAME SPEC variable, first using frame name
            try:
                if spiceypy.dtpool(frame_name_spec):
                    frame_spec = frame_name_spec
                else:
                    all_matrices_ok = False  # Supposed to be not reachable

            except Exception as e:
                print("Error: " + frame_id_spec + " or " + frame_name_spec +
                      " not found in frame definition: "
                      + str(frame_id) + " - " + frame_name + ", err: " + str(e))

                all_matrices_ok = False

                continue

        try:
            # Check if TKFRAME SPEC is MATRIX
            if str(spiceypy.gcpool(frame_spec, 0, 80)[0]) == 'MATRIX':

                # Find the TKFRAME MATRIX variable
                frame_matrix = "TKFRAME_{}_MATRIX".format(frame_id) \
                                    if frame_spec == frame_id_spec \
                                    else "TKFRAME_{}_MATRIX".format(frame_name)

                try:
                    if spiceypy.dtpool(frame_matrix):

                        try:

                            rot_matrix = spiceypy.gdpool(frame_matrix, 0, 9)
                            rot_matrix = np.asarray(rot_matrix).reshape(3, 3)
                            rot_matrix = np.transpose(rot_matrix)
                            rot_matrix = np.asarray(rot_matrix.tolist())  # To avoid error: strided arrays not supported
                            if spiceypy.isrot(rot_matrix, isrot_ntol, isrot_dtol):
                                print(str(frame_id) + " - " + frame_name + " -> OK!")
                            else:
                                print(str(frame_id) + " - " + frame_name + " -> FAIL!")
                                all_matrices_ok = False

                        except Exception as e:
                            print("Error: Cannot obtain the value from " + frame_matrix
                                  + ", err: " + str(e))

                            all_matrices_ok = False

                            pass
                    else:
                        all_matrices_ok = False  # Supposed to be not reachable

                except Exception as e:
                    print("Error: " + frame_matrix +
                          " not found in frame definition: "
                          + str(frame_id) + " - " + frame_name + ", err: " + str(e))

                    all_matrices_ok = False

                    pass

        except Exception as e:
            print(e)

            all_matrices_ok = False

            pass

    if all_matrices_ok:
        print("ALL ROTATION MATRICES OK!")
    else:
        print("SOME ROTATION MATRICES ARE WRONG!")

    return all_matrices_ok


def check_frame_chain(start_time, end_time, num_samples):

    # Tests that any frame defined can be transformed to J2000
    # from start to end with given interval
    all_frames_ok = True

    et_start = spiceypy.utc2et(start_time)
    et_finish = spiceypy.utc2et(end_time)
    times = np.linspace(et_start, et_finish, num_samples)

    ref_frame = "J2000"
    frame_ids = spiceypy.kplfrm(-1)
    failed_frames = []

    for et in times:
        for frame_id in frame_ids:
            if frame_id not in failed_frames:
                frame_name = spiceypy.frmnam(frame_id)

                try:
                    cmat = spiceypy.pxform(frame_name, ref_frame, et)
                except Exception as e:
                    print("Error: Cannot transform the frame: " + frame_name + " to J2000 " +
                          "at time: " + spiceypy.et2utc(et, 'ISOC', 2) + ", err: " + str(e))
                    all_frames_ok = False
                    failed_frames.append(frame_id)
                    pass

        if len(frame_ids) == len(failed_frames):
            break

    return all_frames_ok
