#!/usr/bin/env python3

import math
import numpy as np
import spiceypy as spiceypy
import logging
import os
import numpy as np
from spiceypy.utils.support_types import *
from .utils import time
from .utils import plot

from bokeh.io import output_file, show
from bokeh.plotting import figure, output_file, output_notebook, show
from bokeh.models import ColumnDataSource, DatetimeTickFormatter, LabelSet
from math import pi

import spiops
from .utils.time import et_to_datetime

# TODO: change for Bokeh. And try to generalise it. Also put the alternative of
# not using bokeh at some point
import matplotlib.pyplot as plt

"""
The MIT License (MIT)
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

    shape, frame, bsight, n, bounds = spiceypy.getfov(instid, room, shapelen,
                                                    framelen)

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
        ck_ids = spiceypy.ckobj(ck=ck_kernel[0], outCell=ck_ids)

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


def cov_ck_ker(ck, object, support_ker=list(), time_format= 'UTC',
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
    MAXIV = 2000
    WINSIZ = 2 * MAXIV
    MAXOBJ = 10000

    ck_ids = spiceypy.support_types.SPICEINT_CELL(MAXOBJ)
    ck_ids = spiceypy.ckobj(ck, outCell=ck_ids)

    if object_id in ck_ids:

        object_cov = spiceypy.support_types.SPICEDOUBLE_CELL(WINSIZ)
        spiceypy.scard, 0, object_cov
        object_cov = spiceypy.ckcov(ck=ck, idcode=object_id,
                                      needav=False, level='INTERVAL',
                                      tol=0.0, timsys='TDB',
                                      cover=object_cov)

    else:
        print('{} with ID {} is not present in {}.'.format(object,
                                                           object_id, ck))
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

    return (boundaries)



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
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame,  target_frame,et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck1.append(math.degrees(euler[0]))
        eul2_ck1.append(math.degrees(euler[1]))
        eul3_ck1.append(math.degrees(euler[2]))

    spiceypy.unload(ck1)
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


    ck1_filename = ck1.split('/')[-1].split('.')[0]
    ck2_filename = ck2.split('/')[-1].split('.')[0]

    eul1_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul2_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul3_name = '{}_{}'.format(ck1_filename, ck2_filename)

    plot(et_list, [eul1_ck1,eul1_ck2], yaxis_name=['Euler Angle 1 CK1',
                                                   'Euler Angle 1 CK2'],
                                                    title='Euler Angle 1 {}'.format(eul1_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    plot(et_list, [eul2_ck1,eul2_ck2], yaxis_name=['Euler Angle 2 CK1',
                                                   'Euler Angle 2 CK2'],
                                                    title='Euler Angle 2 {}'.format(eul2_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    plot(et_list, [eul3_ck1,eul3_ck2], yaxis_name=['Euler Angle 3 CK1',
                                                   'Euler Angle 3 CK2'],
                                                    title='Euler Angle 3 {}'.format(eul3_name),
                                                    format=plot_style,
                                                    notebook=notebook)

    if report:
        eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=eul1_name)
        eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=eul2_name)
        eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=eul3_name)

    spiceypy.unload(ck2)

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
        plot(et_list, [bsight_ck], yaxis_name=['Degrees', 'Degrees'],
             title='+Z Axis Angle Difference {}'.format(title_name),
             format=plot_style,
             notebook=notebook)


    if report:
        eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=title_name)
        eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=title_name)
        eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=title_name)

    spiceypy.unload(ck2)

    return


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
    angle_ck1 = []
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame,  target_frame,et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck1.append(math.degrees(euler[0]))
        eul2_ck1.append(math.degrees(euler[1]))
        eul3_ck1.append(math.degrees(euler[2]))

        bsight = spiceypy.mxv(rot_mat, boresight)
        bsight_ang = spiceypy.vsep(bsight, boresight)
        bsight_ck1.append(bsight_ang*spiceypy.dpr())

        (rot_axis, rot_angle) = spiceypy.raxisa(rot_mat)
        angle_ck1.append(rot_angle*spiceypy.dpr()*1000)

    spiceypy.unload(ck1)
    spiceypy.furnsh(ck2)

    eul1_ck2 = []
    eul2_ck2 = []
    eul3_ck2 = []
    bsight_ck2 = []
    angle_ck2 = []
    for et in et_list:

        rot_mat = spiceypy.pxform(spacecraft_frame, target_frame, et)
        euler = (spiceypy.m2eul(rot_mat, 1, 2, 3))
        eul1_ck2.append(math.degrees(euler[0]))
        eul2_ck2.append(math.degrees(euler[1]))
        eul3_ck2.append(math.degrees(euler[2]))

        bsight = spiceypy.mxv(rot_mat, boresight)
        bsight_ang = spiceypy.vsep(bsight, boresight)
        bsight_ck2.append(bsight_ang*spiceypy.dpr())

        (rot_axis, rot_angle) = spiceypy.raxisa(rot_mat)
        angle_ck2.append(rot_angle*spiceypy.dpr()*1000)

    eul1_diff = [i - j for i, j in zip(eul1_ck1, eul1_ck2)]
    eul2_diff = [i - j for i, j in zip(eul2_ck1, eul2_ck2)]
    eul3_diff = [i - j for i, j in zip(eul3_ck1, eul3_ck2)]

    bsight_diff = [i - j for i, j in zip(bsight_ck1, bsight_ck2)]

    angle_diff = [abs(i - j) for i, j in zip(angle_ck1, angle_ck2)]

    ck1_filename = ck1.split('/')[-1].split('.')[0]
    ck2_filename = ck2.split('/')[-1].split('.')[0]

    eul1_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul2_name = '{}_{}'.format(ck1_filename, ck2_filename)
    eul3_name = '{}_{}'.format(ck1_filename, ck2_filename)

    if output == 'euler_angles':


        plot(et_list, [eul1_diff, eul2_diff, eul3_diff], yaxis_name=['Degrees','Degrees','Degrees'],
                                                        title='Euler Angle Differences',
                                                        format=plot_style,
                                                        notebook=notebook)

    elif output == 'boresight':
        plot(et_list, bsight_diff, yaxis_name='Degrees',
             title='Boresight Angle Difference',
             format=plot_style,
             notebook=notebook)

    # Attitude Error
    else:
        plot(et_list, angle_diff, yaxis_name='Milidegrees',
             title='Atittude Error',
             format=plot_style,
             notebook=notebook)

    if report:
        eul_angle_report(et_list, eul1_ck1, eul1_ck2, 1, tolerance, name=eul1_name)
        eul_angle_report(et_list, eul2_ck1, eul2_ck2, 2, tolerance, name=eul2_name)
        eul_angle_report(et_list, eul3_ck1, eul3_ck2, 3, tolerance, name=eul3_name)

    spiceypy.unload(ck2)

    return


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

def hga_angles(sc, time):

    if sc == 'MPO':
        hga_zero_frame = sc + '_HGA_APM'
    elif sc == 'MTM':
        return []
    else:
        hga_zero_frame = sc + '_SPACECRAFT'


    hga_el_frame = sc + '_HGA_EL'
    hga_az_frame = sc + '_HGA_AZ'
    hga_frame = sc + '_HGA'

    sc_id = spiceypy.bodn2c(sc)

    try:

        # Get the rotation matrix between two frames
        cmat = spiceypy.pxform(hga_zero_frame, hga_el_frame, time)

        # Get the quaternion elements from the rotation matrix
        quat = spiceypy.m2q(cmat)
        hga_el = np.rad2deg(2. * np.arccos(quat[0]))  # in radians


        cmat = spiceypy.pxform(hga_el_frame, hga_az_frame, time)
        quat = spiceypy.m2q(cmat)
        hga_az = np.rad2deg(2. * np.arccos(quat[0]))  # in radians


        (earth_vec, lt) = spiceypy.spkezp(399, time, hga_frame, 'NONE', sc_id)
        hga_earth = np.rad2deg(spiceypy.vsep([0,0,1], earth_vec))


    except:
        hga_earth = 0
        hga_az = 0
        hga_el = 0

    return([hga_el, hga_az], hga_earth)


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


def solar_array_angle(sa_frame, time):

    sa_zero_frame = sa_frame + '_ZERO'

    try:

        # Get the rotation matrix between two frames
        cmat = spiceypy.pxform(sa_frame, sa_zero_frame, time)

        # Get the quaternion elements from the rotation matrix
        quat = spiceypy.m2q(cmat)
        sa_ang = np.rad2deg(2. * np.arccos(quat[0]))  # in radians

    except:

        #print('No CK information for {}'.format(time))
        sa_ang = 0

    return(sa_ang)


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


def ck_coverage_timeline(metakernel, sc, notebook=True, html_file_name='test',
                         plot_width=975,plot_height=500):

    if notebook:
        output_notebook()
        plot_width = 975
        plot_height = 300
    else:
        output_file(html_file_name + '.html')
        plot_width = plot_width
        plot_height = plot_height

    cov_start = []
    cov_finsh = []
    kernels = []

    with open(metakernel, 'r') as f:
        for line in f:
            if '/ck/' in line and '{}_sc'.format(sc.lower()) in line:
                kernels.append(line.split('/ck/')[-1].strip().split("'")[0])
            if 'PATH_VALUES' in line and '=' in line:
                path = line.split("'")[1] + '/ck/'

    kernels = list(reversed(kernels))
    ck_kernels = []

    for kernel in kernels:

        cov = cov_ck_ker(path+kernel, '{}_SPACECRAFT'.format(sc.upper()), support_ker=metakernel,
                                time_format='TDB')
        if cov:
            cov_start.append(cov[0])
            cov_finsh.append(cov[-1])
            ck_kernels.append(kernel)


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

    p = figure(y_range=ck_kernels, plot_height=plot_height ,plot_width=plot_width,
                title="CK Kernels Coverage", )
    p.hbar(y=ck_kernels, height=0.2, left=start_dt, right=finsh_dt, color="red")

    labels = LabelSet(x='start_dt', y='ck_kernels', text='ck_kernels', level='glyph',
                  x_offset=5, y_offset=-5, source=source)

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


def spk_coverage_timeline(metakernel, sc, notebook=True, html_file_name='test',
                         plot_width=975, plot_height=500):

    if notebook:
        output_notebook()
        plot_width = 975
        plot_height = 300
    else:
        output_file(html_file_name + '.html')
        plot_width = plot_width
        plot_height = plot_height

    cov_start = []
    cov_finsh = []
    kernels = []

    with open(metakernel, 'r') as f:
        for line in f:
            if '/spk/' in line and sc.lower() in line:
                kernels.append(line.split('/spk/')[-1].strip().split("'")[0])
            if 'PATH_VALUES' in line and '=' in line:
                path = line.split("'")[1] + '/spk/'

    kernels = list(reversed(kernels))
    spk_kernels = []

    for kernel in kernels:

        cov = cov_spk_ker(path+kernel, sc.upper(), support_ker=metakernel,
                                time_format='TDB')
        if cov:
            cov_start.append(cov[0][0])
            cov_finsh.append(cov[0][-1])
            spk_kernels.append(kernel)


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

    p = figure(y_range=spk_kernels, plot_height=plot_height ,plot_width=plot_width,
                title="SPK Kernels Coverage", )
    p.hbar(y=spk_kernels, height=0.2, left=start_dt, right=finsh_dt, color="red")

    labels = LabelSet(x='start_dt', y='spk_kernels', text='spk_kernels', level='glyph',
                  x_offset=5, y_offset=-5, source=source)

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




