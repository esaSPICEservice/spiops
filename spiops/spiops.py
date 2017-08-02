#!/usr/bin/env python3

import math
import spiceypy as cspice
from spiceypy.utils.support_types import *
from .utils import time

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

def fov_illum(mk, sensor, time=None, angle='DEGREES', abcorr='LT+S',
              report=False, unload=False):
    """
    Determine the Illumination of a given FoV (for light scattering computations
    for example). This function is based on  the following SPICE APIs:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/getfov_c.html
    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezp_c.html


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

    cspice.furnsh(mk)

    if time:
        time = cspice.utc2et(time)
    else:
        time = cspice.utc2et('2016-08-10T00:00:00')

    if angle != 'DEGREES' and angle != 'RADIANS':
        print('angle should be either degrees or radians')

    if isinstance(sensor, str):
        instid = cspice.bodn2c(sensor)
    else:
        instid = sensor

    shape, frame, bsight, n, bounds = cspice.getfov(instid, room, shapelen,
                                                    framelen)

    rotation = cspice.pxform(frame, 'J2000', time)

    bsight = cspice.mxv(rotation, bsight)

    # The following assumes that the IDs of the given S/C FK have been defined
    # according to the NAIF/ESS standards:
    #
    #    -NXXX
    #
    #       where:
    #          N is the SC id and can consist on a given number of digits
    #          XXX are three digits that identify the sensor
    sc_id = int(str(instid)[:-3])

    ptarg, lt = cspice.spkezp(10, time, 'J2000', abcorr, sc_id)

    fov_illumination = cspice.vsep(bsight, ptarg)

    if unload:
        cspice.unload(mk)

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

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html

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
    cspice.furnsh(mk)
    boundaries_list = []
    et_boundaries_list = []

    object_id = cspice.bodn2c(object)
    maxwin = 2000
    spk_count = cspice.ktotal('SPK') - 1

    while spk_count >= 0:

        spk_kernel = cspice.kdata(spk_count, 'SPK', 155, 155, 155)

        spk_ids = cspice.spkobj(spk_kernel[0])

        for id in spk_ids:

            if id == object_id:

                object_cov = SPICEDOUBLE_CELL(maxwin)
                cspice.spkcov(spk_kernel[0], object_id, object_cov)

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
                str(cspice.bodc2n(object_id)), time_format, boundaries_list[0],
                boundaries_list[1]))


    if unload:
        cspice.unload(mk)

    return boundaries_list


def cov_spk_ker(spk, support_ker, object='ALL', time_format= 'TDB',
                report=False, unload=False):
    """
    Provides time coverage summary for a given object for a given SPK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html

    The NAIF utility BRIEF can be used for the same purpose.

    :param spk: SPK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least it should be a leapseconds kernel (LSK) and optionally a meta-kernel (MK)
    :type support_ker: Union[str, list]
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
    cspice.furnsh(spk)
    maxwin = 2000

    cspice.furnsh(support_ker)

    boundaries_list = []
    boundaries = []

    spk_ids = cspice.spkobj(spk)

    if object == 'ALL':
        object_id = spk_ids

    else:
        object_id = cspice.bodn2c(object)

    for id in spk_ids:
        if id == object_id:
            object_cov = SPICEDOUBLE_CELL(maxwin)
            cspice.spkcov(spk, object_id, object_cov)

            boundaries = time.cov_int(object_cov=object_cov,
                                      object_id=object_id,
                                      kernel=spk,
                                      time_format=time_format, report=report)

        boundaries_list += boundaries

        if unload:
            cspice.unload(spk)

    return (boundaries_list)


def cov_ck_obj(mk, object, time_format= 'UTC', global_boundary=False,
               report=False, unload=False):
    """
    Provides time coverage summary for a given object for a list of
    binary CK files provided in a meta-kernel. Several options are
    available. This function is based on the following SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckcov_c.html

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
    cspice.furnsh(mk)
    boundaries_list = []
    et_boundaries_list = []

    object_id = cspice.namfrm(object)
    MAXIV = 2000
    ck_count = cspice.ktotal('CK') - 1
    WINSIZ = 2 * MAXIV
    MAXOBJ = 10000

    while ck_count >= 0:

        ck_ids = cspice.support_types.SPICEINT_CELL(MAXOBJ)
        ck_kernel = cspice.kdata(ck_count, 'CK', 155, 155, 155)
        ck_ids = cspice.ckobj(ck=ck_kernel[0], outCell=ck_ids)

        for id in ck_ids:
            if id == object_id:
                object_cov = cspice.support_types.SPICEDOUBLE_CELL(WINSIZ)
                object_cov = cspice.ckcov(ck=ck_kernel[0], idcode=object_id,
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
                body_name = cspice.bodc2n(object_id)
            except:
                body_name = cspice.frmnam(object_id, 60)

            print("Global Coverage for {} [{}]: {} - {}".format(
                   body_name, time_format, boundaries_list[0],
                boundaries_list[1]))

    if unload:
        cspice.unload(mk)

    return boundaries_list


def cov_ck_ker(ck, support_ker, object='ALL', time_format= 'UTC',
               report=False, unload=False):
    """
    Provides time coverage summary for a given object for a given CK file.
    Several options are available. This function is based on the following
    SPICE API:

    http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckcov_c.html

    The NAIF utility CKBRIEF can be used for the same purpose.

    :param ck: CK file to be used
    :type mk: str
    :param support_ker: Support kernels required to run the function. At least it should be a leapseconds kernel (LSK) and a Spacecraft clock kernel (SCLK) optionally a meta-kernel (MK) which is highly recommended.
    :type support_ker: Union[str, list]
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
    cspice.furnsh(ck)

    if isinstance(support_ker, str):
        support_ker = [support_ker]

    for ker in support_ker:
        cspice.furnsh(ker)

    boundaries_list = []
    boundaries = []

    object_id = cspice.namfrm(object)
    MAXIV = 2000
    WINSIZ = 2 * MAXIV
    MAXOBJ = 10000

    if object == 'ALL':
        ck_ids = cspice.support_types.SPICEINT_CELL(MAXOBJ)
        ck_ids = cspice.ckobj(ck, outCell=ck_ids)
    else:
        ck_ids = [cspice.namfrm(object)]

    for id in ck_ids:
        if id == object_id:
            object_cov = cspice.support_types.SPICEDOUBLE_CELL(WINSIZ)
            cspice.scard, 0, object_cov
            object_cov = cspice.ckcov(ck=ck, idcode=object_id,
                                      needav=False, level='SEGMENT',
                                      tol=0.0, timsys='TDB',
                                      cover=object_cov)

            boundaries = time.cov_int(object_cov=object_cov,
                                      object_id=object_id,
                                      kernel=ck,
                                      time_format=time_format, report=report)

    if boundaries != []:
        boundaries_list += boundaries

    if unload:
        cspice.unload(ck)

    return (boundaries_list)


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
    :type body: str
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

    cspice.furnsh(pck)

    #
    # This can actually be a list of bodies.
    #
    cspice.furnsh(body_spk)

    #
    # Get instantaneous Jupiter state at J2000 and compute instantaneous
    # orbital normal.
    #
    state, lt = cspice.spkezr(body, 0.0, 'J2000', 'NONE', 'SUN')
    normal = cspice.ucrss(state[0:3:1], state[3:6:1])

    #
    # Get J2000 -> IAU_{BODY} rotation at J2000 and compute Body pole
    # direction in J2000 at J2000.
    #
    mat = cspice.pxform('IAU_{}'.format(body), 'J2000', 0.0)
    z = cspice.vpack(0.0, 0.0, 1.0)
    pole = cspice.mxv(mat, z)

    #
    # Compute direction Body orbit's ascending node on Body equator at
    # J2000 in J2000 and print it and Body pole as RA/DEC in J2000 in
    # degrees
    #
    ascnod = cspice.ucrss(pole, normal)
    r, ra, dec = cspice.recrad(pole)

    if report:
        print('POLE RA/DEC = {}/{}'.format(ra*cspice.dpr(), dec*cspice.dpr()))

    r, ra, dec = cspice.recrad(ascnod)

    if report:
        print('ASCNOD RA/DEC = {}/{}'.format(ra * cspice.dpr(), dec * cspice.dpr()))

    #
    # Build two vector from a with POLE as Z and ASNOD as X and print rotation
    # from that frame to J200 as Euler angles.
    #
    mat = cspice.twovec(pole, 3, ascnod, 1)
    matxp = cspice.xpose(mat)
    r3, r2, r1 = cspice.m2eul(matxp, 3, 2, 3)

    if file:
      body_id = cspice.bodn2c(body)
      with open('{}_{}_IF_J2000.tf'.format(mission, body), 'w+') as f:

         f.write(r"\begindata")
         f.write('\n \n')
         f.write('    FRAME_{}_{}_IF_J2000   = {}\n'.format(mission, body,
                                                            frame_id))
         f.write("    FRAME_{}_NAME              = '{}_{}_IF_J2000'\n".format(
                 frame_id, mission, body))
         f.write('    FRAME_{}_CLASS             =  4\n'.format(frame_id))
         f.write('    FRAME_{}_CLASS_ID          = {}\n'.format(frame_id,
                                                              body_id))
         f.write('    FRAME_{}_CENTER            =  {}\n'.format(frame_id,
                                                                 body_id))
         f.write('\n')
         f.write("    TKFRAME_{}_SPEC            = 'ANGLES'\n".format(frame_id))
         f.write("    TKFRAME_{}_RELATIVE        = 'J2000'\n".format(frame_id))
         f.write('    TKFRAME_{}_ANGLES          = (\n'.format(frame_id))
         f.write('                                        {}\n'.format(r3 *
                                                                     cspice.dpr()))
         f.write('                                        {}\n'.format(r2 *
                                                                     cspice.dpr()))
         f.write('                                        {}\n'.format(r1 *
                                                                     cspice.dpr()))
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
            r3 * cspice.dpr(),
            r2 * cspice.dpr(),
            r1 * cspice.dpr())

    if unload:
        cspice.unload(pck)
        cspice.unload(body_spk)

    return

