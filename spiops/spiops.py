import spiceypy as cspice
import math


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

def fov_illum(mk, sensor, time=None, angle='DEGREES'):

    cspice.furnsh(mk)

    angle = angle.upper()

    room = 99
    shapelen = 1000
    framelen = 1000

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

    shape, frame, bsight, n, bounds = cspice.getfov(instid, room, shapelen, framelen)

    rotation = cspice.pxform(frame, 'J2000', time)

    bsight = cspice.mxv(rotation, bsight)

    # The following assumes that the IDs of the given S/C FK have been defined according to the NAIF/ESS
    # standards:
    #
    #    -NXXX
    #
    #       where:
    #          N is the SC id and can consist on a given number of digits
    #          XXX are three digits that identify the sensor

    sc_id = int(str(instid)[:-3])

    ptarg, lt = cspice.spkezp(10, time, 'J2000', 'LT+S', sc_id)

    fov_illumination = cspice.vsep(bsight, ptarg)

    cspice.kclear()

    if angle == 'DEGREES':
        return math.degrees(fov_illumination)
    else:
        return fov_illumination
