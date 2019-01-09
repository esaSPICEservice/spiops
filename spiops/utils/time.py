import spiceypy as cspice
from datetime import datetime


#TODO: Function extracted from pyops, need to re-write because of out of date bokeh
def et_to_datetime(et, scale='TDB'):
    """
    convert a SPICE ephemerides epoch (TBD seconds) to a python datetime
    object. The default time scale returned will be TDB but can be set
    to any of the accepted SPICE time scales.

    Args:
        et (float): SPICE ephemerides sceonds (TBD)
        scale (str, optional): time scale of output time (default: TDB)

    Returns:
        datetime: python datetime
    """
    t = cspice.timout(et, 'YYYY-MON-DD HR:MN:SC.### ::{}'.format(scale), 41)
    return datetime.strptime(t, '%Y-%b-%d %H:%M:%S.%f')


def et2cal(time, format='UTC', support_ker=False, unload=False):
    """
    Converts Ephemeris Time (ET) into UTC or Calendar TDB (CAL) time. Accepts
    a single time or a lists of times. This function assumes that the support
    kernels (meta-kernel or leapseconds kernel) has been loaded.

    :param time: Input ET time
    :type time: Union[float, list]
    :param format: Desired output format; 'UTC' or 'CAL'
    :type format: str
    :param unload: If True it will unload the input meta-kernel
    :type unload: bool
    :return: Output time in 'UTC', 'CAL' or 'TDB'
    :rtype: Union[str, list]
    """
    timlen = 62
    out_list = []

    if support_ker:
        cspice.furnsh(support_ker)

    if isinstance(time, float) or isinstance(time, str):
        time = [time]

    for element in time:

        if format == 'UTC':
            out_elm = cspice.et2utc(element, 'ISOC', 3)

        elif format == 'CAL':
            out_elm = cspice.timout(element, "YYYY-MM-DDTHR:MN:SC.###::TDB", timlen)
        else:
            out_elm = element

        out_list.append(out_elm)

    if len(out_list) == 1:
        out_time = out_list[0]
    else:
        out_time = out_list

    if unload:
        cspice.unload(support_ker)

    return out_time


def cal2et(time, format='UTC', support_ker=False, unload=False):
    """
    Converts UTC or Calendar TDB (CAL) time to Ephemeris Time (ET). Accepts
    a single time or a lists of times. This function assumes that the support
    kernels (meta-kernel or leapseconds) kernel has been loaded.

    :param time: Input UTC or CAL time
    :type time: Union[float, list]
    :param format: Input format; 'UTC' or 'CAL'
    :type format: str
    :param unload: If True it will unload the input meta-kernel
    :type unload: bool
    :return: Output ET
    :rtype: Union[str, list]
    """
    out_list = []

    if isinstance(time, str):
        time = [time]

    #
    # We need to specify that is Calendar format in TDB. If it is UTC we need
    # to load the support kernels
    #
    if support_ker:
        cspice.furnsh(support_ker)


    if format == 'CAL':
            time[:] = [x.replace('T', ' ') for x in time]
            time[:] = [x + ' TDB' for x in time]

    for element in time:

        try:
            if format == 'UTC':
                out_elm = cspice.utc2et(element)

            elif format == 'CAL':
                out_elm = cspice.str2et(element)
            else:
                out_elm = element
        except:
                out_elm = cspice.str2et(element)

        out_list.append(out_elm)

    if len(out_list) == 1:
        out_time = out_list[0]
    else:
        out_time = out_list

    if unload:
        cspice.unload(support_ker)

    return out_time


def cov_int(object_cov, object_id, kernel, time_format='TDB',
            global_boundary=False, report=False):
    """
    Generates a list of time windows out of a SPICE cell for which either
    the SPICE API spkcov_c or ckcov_c have been run.


    :param object_cov: SPICE
    :type object_cov:
    :param object_id: Object ID or Name for which we provide the coverage
    :type object_id: Union[str, int]
    :param kernel: Kernel name for which the coverage is being checked
    :type kernel: str
    :param time_format: Desired output format; 'UTC' or 'CAL'
    :type time_format: str
    :param global_boundary: Boolean to indicate whether if we want all the coverage windows or only the absolute start and finish coverage times
    :type global_boundary: bool
    :param report: If True prints the resulting coverage on the screen
    :type report: bool
    :return: Time Windows in the shape of a list
    :rtype: list
    """
    boundaries = False

    if '/' in kernel:
        kernel = kernel.split('/')[-1]

    #
    # Reporting should only be activated if we are not asking for global
    # boundaries.
    #
    if report and not global_boundary:

        try:
            body_name = cspice.bodc2n(object_id)
        except:
            body_name = cspice.frmnam(object_id, 60)

        print("Coverage for {} in {} [{}]:".format(body_name, kernel,
                                                   time_format))

    number_of_intervals = list(range(cspice.wncard(object_cov)))
    interval_start_list = []
    interval_finish_list = []
    coverage = []

    for element in number_of_intervals:
        et_boundaries = cspice.wnfetd(object_cov, element)

        if time_format == 'CAL' or time_format == 'UTC':
            boundaries = et2cal(et_boundaries, format=time_format)
        else:
            boundaries = et_boundaries

        interval_start = boundaries[0]
        interval_finish = boundaries[1]


        if report and not global_boundary:

            print("Interval {}: {} - {}\n".format(element,
                                                  boundaries[0],
                                                  boundaries[1]))

        coverage.append(interval_start)
        coverage.append(interval_finish)
        interval_start_list.append(interval_start)
        interval_finish_list.append(interval_finish)


    #
    # If the global_boundary parameter is set the only output is the global
    # coverage start and finish
    #
    if global_boundary:

        start_time = min(interval_start)
        finish_time = max(interval_finish)

        coverage = et2cal([start_time, finish_time], format=time_format)

    return coverage


def mjd20002et(mjd2000, support_ker=False, unload=False):
    """
    Given a date in MJD2000 (Modified Julian Date 2000) returns the Ephemeris
    time (ET which in SPICE is equivalent to TDB). Accepts a single time entry
    or a list of times.

    :param mjd2000: Date in MJD200
    :type mjd2000: Union[float, list]

    :param support_ker: Support kernels required to run the function. At least it should be a leapseconds kernel (LSK) and optionally a meta-kernel (MK)
    :type support_ker: Union[str, list]
    :param unload: If True it will unload the input support kernel
    :type unload: bool
    :return: Date in ET/TDB
    :rtype: Union[float, list]
    """
    tdb = []

    if support_ker:
        cspice.furnsh(support_ker)

    if not isinstance(mjd2000, list):
        mjd2000 = [mjd2000]

    for time in mjd2000:

        mjd2000 = float(time)
        mjd = mjd2000 + 51544
        jd = mjd + 2400000.5
        jd = str(jd) + ' JD'
        tdb.append(cspice.str2et(jd))

    if unload:
        cspice.unload(support_ker)

    if len(tdb) == 1:
        return tdb[0]
    else:
        return tdb
