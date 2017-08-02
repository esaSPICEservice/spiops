import spiceypy as cspice


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

    if isinstance(time, float):
        time = [time]

    for element in time:

        if format == 'UTC':
            out_elm = cspice.et2utc(element, 'ISOC', 3)

        elif format == 'CAL':
            out_elm = cspice.timout(element,
                                    "YYYY-MM-DDTHR:MN:SC.###::TDB", timlen)
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

        if format == 'UTC':
            out_elm = cspice.utc2et(element)

        elif format == 'CAL':
            out_elm = cspice.str2et(element)
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

    index = 0
    number_of_intervals = range(cspice.wncard(object_cov))
    interval_start = []
    interval_finish = []

    for element in number_of_intervals:
        et_boundaries = cspice.wnfetd(object_cov, index)

        interval_start.append(et_boundaries[0])
        interval_finish.append(et_boundaries[1])

        index += 1

        boundaries = et2cal(et_boundaries, format=time_format)

        interval_start.append(et_boundaries[0])
        interval_finish.append(et_boundaries[1])

        if report and not global_boundary:
            print("Interval {}: {} - {}\n".format(index,
                                                  boundaries[0],
                                                  boundaries[1]))

    #
    # If the global_boundary parameter is set the only output is the global
    # coverage start and finish
    #
    if global_boundary:

        start_time = min(interval_start)
        finish_time = max(interval_finish)

        boundaries = et2cal([start_time, finish_time], format=time_format)

    return boundaries
