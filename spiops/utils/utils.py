from spiops.utils.time import cal2et
from spiops.utils.time import et_to_datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiceypy
import numpy as np
from bokeh.plotting import figure, output_file, output_notebook, show
from bokeh.models import HoverTool
from bokeh.models import ColumnDataSource
from bokeh.models import DatetimeTickFormatter
from bokeh.models import LabelSet
from bokeh.models import Range1d
from tempfile import mkstemp
from shutil import move
import os
import glob
import platform
from os import fdopen, chmod
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from spiops.data import images

MISSION_SPACECRAFTS = {
    'ROSETTA': ["ROS"],
    'VENUS-EXPRESS': ["VEX"],
    'MARS-EXPRESS': ["MEX"],
    'ExoMars2016': ["TGO", "EDM"],
    'BEPICOLOMBO': ["MMO", "MPO", "MTM"],
    'JUICE': ["JUICE"],
    'SOLAR-ORBITER': ["SOLO"],
    'ExoMarsRSP': ["SP", "RM"]
}

def valid_url(html_file_name):
    """
    This function returns a valid URL for an HTML given a filename. 
    The filename is checked in such way that URL non valid characters
    are replaced by other characters. 
    
    This was used due to the fact that we were using the following string:
    
       '67P/CG' -> '67P-CG'
    
    as part of an URL for 2D plotting and it would not work
    
    :param html_file_name: Input filename 
    :type html_file_name: str
    :return: Corrected Input filename without URL non valid characters
    :rtype: str
    """""

    for element in ['$', '_', '.', '+', '!', '*', '(', ')', '/', '\\']:
        
        if element in ['/', '\\', '!', '*', '$']:
            replacement = '-'
        else:
            replacement = '_'
        
        html_file_name = html_file_name.replace(element, replacement)

    return html_file_name


def get_exe_dir():
    if platform.system() == 'Darwin':
        if platform.machine() == 'x86_64':
            return '/exe/macintel_osx_64bit'
        elif platform.machine() == 'arm64':
            return '/exe/macintel_m1_64bit'
    else:
        return '/exe/pc_linux_64bit'

    raise Exception("Unsupported platform! Unable to determine exe directory.")


def get_skd_path(kernel_path):
    if "former_versions" in kernel_path:
        return '/'.join(kernel_path.split('/')[:-3])
    else:
        return '/'.join(kernel_path.split('/')[:-2])


def convert_ESOCorbit2data(orbit_file, support_ker=''):
    time_list = []
    distance_list = []

    with open(orbit_file, 'r') as f:

        read_data = False
        for line in f:

            if read_data:
                line = line.split()

                time = cal2et(line[0], 'CAL', support_ker=support_ker)
                # TODO: Not including velocity at this point; only distance
                distance = np.sqrt(float(line[1])*float(line[1]) +
                                                 float(line[2])*float(line[2]) +
                                                 float(line[3])*float(line[3]))

                time_list.append(time)
                distance_list.append(distance)

            if 'META_STOP' in line:
                read_data = True

            if 'META_START' in line:
                read_data = False

    return [time_list, distance_list]


def convert_OEM2data():
    return


def plot(xaxis, yaxis, xaxis_name='Date', yaxis_name='', title='', format='line',
         external_data=[], notebook=False, mission='', target='', yaxis_units='',
         date_format='TDB', plot_width=975, plot_height=300,
         fill_color=[], fill_alpha=0, background_image=False,
         line_width=2):

    if not isinstance(yaxis_name, list):
        yaxis_name = [yaxis_name]
        yaxis = [yaxis]

    if not title:
        title = '{} {}'.format(mission, yaxis_name).title().upper()

        html_file_name = 'plot_{}_{}_{}-{}.html'.format('Time', yaxis_name,
                                                        mission,
                                                        target)

        html_file_name = valid_url(html_file_name)

    else:
        title = title.upper()

        if ' ' in title:
            html_file_name = title.replace(' ', '_').upper()
        else:
            html_file_name = title
        html_file_name = valid_url(html_file_name)

        # TODO: Move this to the time object (convert to datatime)
        # Function needs to be vectorised
        # x = self.time.window

    if xaxis_name == 'Date':
        window_dt = []
        window = xaxis
        for element in window:
            window_dt.append(et_to_datetime(element, date_format))

        x = window_dt
    else:
        x = xaxis

    y = yaxis

    if notebook:
        output_notebook()
    else:
        output_file(html_file_name + '.html')

    if xaxis_name == 'Date':
        x_axis_type = "datetime"
    else:
        x_axis_type = "auto"

    p = figure(title=title,
               plot_width=plot_width,
               plot_height=plot_height,
               x_axis_label=xaxis_name.upper(),
               y_axis_label=yaxis_units,
               x_axis_type=x_axis_type)

    if xaxis_name == 'Date':

        p.xaxis.formatter = DatetimeTickFormatter(
                                                    seconds=["%Y-%m-%d %H:%M:%S"],
                                                    minsec=["%Y-%m-%d %H:%M:%S"],
                                                    minutes=["%Y-%m-%d %H:%M:%S"],
                                                    hourmin=["%Y-%m-%d %H:%M:%S"],
                                                    hours=["%Y-%m-%d %H:%M:%S"],
                                                    days=["%Y-%m-%d %H:%M:%S"],
                                                    months=["%Y-%m-%d %H:%M:%S"],
                                                    years=["%Y-%m-%d %H:%M:%S"],
                                                )

    hover = HoverTool(
                 tooltips=[(xaxis_name, '@x{0.000}'),
                           (title, '@y{0.000}')],
                 formatters={xaxis_name: 'numeral',
                             title: 'numeral'})
    p.add_tools(hover)

    if external_data:

        window_dt = []
        window = external_data[0]
        for element in window:
            window_dt.append(et_to_datetime(element, 'TDB'))

        x_ext = window_dt
        y_ext = external_data[1]

        if format == 'circle':
            p.circle(x_ext, y_ext, size=5, color='red')
        elif format == 'line':
            p.line(x_ext, y_ext, line_width=2,
                   color='red')

    # add a line renderer with legend and line thickness
    color_list = ['red', 'green', 'blue', 'orange', "black", 'darkgoldenrod', 'chocolate', 'aqua', 'coral',
                  'darkcyan', 'cornflowerblue' 'aquamarine', 'darkturquoise', 'cornsilk']
    index = 0
    color_idx = 0

    if background_image:
        if 'TGO' in mission.upper() or 'MEX' in mission.upper():
            image = 'Mars_Viking_MDIM21_ClrMosaic_global_1024.jpg'
        else:
            image = 'Earth_Contemporary_Basic.png'
        p.image_url(url=[os.path.join(os.path.dirname(images.__file__), image)], x=-180, y=-90,
                    w=360, h=180, anchor="bottom_left", global_alpha=0.6)
        left, right, bottom, top = -180, 180, -90, 90
        p.x_range = Range1d(left, right)
        p.y_range = Range1d(bottom, top)

    if format == 'scatter':
        is_multi_legend = len(y) <= len(yaxis_name)
        for idx in range(len(y)):
            p.circle([x[idx]], [y[idx]], size=3,
                     color=color_list[color_idx] if is_multi_legend else color_list[0],
                     legend=str(yaxis_name[idx]).upper() if is_multi_legend else str(yaxis_name[0]).upper())
            color_idx = idx % len(color_list)
    else:
        for element in y:

            legend = str(yaxis_name[index]).upper()

            if format == 'circle':
                p.line(x, element, line_width=line_width, color=color_list[color_idx], legend=legend)
                p.circle(x, element, fill_color="white", size=8)

            if format == 'circle_only':
                p.circle(x, element, size=3, color=color_list[color_idx], legend=legend)

            elif format == 'line':
                p.line(x, element, line_width=line_width, color=color_list[color_idx], legend=legend)

            index += 1
            color_idx = index % len(color_list)

    p.legend.click_policy = "hide"

    # show the results
    show(p)

    return


def plot3d(data, observer, target):

    x, y, z, = [], [], []

    for element in data:
        x.append(element[0])
        y.append(element[1])
        z.append(element[2])

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)

    label = observer.name + ' w.r.t. ' + target.name + ' on ' + observer.trajectory_reference_frame + ' [km]'
    ax.plot(x, y, z, label=label)
    ax.legend()

    # Make data
    u = np.linspace(0, 2 * np.pi, 360)
    v = np.linspace(0, np.pi, 360)
    x = target.radii[0] * np.outer(np.cos(u), np.sin(v))
    y = target.radii[1] * np.outer(np.sin(u), np.sin(v))
    z = target.radii[2] * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z, color='r')

    plt.show()

    return


def plot_attitude_error(error, max_ang_error, title, plot_style, notebook):

    print('Avg QX error: ', np.mean(error[:, 1]))
    print('Avg QY error: ', np.mean(error[:, 2]))
    print('Avg QZ error: ', np.mean(error[:, 3]))
    print('Avg QW error: ', np.mean(error[:, 4]))
    print('Max angular error [mdeg]: ' + str(max_ang_error))

    plot(error[:, 0],
         [error[:, 1], error[:, 2], error[:, 3], error[:, 4]],
         yaxis_name=['QX', 'QY', 'QZ', 'QW'],
         title=title,
         format=plot_style,
         yaxis_units='Q [-]',
         notebook=notebook)


def replace(file_path, pattern, subst):

    replaced = False

    # Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:

                updated_line = line.replace(pattern, subst)
                new_file.write(updated_line)
                # flag for replacing having happened
                if updated_line != line:
                    replaced = True

    if replaced:
        # Update the permissions
        chmod(abs_path, 0o644)
        # Move new file
        if file_path.isupper():
            move(abs_path, file_path.split('.')[0]+'_LOCAL.TM')
        else:
            move(abs_path, file_path.split('.')[0] + '_local.tm')

        return True

    return False


def get_latest_kernel(kernel_type, path, pattern, dates=False,
                      excluded_kernels=False):

    kernels = []
    kernel_path = os.path.join(path, kernel_type)

    #
    # Get the kernels of type ``type`` from the ``path``/``type`` directory.
    #
    kernels_with_path = glob.glob(kernel_path + '/' + pattern)

    #
    # Include kernels in former_versions if the directory exists except for
    # meta-kernel generation
    #
    if os.path.isdir(kernel_path + '/former_versions'):
        kernels_with_path += glob.glob(kernel_path + '/former_versions/' + pattern)

    for kernel in kernels_with_path:
        kernels.append(kernel.split('/')[-1])

    #
    # Put the kernels in order
    #
    kernels.sort()

    #
    # We remove the kernel if it is included in the excluded kernels list
    #
    if excluded_kernels:
        for kernel in excluded_kernels:
            if kernel in kernels:
                kernels.remove(kernel)

    if not dates:
        #
        # Return the latest kernel
        #
        return kernels.pop()
    else:
        #
        # Return all the kernels with a given date
        #
        previous_kernel = ''
        kernels_date = []
        for kernel in kernels:
            if previous_kernel and previous_kernel.upper().split('_V')[0] == kernel.upper().split('_V')[0]:
                kernels_date.remove(previous_kernel)

            previous_kernel = kernel
            kernels_date.append(kernel)

        return kernels_date


def get_sc(kernel):
    if 'ROSETTA' in kernel.upper():
        return 'ROS'
    if 'VENUS-EXPRESS' in kernel.upper():
        return 'VEX'
    if 'MARS-EXPRESS' in kernel.upper():
        return 'MEX'
    if 'EXOMARS2016' in kernel.upper():
        if 'edm' in kernel:
            return 'em16_edm'
        else:
            return 'em16_tgo'
    if 'BEPICOLOMBO' in kernel.upper():
        if 'mmo' in kernel:
            return 'bc_mmo'
        else:
            return 'bc_mpo'
    if 'JUICE' in kernel.upper():
        return 'juice'
    if 'SOLAR-ORBITER' in kernel.upper():
        return 'solo'
    if 'EXOMARSRSP' in kernel.upper():
        if '_sp_' in kernel:
            return 'emrsp_sp'
        else:
            return 'emrsp_rm'


def get_mission(sc):
    for mission in MISSION_SPACECRAFTS:
        if sc.upper() in MISSION_SPACECRAFTS[mission]:
            return mission
    return None


def target2frame(target):

    if target == '67P/C-G':
        target_frame = '67P/C-G_CK'
    elif target == '21 LUTETIA':
        target_frame = 'LUTETIA_FIXED'
    elif target == 'DIDYMOS' or target == 'DIDYMOON':
        target_frame = '{}_FIXED'.format(target)
    else:
        try:
            target_frame = 'IAU_' + target.upper()
            target_frame_id = spiceypy.namfrm(target_frame)
            if target_frame_id == 0:
                raise Exception
        except:
            try:
                target_id = str(spiceypy.bodn2c(target))
                target_frame = spiceypy.frmnam(int(target_id))
            except:
                target_id += '000'
                target_frame = spiceypy.frmnam(int(target_id))

    return target_frame


def findIntersection(x1, y1, x2, y2, x3, y3, x4, y4):
    px = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    py = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    return [px, py]


def findNearest(array, value):
    """
    Determine for a given value, the element of the array which is closest to
    this value.

    @param array: Input N-Dimensional Array with values
    @type array: Numpy array
    @param value: Value that we want to match with one value of the array
    @type value: float
    @return: Index and value in array closer to value
    @rtype: list
    """

    array = np.asarray(array)
    idx = np.unravel_index(np.argmin(np.abs(array - value)), array.shape)

    return idx, array[idx]


def get_ck_kernel_color(kernel, frame):
    color = "lawngreen"
    ck_type = 'xxx'

    if 'MPO' in frame \
            or 'MMO' in frame \
            or 'MTM' in frame \
            or 'TGO' in frame:
        ck_type = kernel.split('_')[3]

    elif 'JUICE' in frame:
        ck_type = kernel.split('_')[2]

    if ck_type[2] == 'p' \
            or ck_type == 'attc':
        color = 'orange'

    elif ck_type[2] == 'r':
        color = 'green'

    elif ck_type[2] == 't' \
            or ck_type == 'crema':
        color = 'red'

    elif ck_type[2] == 'c':
        color = 'purple'

    elif ck_type[2] == 'm' \
            or ck_type == 'meas':
        color = 'blue'

    return color


def get_plot_style(plot_height, num_rows):
    if plot_height is None:
        empty_plot_height = 50  # px margin for the plot title and X axis scale and labels
        row_height = 40  # px per kernel row
        plot_height = (num_rows * row_height) + empty_plot_height

    hbar_height = 0.15
    lbl_y_offset = int((plot_height / num_rows) * (hbar_height * 0.75))

    return plot_height, hbar_height, lbl_y_offset


def prepare_coverage_plot(p, source_dict, x, y, lbl_y_offset):
    text_font_size = "10pt"

    source = ColumnDataSource(data=source_dict)
    labels = LabelSet(x=x, y=y, text=y, level='glyph', x_offset=-2, y_offset=lbl_y_offset,
                      source=source, render_mode='canvas', text_font_size=text_font_size)
    p.add_layout(labels)

    p.xaxis.formatter = DatetimeTickFormatter(seconds=["%Y-%m-%d %H:%M:%S"],
                                              minsec=["%Y-%m-%d %H:%M:%S"],
                                              minutes=["%Y-%m-%d %H:%M:%S"],
                                              hourmin=["%Y-%m-%d %H:%M:%S"],
                                              hours=["%Y-%m-%d %H:%M:%S"],
                                              days=["%Y-%m-%d %H:%M:%S"],
                                              months=["%Y-%m-%d %H:%M:%S"],
                                              years=["%Y-%m-%d %H:%M:%S"])

    p.xaxis.major_label_orientation = 0  # pi/4
    p.yaxis.visible = False
    p.xaxis.axis_label_text_font_size = text_font_size


def get_exclude_intervals(mission_config=None, method="all"):
    exclude_intervals_str = []
    if mission_config is not None:
        if "exclude_intervals" in mission_config:
            if method in mission_config["exclude_intervals"]:
                exclude_intervals_str = mission_config["exclude_intervals"][method]
            elif "all" in mission_config["exclude_intervals"]:
                exclude_intervals_str = mission_config["exclude_intervals"]["all"]

    exclude_intervals = []
    for ex_int_str in exclude_intervals_str:
        exclude_intervals.append([spiceypy.str2et(ex_int_str[0]), spiceypy.str2et(ex_int_str[1])])

    return exclude_intervals


def is_excluded(et, exclude_intervals, curr_interval_idx=-1):
    if len(exclude_intervals):

        if curr_interval_idx < 0:
            curr_interval_idx = 0

        if curr_interval_idx < len(exclude_intervals):
            curr_int = exclude_intervals[curr_interval_idx]
            if et < curr_int[0]:
                # Pending to enter in current interval
                return False, curr_interval_idx
            elif et <= curr_int[1]:
                # Inside current interval
                return True, curr_interval_idx
            else:
                # Passed current interval, goto next
                return False, curr_interval_idx + 1
        else:
            return False, len(exclude_intervals)

    return False, -1