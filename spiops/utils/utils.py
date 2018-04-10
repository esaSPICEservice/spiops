#from spiops import data as data
from .time import cal2et
from .time import et_to_datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from bokeh.plotting import figure, output_file, output_notebook, show
from bokeh.models import HoverTool
from bokeh.models import DatetimeTickFormatter

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
        
        html_file_name = html_file_name.replace(element,replacement)

    return html_file_name


def convert_ESOCorbit2data(orbit_file, support_ker=''):

    #orbit = data.Data()
    orbit_data = []
    time_list = []
    distance_list = []

    with open(orbit_file, 'r') as f:

        read_data = False
        for line in f:

            if read_data:
                line = line.split()

                time = cal2et(line[0], 'CAL', support_ker=support_ker)
                #TODO: Not including velocity at this point; only distance
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


def plot(time_list, yaxis, yaxis_name='', title='', format='line',
         external_data=[], notebook=False, mission='', target='',
         date_format='TDB'):

    if not isinstance(yaxis_name, list):
        yaxis_name = [yaxis_name]
        yaxis = [yaxis]



    if not title:
        title = '{} {}'.format(mission, yaxis_name).title()

        html_file_name = 'plot_{}_{}_{}-{}.html'.format('Time', yaxis_name,
                                                        mission,
                                                        target)

        html_file_name = valid_url(html_file_name)

    else:
        title = title

        if ' ' in title:
            html_file_name = title.replace(' ', '_').lower()
        else:
            html_file_name = title
        html_file_name = valid_url(html_file_name)

        # TODO: Move this to the time object (convert to datatime)
        # Function needs to be vectorised
        # x = self.time.window
    window_dt = []
    window = time_list
    for element in window:
        window_dt.append(et_to_datetime(element, date_format))

    x = window_dt
    y = yaxis

    if notebook:
        output_notebook()
        plot_width = 975
        plot_height = 500
    else:
        output_file(html_file_name + '.html')
        plot_width = 1000
        plot_height = 1000

    p = figure(title=title,
                plot_width=plot_width,
                plot_height=plot_height,
                x_axis_label='Date in {}'.format(date_format),
                y_axis_label=title,
                x_axis_type="datetime")

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
                 tooltips=[ ('Date', '@x{0.000}'),
                            (title, '@y{0.000}'),
                          ],
                 formatters={'Date': 'numeral',
                             title: 'numeral',
                            })

    p.add_tools(hover)

    if external_data:

        window_dt = []
        window = external_data[0]
        for element in window:
            window_dt.append(et_to_datetime(element, 'TDB'))

        x_ext = window_dt
        y_ext = external_data[1]

        if format == 'circle':
            p.circle(x_ext, y_ext, legend='External Data', size=5, color='red')
        elif format == 'line':
            p.line(x_ext, y_ext, legend='External Data', line_width=2,
                   color='red')

    # add a line renderer with legend and line thickness
    color_list = ['red', 'blue', 'green']
    index = 0
    for element in y:
        if format == 'circle':
            p.circle(x, element, legend=yaxis_name[index],
                     size=2, color=color_list[index])

        elif format == 'line':
            p.line(x, element, legend=yaxis_name[index],
                   line_width=2, color=color_list[index])
        index += 1

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
    theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)

    ax.plot(x, y, z, label= observer.name + ' w.r.t. ' + target +
            ' on ' + observer.trajectory_reference_frame + ' [km]')

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