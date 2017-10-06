import spiceypy as cspice
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from bokeh.plotting import figure, output_file, output_notebook, show
from bokeh.models import HoverTool
from bokeh.models import DatetimeTickFormatter
#from pyops import time as pytime



class Body(object):
    def __init__(self, body, time=object(), target=None):

        if isinstance(body, str):
            name = body
            id = cspice.bodn2c(body)
        else:
            id = body
            name = cspice.bodc2n(body)

        if target:
            self.target = target

        self.name = name
        self.id = id
        self.time = time

        #
        # Parameters for the Geometry Computation
        #
        self.previous_tw = []
        self.geometry_flag = False


    def __getattribute__(self, item):

        if item in ['altitude', 'distance', 'zaxis_target_angle']:
            self.__Geometry()
            return object.__getattribute__(self, item)
        else:
            return object.__getattribute__(self, item)


    def __getattr__(self, item):

        if item in ['state_in_window']:
            self.__StateInWindow()
            return object.__getattribute__(self, item)


    def State(self, target=False, reference_frame=False, current=False):

        if self.target and not target and not reference_frame:
            target = self.target.name
            reference_frame = self.target.frame

        if not self.target:
            if target is False: target = 'J2000'
            if reference_frame is False: reference_frame = 'J2000'

        self.trajectory_reference_frame = reference_frame

        if not current:
            current = self.time.current

        state, lt = cspice.spkezr(target, current, reference_frame,
                                  self.time.abcorr, self.name)

        return state


    def Orientation(self, frame='', target_frame='', current=False,
                    format='msop quaternions'):

        if self.target and not target_frame:
            target_frame = self.target.frame

        if not self.target and not target_frame:
            target_frame = 'J2000'

        if not frame:
            frame = self.frame

        if not current:
            current = self.time.current
        else:
            #TODO: need to add a time conversion here
            current = current

        rot_mat = cspice.pxform(target_frame, frame, current)

        if format == 'spice quaternions':
            orientation = cspice.m2q(rot_mat)

        if format == 'msop quaternions':
            quaternions = cspice.m2q(rot_mat)
            orientation = [-quaternions[1],
                           -quaternions[2],
                           -quaternions[3],
                            quaternions[0]]

        elif format == 'euler angles':
            orientation = (cspice.m2eul(rot_mat, 3, 2, 1))

        elif format == 'rotation matrix':
            orientation = rot_mat

        return orientation


    def __StateInWindow(self, target=False, reference_frame=False, start=False,
                     finish=False):

        state_in_window = []
        for et in self.time.window:
            state_in_window.append(self.State(target, reference_frame, et))

        self.state_in_window = state_in_window

        return


    def __Geometry(self):

        #if self.geometry_flag is True and \
        #                self.time.window.all() == self.previous_tw.all():
        #    return

        distance = []
        altitude = []
        subpoint_xyz = []
        subpoint_pgc = []
        subpoint_pcc = []
        zaxis_target_angle = []

        tar = self.target
        time = self.time

        for et in time.window:

            #
            # Compute the distance
            #
            ptarg, lt = cspice.spkpos(tar.name, et, tar.frame, time.abcorr,
                                      self.name)
            vout, vmag = cspice.unorm(ptarg)
            distance.append(vmag)

            #
            # Compute the geometric sub-observer point.
            #
            spoint, trgepc, srfvec = cspice.subpnt(tar.method, tar.name, et,
                                                   tar.frame, time.abcorr,
                                                   self.name)
            subpoint_xyz.append(spoint)

            #
            # Compute the observer's altitude from SPOINT.
            #
            dist = cspice.vnorm(srfvec)
            altitude.append(dist)


            #
            # Convert the sub-observer point's rectangular coordinates to
            # planetographic longitude, latitude and altitude.
            #
            spglon, spglat, spgalt = cspice.recpgr(tar.name, spoint,
                                                   tar.radii_equ, tar.flat)

            #
            # Convert radians to degrees.
            #
            spglon *= cspice.dpr()
            spglat *= cspice.dpr()

            subpoint_pgc.append([spglon, spglat, spgalt])

            #
            #  Convert sub-observer point's rectangular coordinates to
            #  planetocentric radius, longitude, and latitude.
            #
            spcrad, spclon, spclat = cspice.reclat(spoint)

            #
            # Convert radians to degrees.
            #
            spclon *= cspice.dpr()
            spclat *= cspice.dpr()

            subpoint_pcc.append([spcrad, spclon, spclat])

            #
            # Compute the angle between the observer's Z axis and the geometric
            # sub-observer point
            #
            obs_tar, ltime = cspice.spkpos(tar.name, et,
                                                   'J2000', time.abcorr,
                                                   self.name)
            obs_zaxis = [0,0,1]
            matrix = cspice.pxform(self.frame, 'J2000', et)
            vecout = cspice.mxv(matrix, obs_zaxis)

            zax_target_angle = cspice.vsep(vecout, obs_tar)
            zax_target_angle *= cspice.dpr()
            zaxis_target_angle.append(zax_target_angle)


        self.distance = distance
        self.altitude = altitude
        self.subpoint_xyz = subpoint_xyz
        self.subpoint_pgc = subpoint_pgc
        self.subpoint_pcc = subpoint_pcc
        self.zaxis_target_angle = zaxis_target_angle

        self.geometry_flag = True
        self.previous_tw = self.time.window



    def Plot(self, yaxis, xaxis='time', title='', external_data=[],
             notebook=False):
        import spiops.utils.utils as utils
        import spiops.utils.time as utime

        self.__Geometry()

        if not title:
            title ='{} {}'.format(self.name, yaxis).title()

            html_file_name = 'plot_{}_{}_{}-{}.html'.format(xaxis, yaxis,
                                                            self.name,
                                                            self.target.name)

            html_file_name = utils.valid_url(html_file_name)

        else:
            title=title

            html_file_name = title
            html_file_name = utils.valid_url(html_file_name)

        #TODO: Move this to the time object (convert to datatime)
        # Function needs to be vectorised
        #x = self.time.window
        window_dt = []
        window = self.time.window
        for element in window:
            window_dt.append(utime.et_to_datetime(element, 'TDB'))

        x = window_dt
        y = self.__getattribute__(yaxis)

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
                   x_axis_label='Date in TBD',
                   y_axis_label='{}'.format(yaxis).title(),
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
                               ('{}'.format(yaxis).title(), '@y{0.000}'),
                             ],
                    formatters={'Date': 'numeral',
                                '{}'.format(yaxis).title(): 'numeral',
                               })

        p.add_tools(hover)

        if external_data:

            window_dt = []
            window = external_data[0]
            for element in window:
                window_dt.append(utime.et_to_datetime(element, 'TDB'))

            x_ext = window_dt

            #x_ext = external_data[0]
            y_ext = external_data[1]
            p.circle(x_ext, y_ext, legend='External Data', size=5, color='red')

        # add a line renderer with legend and line thickness
        p.line(x, y, legend='{}'.format(yaxis).title(), line_width=2)

        # show the results
        show(p)


        return



    def Plot3D(self, data='trajectory', reference_frame=False):


        #TODO: Arrange the reference frame flow
        if not self.state_in_window:
            self.__StateInWindow(reference_frame=reference_frame)


        data = self.state_in_window

        x, y, z, = [], [], []

        for element in data:
           x.append(element[0])
           y.append(element[1])
           z.append(element[2])


        mpl.rcParams['legend.fontsize'] = 10

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)

        ax.plot(x, y, z, label= self.name + ' w.r.t. ' + self.target.name +
                ' on ' + self.trajectory_reference_frame + ' [km]')

        ax.legend()

        # Make data
        u = np.linspace(0, 2 * np.pi, 360)
        v = np.linspace(0, np.pi, 360)
        x = self.target.radii[0] * np.outer(np.cos(u), np.sin(v))
        y = self.target.radii[1] * np.outer(np.sin(u), np.sin(v))
        z = self.target.radii[2] * np.outer(np.ones(np.size(u)), np.cos(v))

        # Plot the surface
        ax.plot_surface(x, y, z, color='r')

        plt.show()

        return


class Target(Body):

    def __init__(self, body, time=object(),
                 target=False, frame='',
                 method='INTERCEPT/ELLIPSOID'):
        """

        :param body:
        :type body:
        :param time:
        :type time:
        :param target: It no target is provided the default is 'SUN'
        :type target:
        :param frame:
        :type frame:
        :param method:
        :type method:
        """

        #
        # In the target parameter for the following if statement we need to use
        # the class object(), which is empty, to avoid A Recursion Error.
        #
        if not target:
            target = Target('SUN', time=time, target=object())

        super(Target, self).__init__(body, time=time, target=target)

        if not frame:
            self.frame = 'IAU_{}'.format(self.name)
        else:
            self.frame = frame

        self.method = method

        self.__getRadii()

    def __getRadii(self):
        self.radii = cspice.bodvar(self.id, 'RADII', 3)

        self.radii_equ = self.radii[0]
        self.radii_pol = self.radii[2]
        self.flat = (self.radii_equ - self.radii_pol) / self.radii_equ


class Observer(Body):
    def __init__(self, body, time=object(), target=False, frame=''):

        super(Observer, self).__init__(body, time=time, target=target)

        if not frame:
            self.frame = '{}_SPACECRAFT'.format(self.name)
            if cspice.namfrm(self.frame) == 0:
                self.frame = self.name
            if cspice.namfrm(self.frame) == 0:
                print('The frame name has not been able to be built; please introduce it manually')
        else:
            self.frame = frame


