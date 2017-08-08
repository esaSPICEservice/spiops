import spiceypy as cspice
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


class Body(object):
    def __init__(self, body, time=object(), target=False):

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

        if item in ['altitude', 'distance']:
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

        tar = self.target
        time = self.time

        for et in time.window:

            #
            # Compute the geometric sub-observer point.
            #
            spoint, trgepc, srfvec = cspice.subpnt(tar.method, tar.name, et,
                                                   tar.frame, time.abcorr,
                                                   self.name)
            subpoint_xyz.append(spoint)

            #
            # Compute the observer's distance from SPOINT.
            #
            dist = cspice.vnorm(srfvec)
            distance.append(dist)

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

            altitude.append(dist - spcrad)


        self.altitude = altitude
        self.distance = distance
        self.subpoint_xyz = subpoint_xyz
        self.subpoint_pgc = subpoint_pgc
        self.subpoint_pcc = subpoint_pcc

        self.geometry_flag = True
        self.previous_tw = self.time.window



    def Plot(self, yaxis, xaxis='time'):

        self.__Geometry()

        x = self.time.window
        y = self.__getattribute__(yaxis)


        fig, ax = plt.subplots()


        ax.plot(x, y, '-')
        ax.set_xlabel('TDB [sec]')
        ax.set_ylabel('{}'.format(yaxis))
        ax.set_title('Spacecraft {}'.format(yaxis))

        plt.show()

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
                                ' on ' + self.trajectory_reference_frame +
                                ' [km]')
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
    def __init__(self, body, time=object(), target=False, frame=False,
                 method='INTERCEPT/ELLIPSOID'):

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
    def __init__(self, body, time=object(), target=False):

        super(Observer, self).__init__(body, time=time, target=target)

        self.frame = '{}_SPACECRAFT'.format(self.name)
