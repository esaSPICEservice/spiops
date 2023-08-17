import spiceypy as spiceypy
import numpy as np
from spiops import spiops
from spiops.utils import utils
from spiops.utils import naif
from bokeh.layouts import row
from bokeh.plotting import figure, show, output_notebook
from bokeh.models.glyphs import Ellipse
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid


class Body(object):
    def __init__(self, body, time=object(), target=None, mission_config=None):

        if isinstance(body, str):
            name = body
            id = spiceypy.bodn2c(body)
        else:
            id = body
            name = spiceypy.bodc2n(body)

        if target:
            self.target = target

        if mission_config:
            self.mission_config = mission_config

        self.name = name
        self.id = id
        self.time = time

        #
        # Parameters for the Geometry Computation
        #
        self.previous_tw = []
        self.geometry_flag = False


    def __getattribute__(self, item, boresight=''):

        if item in ['altitude',
                    'distance',
                    'velocity',
                    'zaxis_target_angle',
                    'myaxis_target_angle',
                    'groundtrack',
                    'boresight_groundtrack'
                    'trajectory']:
            self.__Geometry(boresight=boresight)
            return object.__getattribute__(self, item)
        elif item in ['sa_ang_p',
                      'sa_ang_n',
                      'sa_ang',
                      'saa_sa',
                      'saa_sc',
                      'hga_earth',
                      'hga_angles',
                      'mga_earth',
                      'mga_angles',
                      'roll_angles']:
            self.__Structures()
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

        state, lt = spiceypy.spkezr(target, current, reference_frame,
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

        rot_mat = spiceypy.pxform(target_frame, frame, current)

        if format == 'spice quaternions':
            orientation = spiceypy.m2q(rot_mat)

        if format == 'msop quaternions':
            quaternions = spiceypy.m2q(rot_mat)
            orientation = [-quaternions[1],
                           -quaternions[2],
                           -quaternions[3],
                            quaternions[0]]

        elif format == 'euler angles':
            orientation = (spiceypy.m2eul(rot_mat, 3, 2, 1))

        elif format == 'rotation matrix':
            orientation = rot_mat

        return orientation

    def __ClockDrift(self, enddate=False):

        try:

            skd_path = None
            sclk_start_coeff_idx = 1
            if self.mission_config is not None:
                skd_path = self.mission_config["skd_path"]
                if "sclk_start_coeff_idx" in self.mission_config:
                    sclk_start_coeff_idx = self.mission_config["sclk_start_coeff_idx"][self.name.lower()]

            sclk_path = naif.get_latest_step_sclk(self.name, skd_path=skd_path)

            coeffs = naif.read_sclk_coefficiends(sclk_path)

            sclk_start = int(coeffs[sclk_start_coeff_idx][0])
            sclk_end = int(coeffs[-1][0])

        except Exception as ex:
            print("Error: Could not obtain SCLK time bounds. Error: " + str(ex))
            return

        ticks_per_second = spiceypy.gdpool('SCLK01_MODULI_{}'.format(str(-1*self.id)), 0, 1000)[1]
        step = int(((sclk_end - sclk_start)/10000))  # 10000 points in the plot

        if not enddate:
            et_end = self.time.getTime('finish', 'utc')
        else:
            et_end = spiceypy.utc2et(enddate)

        sclk = []
        ephtime = []
        for clk in range(sclk_start, sclk_end, step):
            sclk.append(clk)
            et = spiceypy.sct2e(self.id, clk)
            ephtime.append(et)

        dates = []
        drift = []
        for j in range(0, len(ephtime), 1):
            if ephtime[j] >= et_end:
                break
            drift.append((sclk[j] - sclk[0]) / ticks_per_second - (ephtime[j] - ephtime[0]))
            dates.append(ephtime[j])

        self.clock_dates = dates
        self.clock_drift = drift

        return


    def __StateInWindow(self, target=False, reference_frame=False, start=False,
                     finish=False):

        state_in_window = []
        for et in self.time.window:
            state_in_window.append(self.State(target, reference_frame, et))

        self.state_in_window = state_in_window

        return


    def __Structures(self):

        if self.structures_flag is True and \
                self.time.window.all() == self.previous_tw.all():
            return

        time = self.time

        #
        # We determine the mission
        #
        if self.name == 'TGO':
            plus_array = 'TGO_SA+Z'
            minus_array = 'TGO_SA-Z'
        elif self.name == 'MPO':
            plus_array = 'MPO_SA'
            minus_array = ''
        elif self.name == 'MTM':
            plus_array = 'MTM_SA+X'
            minus_array = 'MTM_SA-X'
        else:
            plus_array = '{}_SA+Y'.format(self.name.upper())
            minus_array = '{}_SA-Y'.format(self.name.upper())

        #
        # Solar Arrays
        #
        sa_ang1_p_list = []
        sa_ang1_n_list = []
        sa_ang2_p_list = []
        sa_ang2_n_list = []
        sa_ang3_p_list = []
        sa_ang3_n_list = []

        saa_sa_p_list = []
        saa_sa_n_list = []

        saa_sc_x_list = []
        saa_sc_y_list = []
        saa_sc_z_list = []

        hga_earth = []
        hga_angles_el = []
        hga_angles_az = []

        mga_earth = []
        mga_angles_el = []
        mga_angles_az = []

        roll_angle_1 = []
        roll_angle_2 = []
        roll_angle_3 = []

        #
        # High Gain Antennas
        #
        for et in time.window:

            #
            # SA mechanisms
            #
            try:
                # Of course we need to include all possible cases including only one Solar Array

                (sa_ang1_p, sa_ang2_p, sa_ang3_p) = spiops.solar_array_angles(plus_array, et)
                if minus_array:
                    (sa_ang1_n, sa_ang2_n, sa_ang3_n) = spiops.solar_array_angles(minus_array, et)
                saa = spiops.solar_aspect_angles(self.name, et)

                sa_ang1_p_list.append(sa_ang1_p)
                sa_ang2_p_list.append(sa_ang2_p)
                sa_ang3_p_list.append(sa_ang3_p)
                if minus_array:
                    sa_ang1_n_list.append(sa_ang1_n)
                    sa_ang2_n_list.append(sa_ang2_n)
                    sa_ang3_n_list.append(sa_ang3_n)

                saa_sa_p_list.append(saa[0][0])
                if minus_array:
                    saa_sa_n_list.append(saa[0][1])

                saa_sc_x_list.append(saa[1][0])
                saa_sc_y_list.append(saa[1][1])
                saa_sc_z_list.append(saa[1][2])

            except:
                sa_ang1_p_list.append(0)
                sa_ang2_p_list.append(0)
                sa_ang3_p_list.append(0)
                if minus_array:
                    sa_ang1_n_list.append(0)
                    sa_ang2_n_list.append(0)
                    sa_ang3_n_list.append(0)

                saa_sa_p_list.append(0)
                if minus_array:
                    saa_sa_n_list.append(0)

                saa_sc_x_list.append(0)
                saa_sc_y_list.append(0)
                saa_sc_z_list.append(0)

            #
            # HGA mechanisms
            #
            if self.name != 'MTM' and self.name != 'JUICE':
                try:
                    hga_angles_ang, hga_earth_ang = spiops.hga_angles(self.name, et)
                except:
                    hga_angles_ang = [0, 0]
                    hga_earth_ang = 0

                hga_angles_el.append(hga_angles_ang[1])
                hga_angles_az.append(hga_angles_ang[0])
                hga_earth.append(hga_earth_ang)

            #
            # MGA mechanisms
            #
            if self.name != 'MTM':
                try:
                    mga_angles_ang, mga_earth_ang = spiops.mga_angles(self.name, et)
                except:
                    mga_angles_ang = [0, 0]
                    mga_earth_ang = 0

                mga_angles_el.append(mga_angles_ang[1])
                mga_angles_az.append(mga_angles_ang[0])
                mga_earth.append(mga_earth_ang)

            #
            # Roll angle
            #
            try:
                roll_angle = spiops.roll(et)
            except:
                roll_angle = [0, 0, 0]

            roll_angle_1.append(roll_angle[0])
            roll_angle_2.append(roll_angle[1])
            roll_angle_3.append(roll_angle[2])

        self.sa_ang_p = [sa_ang1_p_list, sa_ang2_p_list, sa_ang3_p_list]
        self.sa_ang = [sa_ang1_p_list, sa_ang2_p_list, sa_ang3_p_list]
        if minus_array:
            self.sa_ang_n = [sa_ang1_n_list, sa_ang2_n_list, sa_ang3_n_list]
            self.saa_sa = [saa_sa_p_list, saa_sa_n_list]
        else:
            self.sa_ang_p = [sa_ang1_p_list, sa_ang2_p_list, sa_ang3_p_list]
            self.saa_sa = saa_sa_p_list

        self.saa_sc = [saa_sc_x_list, saa_sc_y_list, saa_sc_z_list]

        self.hga_earth = hga_earth
        self.hga_angles = [hga_angles_el, hga_angles_az]

        self.mga_earth = mga_earth
        self.mga_angles = [mga_angles_el, mga_angles_az]

        self.roll_angles = [roll_angle_1, roll_angle_2, roll_angle_3]

        self.structures_flag = True
        self.previous_tw = self.time.window

        return


    def __Geometry(self, boresight=''):

        #if self.geometry_flag is True and \
        #                self.time.window.all() == self.previous_tw.all():
        #    return

        distance = []
        altitude = []
        boresight_latitude = []
        boresight_longitude = []
        latitude = []
        longitude = []
        subpoint_xyz = []
        subpoint_pgc = []
        subpoint_pcc = []
        zaxis_target_angle = []
        myaxis_target_angle = []
        yaxis_target_angle = []
        xaxis_target_angle = []
        beta_angle = []

        qs, qx, qy, qz = [], [], [] ,[]
        x, y, z = [],[],[]


        tar = self.target
        time = self.time

        for et in time.window:

            try:
                #
                # Compute the distance
                #
                ptarg, lt = spiceypy.spkpos(tar.name, et, tar.frame, time.abcorr,
                                          self.name)
                x.append(ptarg[0])
                y.append(ptarg[1])
                z.append(ptarg[2])

                vout, vmag = spiceypy.unorm(ptarg)
                distance.append(vmag)


                #
                # Compute the geometric sub-observer point.
                #
                if tar.frame == 'MARSIAU':
                    tar_frame = 'IAU_MARS'
                else:
                    tar_frame = tar.frame
                spoint, trgepc, srfvec = spiceypy.subpnt(tar.method, tar.name, et,
                                                       tar_frame, time.abcorr,
                                                       self.name)
                subpoint_xyz.append(spoint)

                #
                # Compute the observer's altitude from SPOINT.
                #
                dist = spiceypy.vnorm(srfvec)
                altitude.append(dist)


                #
                # Convert the sub-observer point's rectangular coordinates to
                # planetographic longitude, latitude and altitude.
                #
                spglon, spglat, spgalt = spiceypy.recpgr(tar.name, spoint,
                                                       tar.radii_equ, tar.flat)

                #
                # Convert radians to degrees.
                #
                spglon *= spiceypy.dpr()
                spglat *= spiceypy.dpr()

                subpoint_pgc.append([spglon, spglat, spgalt])

                #
                #  Convert sub-observer point's rectangular coordinates to
                #  planetocentric radius, longitude, and latitude.
                #
                spcrad, spclon, spclat = spiceypy.reclat(spoint)


                #
                # Convert radians to degrees.
                #
                spclon *= spiceypy.dpr()
                spclat *= spiceypy.dpr()

                subpoint_pcc.append([spclon, spclat, spcrad])
                latitude.append(spclat) #TODO: Remove with list extraction
                longitude.append(spclon)  # TODO: Remove with list extraction

                #
                # Compute the geometric sub-boresight point.
                #
                if tar.frame == 'MARSIAU':
                    tar_frame = 'IAU_MARS'
                else:
                    tar_frame = tar.frame


                if boresight:
                    try:
                        id = spiceypy.bodn2c(boresight)
                        (shape,framen, bsight, n, bounds) = spiceypy.getfov(id, 80)
                        mat = spiceypy.pxform(framen,tar_frame,et)
                    except:
                        framen = boresight
                        bsight = 0,0,1
                else:
                    bsight = self.name

                try:
                    if tar.method == 'INTERCEPT/ELLIPSOID':
                        method = 'ELLIPSOID'
                    else:
                        method = tar.method
                    spoint, trgepc, srfvec = spiceypy.sincpt(method, tar.name, et,
                                                           tar_frame, time.abcorr,
                                                           self.name, framen,
                                                           bsight)

                    #
                    # Convert the sub-observer point's rectangular coordinates to
                    # planetographic longitude, latitude and altitude.
                    #
                    spglon, spglat, spgalt = spiceypy.recpgr(tar.name, spoint,
                                                           tar.radii_equ, tar.flat)

                    #
                    # Convert radians to degrees.
                    #
                    spglon *= spiceypy.dpr()
                    spglat *= spiceypy.dpr()


                    #
                    #  Convert sub-observer point's rectangular coordinates to
                    #  planetocentric radius, longitude, and latitude.
                    #
                    spcrad, spclon, spclat = spiceypy.reclat(spoint)


                    #
                    # Convert radians to degrees.
                    #
                    spclon *= spiceypy.dpr()
                    spclat *= spiceypy.dpr()

                    boresight_latitude.append(spclat)
                    boresight_longitude.append(spclon)

                except:
                    pass

                #
                # Compute the angle between the observer's S/C axis and the
                # geometric sub-observer point
                #
                obs_tar, ltime = spiceypy.spkpos(tar.name, et,
                                                       'J2000', time.abcorr,
                                                       self.name)
                obs_zaxis  = [0,  0, 1]
                obs_myaxis = [0, -1, 0]
                obs_yaxis = [0, 1, 0]
                obs_xaxis = [1, 0, 0]

                #
                # We need to account for when there is no CK attitude available.
                #
                try:
                    matrix = spiceypy.pxform(self.frame, 'J2000', et)

                    z_vecout = spiceypy.mxv(matrix, obs_zaxis)
                    zax_target_angle = spiceypy.vsep(z_vecout, obs_tar)
                    zax_target_angle *= spiceypy.dpr()
                    zaxis_target_angle.append(zax_target_angle)

                    my_vecout = spiceypy.mxv(matrix, obs_myaxis)
                    myax_target_angle = spiceypy.vsep(my_vecout, obs_tar)
                    myax_target_angle *= spiceypy.dpr()
                    myaxis_target_angle.append(myax_target_angle)

                    y_vecout = spiceypy.mxv(matrix, obs_myaxis)
                    yax_target_angle = spiceypy.vsep(y_vecout, obs_tar)
                    yax_target_angle *= spiceypy.dpr()
                    yaxis_target_angle.append(yax_target_angle)

                    x_vecout = spiceypy.mxv(matrix, obs_myaxis)
                    xax_target_angle = spiceypy.vsep(x_vecout, obs_tar)
                    xax_target_angle *= spiceypy.dpr()
                    xaxis_target_angle.append(xax_target_angle)


                    quat = spiceypy.m2q(spiceypy.invert(matrix))
                    qs.append(quat[0])
                    qx.append(-1*quat[1])
                    qy.append(-1*quat[2])
                    qz.append(-1*quat[3])

                except:
                    zaxis_target_angle.append(0.0)
                    myaxis_target_angle.append(0.0)
                    yaxis_target_angle.append(0.0)
                    xaxis_target_angle.append(0.0)
                    qs.append(0.0)
                    qx.append(0.0)
                    qy.append(0.0)
                    qz.append(0.0)

                beta_angle.append(spiops.beta_angle(self.name, self.target.name,
                                                    et))
            except:
                boresight_latitude = 0
                boresight_longitude = 0
                distance = 0
                altitude = 0
                latitude = 0
                longitude = 0
                subpoint_xyz = [0,0,0]
                subpoint_pgc =  [0,0,0]
                subpoint_pcc =  [0,0,0]
                zaxis_target_angle = 0
                myaxis_target_angle = 0
                yaxis_target_angle = 0
                xaxis_target_angle = 0
                beta_angle = 0
                (qx, qy, qz, qs) = 0, 0, 0, 0
                (x, y, z) = 0, 0, 0

        self.boresight_latitude = boresight_latitude
        self.boresight_longitude = boresight_longitude
        self.distance = distance
        self.altitude = altitude
        self.latitude = latitude
        self.longitude = longitude
        self.subpoint_xyz = subpoint_xyz
        self.subpoint_pgc = subpoint_pgc
        self.subpoint_pcc = subpoint_pcc
        self.zaxis_target_angle = zaxis_target_angle
        self.myaxis_target_angle = myaxis_target_angle
        self.yaxis_target_angle = yaxis_target_angle
        self.xaxis_target_angle = xaxis_target_angle
        self.beta_angle = beta_angle
        self.quaternions = [qx, qy, qz, qs]
        self.trajectory = [x,y,z]

        self.geometry_flag = True
        self.previous_tw = self.time.window

        return


    def Plot(self, yaxis = 'distance', date_format='TDB', external_data=[],
             notebook=False, boresight=''):

        self.__Geometry(boresight=boresight)
        self.__Structures()

        #
        # Not Time in X axis
        #
        if yaxis == 'groundtrack':
            utils.plot(self.__getattribute__('longitude'),
                       self.__getattribute__('latitude'),
                       notebook=notebook, xaxis_name = 'Longitude',
                       yaxis_name='Latitude', mission=self.name,
                       target=self.target.name, background_image=True,
                       format ='circle_only')

            return

        if yaxis == 'boresight_groundtrack':
            utils.plot(self.__getattribute__('boresight_longitude', boresight=boresight),
                       self.__getattribute__('boresight_latitude', boresight=boresight),
                       notebook=notebook, xaxis_name = 'Longitude',
                       yaxis_name='Latitude', mission=self.name,
                       target=self.target.name, background_image=True,
                       format ='circle_only')

            return

        elif yaxis == 'trajectory':

            if notebook:
                output_notebook()

            # Make data
            x = float(self.target.radii[0])
            y = float(self.target.radii[1])
            z = float(self.target.radii[2])


            # create a new plot
            s1 = figure(width=300, plot_height=300, title='X-Y [KM]')
            s1.ellipse(x=0, y=0, width=x, height=y, color="orange")

            s1.line(self.__getattribute__('trajectory')[0],
                      self.__getattribute__('trajectory')[1],
                      color="red")

            # create another one
            s2 = figure(width=300, height=300, title='X-Z [KM]')
            s2.ellipse(x=0, y=0, width=x, height=z, color="orange")
            s2.line(self.__getattribute__('trajectory')[0],
                        self.__getattribute__('trajectory')[2],
                        color="green")


            # create another one
            s3 = figure(width=300, height=300, title='Y-Z [KM]')
            s3.ellipse(x=0, y=0, width=y, height=z, color="orange")
            s3.line(self.__getattribute__('trajectory')[1],
                        self.__getattribute__('trajectory')[2],
                        color="blue")


            # put all the plots in an HBox
            show(row(s1, s2, s3))

            return

        #
        # Time in X axis
        #

        xaxis = self.time.window
        xaxis_name = 'Date'

        if (yaxis == 'sa_ang') or (yaxis == 'sa_ang_p'):
            yaxis_name = ['sa_ang1_p','sa_ang2_p', 'sa_ang3_p']
            yaxis_units = 'deg'
        elif yaxis == 'sa_ang_n':
            if self.name != 'MPO':
                yaxis_name = ['sa_ang1_n', 'sa_ang2_n', 'sa_ang3_n']
                yaxis_units = 'deg'
            else:
                yaxis_name = ['sa_ang1_p','sa_ang2_p', 'sa_ang3_p']
                yaxis_units = 'deg'
        elif yaxis == 'saa_sc':
            yaxis_name = ['saa_sc_x', 'saa_sc_y', 'saa_sc_z']
            yaxis_units = 'deg'
        elif yaxis == 'saa_sa':
            if self.name != 'MPO':
                yaxis_name = ['saa_sa_p', 'saa_sa_n']
                yaxis_units = 'deg'
            else:
                yaxis_name = 'saa_sa'
                yaxis_units = 'deg'
        elif yaxis == 'hga_angles':
            yaxis_name = ['hga_el', 'hga_az']
            yaxis_units = 'deg'
        elif yaxis == 'mga_angles':
            yaxis_name = ['mga_el', 'mga_az']
            yaxis_units = 'deg'
        elif yaxis == 'roll_angles':
            yaxis_name = ['roll_1', 'roll_2', 'roll_3']
            yaxis_units = 'deg'
        elif yaxis == 'hga_earth':
            yaxis_name = 'hga_earth'
            yaxis_units = 'deg'
        elif yaxis == 'mga_earth':
            yaxis_name = 'mga_earth'
            yaxis_units = 'deg'
        elif yaxis == 'quaternions':
            yaxis_name = ['qx','qy','qz','qs']
            yaxis_units = ''
        elif yaxis == 'clock_drift':
            self.__ClockDrift()
            xaxis = self.clock_dates
            xaxis_name = 'Date'
            yaxis_name = 'Delta Clock Counts SC-Ground'
            yaxis_units = 's'
        elif yaxis == 'zaxis_target_angle':
            yaxis_units = 'deg'
            yaxis_name = yaxis
        elif yaxis == 'yaxis_target_angle':
            yaxis_units = 'deg'
            yaxis_name = yaxis
        elif yaxis == 'xaxis_target_angle':
            yaxis_units = 'deg'
            yaxis_name = yaxis
        else:
            yaxis_name = yaxis
            yaxis_units = 'km'


        utils.plot(xaxis, self.__getattribute__(yaxis), notebook=notebook,
                   external_data=external_data, xaxis_name=xaxis_name,
                   yaxis_name=yaxis_name, mission = self.name, yaxis_units=yaxis_units,
                   target = self.target.name, date_format=date_format)

        return



    def Plot3D(self, data='trajectory', reference_frame=False):

        self.__Geometry()

        #TODO: Arrange the reference frame flow
        if not self.state_in_window:
            self.__StateInWindow(reference_frame=reference_frame)


        data = self.state_in_window

        utils.plot3d(data, self, self.target)

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

        try:
            self.radii = spiceypy.bodvar(self.id, 'RADII', 3)

        except:
            print("Ephemeris object has no radii")
            return

        self.radii_equ = self.radii[0]
        self.radii_pol = self.radii[2]
        self.flat = (self.radii_equ - self.radii_pol) / self.radii_equ

        return


class Observer(Body):
    def __init__(self, body, time=object(), target=False, frame='', mission_config=None):

        super(Observer, self).__init__(body, time=time, target=target, mission_config=mission_config)

        if not frame:
            self.frame = utils.get_frame(self.name)
            if spiceypy.namfrm(self.frame) == 0:
                self.frame = self.name

            if spiceypy.namfrm(self.frame) == 0:
                #TODO: Fix this shit
                self.frame = '{}_LANDER'.format(self.name)
        else:
            self.frame = frame

        if spiceypy.namfrm(self.frame) == 0:
            print('The frame name was not recognized. Frame: ' + str(self.frame))
