from spiops import spiops


def test_ckdiff():

    resolution = 1000000
    tolerance = 0.0001 # deg
    spacecraft_frame = 'ROS_LANDER'
    target_frame = 'J2000'
    mk = 'data/SPICE/ROSETTA/mk/ROS_CKDIFF_TEST.TM'
    ck1 = 'data/SPICE/ROSETTA/mk/ROS_CKDIFF_TEST.TM'
    ck2 = 'data/SPICE/ROSETTA/ck/LATT_EME2LDR_SDL_SONC_V1_1.BC'

    spiops.ckdiff(mk, ck1, ck2, spacecraft_frame, target_frame,
                  resolution, tolerance)


def test_convert_ESOCorbit2data():

    print(spiops.utils.convert_ESOCorbit2data('data/ESOC/fdy/LORB_EME2000_RBD_1_V1_0.ROS'))


def test_valid_url():

    out = spiops.utils.valid_url('67P/CG')

    assert out == '67P-CG'


def test_cal2et():

    time1 = spiops.time.cal2et('2000-01-01T12:00:00', format='CAL',
                               support_ker='MEX_OPS.TM', unload=True)

    time2 = spiops.time.cal2et('2000-01-01T12:00:00', format='UTC',
                               support_ker='MEX_OPS.TM', unload=True)

    assert time1 == 0.0
    assert time2 == 64.18392728473108


def test_et2cal():

    time1 = spiops.time.et2cal(0.0, format='UTC', support_ker='MEX_OPS.TM',
                               unload=True)
    time2 = spiops.time.et2cal(64.18392728473108, format='UTC',
                               support_ker='MEX_OPS.TM',unload=True)

    assert time1 == '2000-01-01T11:58:55.816'
    assert time2 == '2000-01-01T12:00:00.000'


def test_mjd20002et():


    time = spiops.time.mjd20002et(6863.0790, support_ker='naif0012.tls',
                                  unload=True)

    assert time == 592926894.7823608


def test_fov_illum():

    angle = spiops.fov_illum(mk='MEX_OPS.TM',
                             sensor='MEX_VMC',
                             time='2017-04-01', unload=True)

    assert angle == 51.4628080108263


def test_cov_spk_obj():

    cov = spiops.cov_spk_obj(mk='MEX_OPS.TM',
                             object='MEX',
                             time_format='UTC',
                             unload=True)

    assert cov == [['2017-05-31T23:35:23.939', '2017-06-30T20:46:40.995'],
                   ['2017-04-30T23:42:46.000', '2017-05-31T23:35:23.939'],
                   ['2017-03-31T23:50:43.000', '2017-04-30T23:42:46.000'],
                   ['2017-02-28T22:01:45.000', '2017-03-31T23:50:43.000'],
                   ['2017-01-31T22:42:20.000', '2017-02-28T22:01:45.000'],
                   ['2016-12-31T21:29:46.000', '2017-01-31T22:42:20.000']]


def test_cov_spk_ker():

    cov = spiops.cov_spk_ker(spk='/Users/mcosta/Dropbox/SPICE/SPICE_MEX/ftp/data/SPICE/MARS-EXPRESS/kernels'
                            '/spk/ORMM_T19_170601000000_01351.BSP',
                             support_ker='/Users/mcosta/Dropbox/SPICE/SPICE_MEX/ftp/data/SPICE/MARS-EXPRESS/kernels/lsk'
                             '/NAIF0012.TLS',
                             object='MEX',
                             time_format='UTC')

    assert cov == ['2017-05-31T23:35:23.939', '2017-06-30T20:46:40.995']


def test_cov_ck_obj():

    cov = spiops.cov_ck_obj(mk='MEX_OPS.TM',
                            object='MEX_SC_REF',
                            time_format='UTC')

    assert cov == [['2016-12-31T23:58:51.815', '2017-07-12T14:58:50.651'],
                   ['2017-03-31T23:50:43.000', '2017-06-27T06:08:19.973']]


def test_cov_ck_ker():

    cov = spiops.cov_ck_ker(ck='/Users/mcosta/Dropbox/SPICE/SPICE_MEX/ftp/data'
                               '/SPICE/MARS-EXPRESS/kernels/ck/ATNM_VMC_170101_170712_V01.BC',
                            support_ker='MEX_OPS.TM',
                            object='MEX_SC_REF',
                            time_format='UTC')

    assert cov ==  ['2016-12-31T23:58:51.815', '2017-07-12T14:58:50.651']


def test_fk_body_ifj2000():

    transf = spiops.fk_body_ifj2000('JUICE', 'JUPITER',
              '/Users/mcosta/Dropbox/SPICE/SPICE_JUICE/ftp/data/SPICE/JUICE'
              '/kernels/pck/pck00010.tpc',
                 '/Users/mcosta/Dropbox/SPICE/SPICE_JUICE/ftp/data/SPICE/JUICE'
              '/kernels/spk/jup310.bsp', ',-28970', report=False, file=False,
              unload=True)

    assert transf == 'JUPITER_IF->J2000 (3-2-3): -88.05720404270757 - ' \
                     '25.504190046604286 - -48.96872816579391'