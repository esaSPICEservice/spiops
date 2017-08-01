from spiops import spiops

def test_fov_illum():

    angle = spiops.fov_illum(mk='MEX_OPS.TM',
                             sensor='MEX_VMC',
                             time='2017-04-01')

    assert angle == 51.4628080108263


def test_cov_spk_obj():

    cov = spiops.cov_spk_obj(mk='MEX_OPS.TM',
                             object='MEX',
                             time_format='UTC')

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
                   ['2017-03-31T23:50:43.000', '2017-06-27T06:08:19.973'],
                   ['2016-12-31T23:58:51.815', '2017-07-12T14:58:50.651'],
                   ['2017-03-31T23:50:43.000', '2017-06-27T06:08:19.973']]


def test_cov_ck_ker():

    cov = spiops.cov_ck_ker(ck='/Users/mcosta/Dropbox/SPICE/SPICE_MEX/ftp/data'
                               '/SPICE/MARS-EXPRESS/kernels/ck/ATNM_VMC_170101_170712_V01.BC',
                            support_ker='MEX_OPS.TM',
                            object='MEX_SC_REF',
                            time_format='UTC')

    assert cov ==  ['2016-12-31T23:58:51.815', '2017-07-12T14:58:50.651']
