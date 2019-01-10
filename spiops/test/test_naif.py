from spiops import spiops

def test_brief():

    output = spiops.naif.brief(mk='MEX_OPS.TM',
                            object='MEX_SC_REF',
                            time_format='UTC')

    assert cov == [['2016-12-31T23:58:51.815', '2017-07-12T14:58:50.651'],
                   ['2017-03-31T23:50:43.000', '2017-06-27T06:08:19.973']]