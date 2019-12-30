from spiops import spiops

def test_brief():

    output = spiops.naif.brief(mk='MEX_OPS.TM',
                            object='MEX_SC_REF',
                            time_format='UTC')

    print(output)

