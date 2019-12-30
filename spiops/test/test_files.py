from spiops import spiops

def test_brief():

    updated_mk = spiops.update_former_versions('/Users/mcosta/SPICE/ROSETTA/kernels/mk',
                                                 '/Users/mcosta/SPICE/ROSETTA/kernels')

    print(updated_mk)

