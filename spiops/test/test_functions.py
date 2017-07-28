import pytest
import spiops

def test_fov_illum():

    angle = spiops.fov_illum(mk='mex_ops.tm',
                             sensor='MEX_VMC',
                             time='2016-10-27')

    assert angle == 126.98764786304272
