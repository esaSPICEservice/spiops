class Data(object):
    """
    This object is a container of data. Data can be of any kind but not a
    SPICE kernel, although it can be the output of a SPICE kernel.
    """
    #TODO: Should adcsng use data objects as well?

    def __init__(self, start, finish, current=False, resolution=False,
                 abcorr='NONE', format='UTC'):

        return