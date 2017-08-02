import numpy as np

class TimeWindow(object):
    """
    This class defines a Time Window and is the lowest level time entity for
    spiops. The time window is generated with a start and finish times and it
    can also indicate a current instant.
    """

    def __init__(self, start, finish, current=False, resolution=False,
                 format='UTC'):
        """
        This method creates a new Time Window.

        :param start: Start time of the time window
        :type start: Union[str, float]
        :param finish: Finish time of the time window
        :type finish: Union[str, float]
        :param current: Selected time within the time window (used for computations that might require a current time)
        :type current: Union[str, float]
        :param resolution: Resolution of the time window in seconds
        :type resolution: float
        :param format: Calendar format in which we want to see times ('UTC' or 'CAL')
        :type format: str
        """
        import spiops.utils.time as zztime


        if isinstance(start, str):
            start = zztime.cal2et(start, format=format)
        if isinstance(finish, str):
            finish = zztime.cal2et(finish, format=format)

        self.start = start
        self.finish = finish
        self.now = current
        self.res = resolution
        self.format = format
        self.time_set = []


    def __buildWindow(self):
        """
        Internal method to generate the time set of the time window.
        """

        #
        # if no resolution is given a default ten element sample is computed
        #
        if not self.res or self.res >= self.finish - self.start:
            print('No resolution provided or resolution is bigger than the '
                  'time interval. Setting to default')
            self.res = (self.finish - self.start)/10.0

        self.time_set = np.arange(self.start, self.finish, self.res)


    def __getattribute__(self, item):
        """
        Internal method to overwrite get attribute.
        """

        if item in ['time_set']:
            self.__buildWindow()
            return object.__getattribute__(self, item)
        else:
            return object.__getattribute__(self, item)


    def getTime(self, item):
        """
        Method to obtain the Start, Finish or Current time of the Time Window.

        :param item: Time attribute that we want to get from the class: 'start', 'finish' or 'current'
        :type item: str
        :return: Requested time attribute with the specified format
        :rtype: str
        """
        import spiops.utils.time as zztime

        if item in ['start', 'finish', 'current']:
            return zztime.et2cal(object.__getattribute__(self, item),
                                 self.format)


    def getTimeInterval(self, format=False):
        """
        Method to obtain the Time Interval in a given time format.

        :param format: Input format; 'UTC' or 'CAL'
        :type format: str
        :return: Time Interval in the provided time format
        :rtype: list
        """

        import spiops.utils.time as zztime

        if not format:
            format = self.format

        return zztime.et2cal(self.time_set, format)

