"""
Data management script
"""
import collections
from datetime import datetime, timedelta

import pynk


def load_data(data_provider,
              param_names, start, end,
              time_slice_seconds=None,
              mustlink=None,
              **kwargs):
    """
    Retrieve data from MUSTLINK and transform metadata/values lists into object with named tuples Column and Data
    :param data_provider:
    :param param_names:
    :param start:
    :param end:
    :param time_slice_seconds:
    :param mustlink:
    :param kwargs: see at pynk.DataProvider.load_parameters_data:
                 calibrate=None,
                 aggregation_function=None, aggregation_type=None, aggregation_value=None,
                 compression_error=None,
                 chunk_count=None
    :return: a list of tuples for each column, with Column tuple as definition and a list of Data tuples
    """

    Data = collections.namedtuple('Data', 'date, value, calibrated_value')

    dp = pynk.get(mustlink).data_providers.load_data_provider(data_provider)
    if not dp:
        raise ValueError("Data provider %s doesn't exist" % data_provider)

    def retrieve(_start, _end):
        def map_result(api_result):

            #metadata is variable. Different dataproviders supply different metadata parameters, so we can not predefine it
            metadata = list(api_result.keys())
            metadata.remove('data')
            # create the row as a list to ensure that we can populate the namedtuple
            row = [api_result[key] for key in metadata]
           
            Column = collections.namedtuple('Column', metadata)

            #build the namedtuple from the given list
            col = Column(*row)
  
            values = [Data(int(d['date']), d['value'], d['calibratedValue']) for d in api_result['data']]
            
            return col, values

        return [map_result(param_value) for param_value in dp.get_parameters_data(param_names, _start, _end, **kwargs)]

    if time_slice_seconds:
        real_start = start if type(start) == datetime else datetime.strptime(start, pynk.DT_FMT)
        real_end = end if type(end) == datetime else datetime.strptime(end, pynk.DT_FMT)

        return _TimeSlicerIterable2(retrieve, real_start, real_end, time_slice_seconds)
    else:
        return retrieve(start, end)


class _TimeSlicerIterable2(list):
    def __init__(self, retrieve_func, real_start, real_end, time_slice_seconds):
        super().__init__()

        def time_slice_generator():
            _start = real_start - timedelta(seconds=1)
            _end = _start

            while _end < real_end:
                _start = _end + timedelta(seconds=1)  # from and to time boundaries are inclusive
                _end = _start + timedelta(seconds=time_slice_seconds)
                if _end > real_end:
                    _end = real_end

                yield _start, _end

        self.retrieve_func = retrieve_func
        self.time_slices = time_slice_generator

    def __iter__(self):
        return self.__slice_generator()

    def __len__(self):
        return 1

    def __slice_generator(self):
        for _from, _to in self.time_slices():
            buffer = self.retrieve_func(_from, _to)
            for i in range(0, len(buffer)):
                yield buffer[i]
