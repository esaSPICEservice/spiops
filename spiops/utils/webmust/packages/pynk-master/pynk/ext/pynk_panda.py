import pandas as pd
import pynk
from pynk import datapynk
import json


def get_parameters(data_provider,  key=None, value=None, parameter_type=None, mustlink=None):
    dp = pynk.get(mustlink).data_providers.load_data_provider(data_provider)
    params = dp.get_parameters(key=key, value=value, parameter_type=parameter_type)
    dataframe = pd.DataFrame(params)
    dataframe.set_index('Id', inplace=True)
    return dataframe

def get_statistics(data_provider, parameter_names, start, end, mustlink=None, **kwargs):
    dp = pynk.get(mustlink).data_providers.load_data_provider(data_provider)
    stats = dp.get_statistics(parameter_names, start,end, **kwargs)
    dataframe = pd.DataFrame(stats)
    dataframe.set_index('parameter', inplace=True)
    return dataframe

def get_parameter_data(data_provider, parameter_names, start, end,
                                  calibrated=False,
                                  mustlink=None,
                                  **kwargs):

    param_values = datapynk.load_data(data_provider, parameter_names, start, end, mustlink=mustlink, **kwargs)
    return __dataframe(param_values, calibrated)

def get_table_data(data_provider, table, start, end, mode='BRIEF', max_rows=5000, filter_keys=[' '], filter_values=[' '], mustlink=None):
    dp = pynk.get(mustlink).data_providers.load_data_provider(data_provider)
    table = dp.get_table_data(table, start, end, mode=mode, max_rows=max_rows, filter_keys=filter_keys, filter_values=filter_values)
    headers_dict = table['headers']
    data_dict = table['data']
    dataframe = pd.DataFrame(data_dict)
    dataframe = dataframe.rename(columns=headers_dict)
    return dataframe




def __dataframe(param_values, calibrated):
    df = pd.DataFrame([])

    for param_value in param_values:
        dates = []
        values = []
        for element in param_value[1]:
            dates.append(element.date)
            if calibrated:
                values.append(element.calibrated_value)
            else:
                values.append(element.value)

        series = pd.Series(values, index=dates)
        series.index = pd.to_datetime(series.index, unit='ms')
        series.name = param_value[0].name

        df = pd.concat((df, series), axis=1)

    return df
