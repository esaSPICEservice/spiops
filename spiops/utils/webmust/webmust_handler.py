"""
Created June, 2021

@author: Ricardo Valles Blanco (ESAC)

This module allow to get TM info pynk API from MUST database
"""

import pandas as pd
import pynk
from pynk.ext import pynk_panda


class WebmustHandler(object):

    def __init__(self, config=None, mission_phase='BEPICRUISE'):

        if config:

            self.conf = config

        else:
            self.conf = {
                "username": "XXXXXXXX",
                "password": "XXXXXXXX",
                "base_url": "https://bepicolombo.esac.esa.int/webclient-must/mustlink"
            }

        self.data_provider = mission_phase

        pd.options.display.max_rows = 10000
        pd.options.display.max_columns = 5000

        # Initialize internal pynk status with the default Mustlink instance
        pynk.init_defaults(self.conf)

        # Instantiate a MUSTLink object.
        self.mlink = pynk.get(self.conf)

        self.bepic = self.mlink.data_providers.load_data_provider(mission_phase)

    def get_tm(self, param_list, start_date, stop_date, calibrated=False):
        df = pynk_panda.get_parameter_data(self.data_provider, param_list, start_date, stop_date, calibrated=calibrated)
        return df


if __name__ == '__main__':

    print('start test')

    start, stop = ('2020-04-09 03:14:58', '2020-04-09 04:14:58')
    param_list = ['NCADAF41', 'NCADAF42', 'MPO_SA_nudging_ESS_clonned_from_thomas']

    tm = WebmustHandler(mission_phase='BEPICRUISE')
    df_0 = tm.get_tm(param_list, start, stop)
    print(df_0)

