## What is it?

**pynkpandas** is a plugin for **pynk** - a Python wrapper for **MUSTLink**.
**pynkpandas** combines the generic **pynk** data representations with [**pandas**]
(https://github.com/pandas-dev/pandas).

It does support spacecraft mission engineers and data scientists to combine spacecraft data with 
the flexibility of pandas.


## Main Features
**pynkpandas** can support you on:

  - Getting numerical and textual telemetry timeseries data from **ARES** or **MUST** as **pandas** series and dataframes.
  - Getting tabular data like Telecommand History, event data and alphanumeric displays as **pandas** dataframes.
  - Getting the computation result of **MUSE** scripts as Pandas Dataframes.
  - Search and retrieve parameter **metadata**
  
  
## Where to get it
The source code is currently hosted on ESA Gitlab at:
https://gitlab.esa.int/mustlink-dev/pynkpandas

The module can be downloaded and installed from an internal Pypi server available under:
http://damachine.esoc.esa.int:8999/simple/pynkpandas/




```sh
# PyPI
pip install --trusted-host damachine.esoc.esa.int --extra-index-url http://damachine.esoc.esa.int:8999/ pynkpandas
```

## Dependencies

- [pynk](https://gitlab.esa.int/mustlink-dev/pynkpandas):0.0.1 or above
- [pandas](https://github.com/pandas-dev/pandas):0.23.4 or above


## License
TBD
