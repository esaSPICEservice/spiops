from setuptools import setup

setup(
    name='pynkpandas',
    packages=['pynkpandas'],
    description='pynkpandas is a plugin for pynk. It does support data analysis.',
    version='0.0.1',
    url='https://gitlab.esa.int/mustlink-dev/pynkpandas',
    author='esa/esoc/ops-os',
    install_requires=["pandas"],
    python_requires='>=3',
    dependency_links=["http://damachine.esoc.esa.int:8999/simple/pynk"],
    keywords=['pip','mustlink','ares', 'must', 'webmust', 'data analytics', 'muse', 'drmust', 'novelty detection', 'pyares', 'housekeeping', 'telemetry', 'timeseries', 'pynk', 'pandas']
    )