from setuptools import setup, find_packages

setup(
    name='ops_pynk',
    packages=find_packages(include=['pynk', 'pynk.*']),
    description='pynk is a python wrapper for MUSTLink.',
    version='0.4',
    url='https://gitlab.esa.int/mustlink-dev/pynk',
    author='esa/esoc/ops-os',
    install_requires=["requests", 'pandas'],
    python_requires='>=3',
    keywords=['pip', 'mustlink', 'ares', 'must', 'webmust', 'data analytics', 'muse', 'drmust', 'novelty detection',
              'pyares', 'housekeeping', 'telemetry', 'timeseries']
)
