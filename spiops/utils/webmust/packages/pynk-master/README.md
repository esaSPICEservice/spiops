## What is it?

**pynk** is a Python wrapper for **MUSTLink**. A REST-API that provides fast and easy 
access to spacecraft data parameter archives such as **ARES** or **MUST** in a unified way.
 **MUSTLink** also provides a way to interface the service layer on top the parameter
 data archives. It is possible to interface **DrMUST**, **Novelty Detection**, **muse**,
 **user management** and other features like **WebMUST** profiles and data 
 exportation using **Mustlink**.



## Main Features
**pynk** can support you on:

  - Getting telemetry data from **ARES** or **MUST**, including numerical and textual data,
      telecommand history and events.
  - Create, manage and execute Javascript or OL code to create synthetic data.
  - Create, read, update and delete user related data.
  - Create, read, update and delete profiles related data.
  - Interface **DrMUST** and **Novelty Detection** provided services.
  - Export existing data into multiple formats like DARC, CSV and ASCII.
  - Search and navigate in the existing telemetry. 



## Where to get it
The source code is currently hosted on ESA Gitlab at:
https://gitlab.esa.int/mustlink-dev/pynk

The module can be downloaded and installed from an internal Pypi server available under:
http://damachine.esoc.esa.int:8999/simple/pynk/

```sh
# PyPI
pip install --trusted-host damachine.esoc.esa.int --extra-index-url http://damachine.esoc.esa.int:8999/ ops_pynk
```



## Dependencies

- [Requests](https://github.com/requests/requests):2.20.0 or above
- [Pandas](https://github.com/pandas-dev/pandas):0.23.0 or above




## How to use it
**pynk** supply a core class `Mustlink`, in order to manage the communication from and to the MUSTLINK API REST provider, 
and some specialised classes like `Users`, `Groups` and `DataProviders` to map actual API endpoints to python functions. 

All specialised classes accept and return standard python objects; return objects represent JSON response from the API. 
Specialised classes supply only slightly higher logic and perform only basic coherency checks.

Around **pynk** core some high level functions have been developed and supplied by files organized per area of responsibility: 

- `user.py` for user management
- `data.py` for data extraction

**pynk** module supports multiple MUSTLINK instances but has simplified usage when only one default instance is provided: 
to activate the default settings see `pynk.init_defaults` function, after activated is possible to invoke functions from 
specialised classes without specifying the `Mustlink` instance.



**Configure and use pynk with default instance:**

```python
import datetime
import pynk
from pynk import userpynk, datapynk

# configuration dictionary must have the following keys
conf = {
    'base_url': "http://must.site.esa.int:8080/mustlink-api",
    'username': "test_user",
    'password': "test_password"
}

# initialize internal pynk status with the default Mustlink instance
pynk.init_defaults(conf)

# cycle through all users, using the high-level logic user.py script
print("List of all users' login")
for usr in userpynk.list_users():
    print(usr.login)

# load data from data provider SENTINEL2A
some_data_columns = datapynk.load_data('SENTINEL2A', 
                                       ['AJTD0360', 'AJTD0380', 'AJTD0390', 'AJTD0860'], 
                                       "2018-10-01 00:00:00", datetime.datetime(2018, 10, 2, 0, 0, 0), 
                                       aggregation_function="AVG", aggregation_type="Minutes", aggregation_value=30)

for data_col in some_data_columns:
    print(data_col['metadata'])
    # ... 
    # ...

```



**Configure and use pynk with multiple instance:** 

```python
import pynk
from pynk import userpynk, datapynk

conf1 = {
    'base_url': "http://must.site.esa.int:8080/mustlink-api",
    'username': "test_user",
    'password': "test_password"
}
conf2 = {
    'base_url': "http://must.othersite.esa.int:8080/mustlink-api",
    'username': "test_user",
    'password': "test_password"
}

# get 2 different instances of Mustlink, respectively configured with different MUSTLINK site
ml1 = pynk.get(conf1)
ml2 = pynk.get(conf2)

print("List of all users' login from site 1")
for usr in userpynk.list_users(mustlink=ml1):  # mustlink instance must be specified
    print(usr.login)

print("Load data from site 2")
datapynk.load_data('SENTINEL2A', 
                   ['AJTD0360', 'AJTD0380', 'AJTD0390', 'AJTD0860'], 
                   "2018-10-01 00:00:00", "2018-10-02 00:00:00",  
                   aggregation_function="AVG", aggregation_type="Minutes", aggregation_value=30,
                   mustlink=ml2) # mustlink instance must be specified
```



**Access to specialised classes from pynk:**

```python
import pynk

mustlink = pynk.get({
        'base_url': "http://must.site.esa.int:8080/mustlink-api",
        'username': "test_user",
        'password': "test_password"
    })

print("User 'test_user':")
# access to the specialised bridge object 'Users' to request a user JSON representation
print(mustlink.users().load('test_user'))


s2a_provider = mustlink.data_providers().load_data_provider('SENTINEL2A')
s2a_provider.load_parameters_data(['AJTD0360', 'AJTD0380', 'AJTD0390', 'AJTD0860'], 
                                  "2018-10-01 00:00:00", "2018-10-02 00:00:00")


```


**Obtain a pandas.DataFrame from the pynk_panda extension script:**

```python
import pynk
from pynk.ext import pynk_panda

pynk.init_defaults({
    'base_url': "http://must.site.esa.int:8080/mustlink-api",
    'username': "test_user",
    'password': "test_password"
})

print("Print first 100 values ")
param_df = pynk_panda.get_aggregated_parameter_data("SENTINEL2A", ["AJTD0360", "AJTD0380", "AJTD0390", "AJTD0860"],
                                                    "2018-10-01 00:00:00", "2018-10-02 00:00:00", calibrated=False)
print(param_df.head(100))
```


**A complete example using PYNK and PYNK_PANDA**

```python

import pynk
from pynk.ext import pynk_panda

config = {
    'base_url': "http://your.mustlink.installation/mustlink",
    'username': "<your username>",
    'password': "<your password>"
}

#Instantiate a MUSTLink object. 
mlink = pynk.get(config)

#Search for all the dataprovider objects available on a Mustlink instance
print(mlink.data_providers().search('Int'))

#Return a dataprovider object
cryosat2 = mlink.data_providers().load_data_provider('Cryosat2')
#xmm = mlink.data_providers().load_data_provider('XMM')
#mex = mlink.data_providers().load_data_provider('MEX')


#Get all the parameters available on a given dataprovider
cryosat2.get_parameters()

#Get all the parameters available on a given dataprovider, matching the search_text on a dataframe
pynk_panda.get_parameters('Cryosat2', 'AST5011', mustlink=mlink)


#Get all the details for the given parameter
cryosat2.get_parameter_details('AST50116')


#Get all the aggregations available for a given dataprovider
cryosat2.get_aggregation_periods()
cryosat2.get_aggregation_functions()

#Get all the available statistics for the given parameter list
cryosat2.get_statistics(['AST50115','AST50116'],'2018-11-12 00:00:00', '2018-11-12 08:00:00')

#Get all the available statistics for the given parameter list as a dataframe
pynk_panda.get_statistics('Cryosat2',['AST50115','AST50116'],'2018-11-12 00:00:00', '2018-11-12 08:00:00', mustlink=mlink)

# The following approach lead to a one timer request that only retrieves when all the data is available. Very long requests may lead to timeout or server failiure 
cryosat2.get_parameters_data(['AST50116'],'2018-11-12 00:00:00', '2018-11-12 02:00:00')

cryosat2.get_parameters_data(['AST50116', 'AST50115'],'2018-11-12 00:00:00', '2018-11-12 02:00:00',
                                   aggregation_function="AVG", aggregation_type="Hours", aggregation_value=1)

# The following approach retrieves the data from the API following a stream fashioned way that splits the request in smaller chunks making it safer for longer requests
datapynk.load_data('Cryosat2', ['AST50116'],'2018-11-08 00:00:00', '2018-11-12 00:00:10', mustlink=mlink)


# The following approach provide you the requested data as a pandas dataframe. The request is safe since it is using the datapynk specialized script
pynk_panda.get_parameter_data('Cryosat2',['AST50115','AST50116'],'2018-11-12 00:00:00', '2018-11-12 08:00:00', mustlink=mlink)


#Check what tables are available. The dataType field is used to request tables related data
cryosat2.get_tables()

#Get the metadata of a single table. It does include the headers and the table title
cryosat2.get_table_metadata('TC')

#Get Tables Data. The table mode can be BRIEF or FULL. BRIEF is taken by default.
cryosat2.get_table_data('TC','2018-11-04 00:00:00', '2018-11-04 00:30:00', mode='FULL', max_rows=5000)

#Get Tables data as a Pandas Dataframe
pynk_panda.get_table_data('Cryosat2','TC','2018-11-12 00:00:00', '2018-11-12 08:00:00', mode='FULL', max_rows=2000, mustlink=mlink)
```


## License
TBD
