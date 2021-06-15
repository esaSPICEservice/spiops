import json
import logging
import re
from datetime import datetime

import requests

log = logging.getLogger(__name__)
RE_ALL = re.compile(".*")

DT_FMT = '%Y-%m-%d %H:%M:%S'

DATA_PROVIDERS_URL = '/dataproviders'

USER_MANAGEMENT_URL = '/usermanagement'
USER_URL = '/usermanagement/users'
GROUPS_URL = '/usermanagement/groups'
SCRIPTS_URL = '/scripts'
GROUP_URL = '/usermanagement/group'
PROJECTS_URL = '/usermanagement/projects'
ROLES_URL = '/usermanagement/roles'
DRMUST_URL = '/analysis/drmust'
NOVELTY_URL = '/analysis/novelty'

# default mustlink instance
__mustlink = None


class MustlinkException(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class Mustlink:
    """
    Wrapping class for MUSTLINK server.
    Manages authentication and supply methods for composing parameters.
    """

    def __init__(self, base_url, username, password):
        self.base_url = base_url
        self.username = username
        self.password = password

        self.headers = {
            'Content-type': 'application/json',
        }

        self.users = Users(self)
        self.groups = Groups(self)
        self.projects = Projects(self)
        self.roles = Roles(self)
        self.data_providers = DataProviders(self)
        self.exports = Exports(self)
        self.drmust = DrMust(self)
        self.novelty = Novelty(self)
        self.scripts = Scripts(self)

    def __enter__(self):
        self.__authenticate()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def get(self, *endpoint, parameters=None):
        return self.__do_request(requests.get,
                                 url=self.__prepare_url(endpoint),
                                 params=self.__prepare_params(parameters),
                                 headers=self.headers)

    def put(self, *endpoint, parameters=None, body=None):
        return self.__do_request(requests.put,
                                 url=self.__prepare_url(endpoint),
                                 params=self.__prepare_params(parameters),
                                 data=json.dumps(body) or "",
                                 headers=self.headers)

    def post(self, *endpoint, parameters=None, body=None):
        return self.__do_request(requests.post,
                                 url=self.__prepare_url(endpoint),
                                 params=self.__prepare_params(parameters),
                                 data=json.dumps(body) or "",
                                 headers=self.headers)

    def delete(self, *endpoint):
        return self.__do_request(requests.delete,
                                 url=self.__prepare_url(endpoint),
                                 headers=self.headers)

    def __do_request(self, requests_method, **param):
        # double loop if authentication expired
        for force_flag in [False, True]:
            self.__authenticate(force=force_flag)

            response = requests_method(**param)

            if 200 <= response.status_code < 300:
                return response.json() if response.content else None
            elif response.status_code == 401:
                # tries one again to authenticate, after removing previous auth token
                continue
            else:
                exception_content = 'Response error ' + str(response.status_code)
                try:
                    if response.json():
                        exception_content += " :: " + str(response.json())
                except json.decoder.JSONDecodeError:
                    exception_content += " :: " + str(response.content)

                raise MustlinkException(exception_content)

    def __prepare_url(self, endpoint):
        all_segments = []

        def _as_list(x):
            return x if isinstance(x, (list, tuple)) else (x,)

        for segment in _as_list(endpoint):
            all_segments.extend(_as_list(segment))

        endpoint_path = re.sub("/+", '/', "/".join([str(e) for e in all_segments]))

        return "{base_url}{endpoint}".format(base_url=self.base_url, endpoint=endpoint_path)

    @staticmethod
    def __prepare_params(parameters):
        if parameters:
            params = {}
            for k, v in parameters.items():
                if v:
                    if type(v) == datetime:
                        params[k] = v.strftime(DT_FMT)
                    else:
                        params[k] = v

            return params

    def __authenticate(self, force=False):
        auth_header = 'Authorization'
        if force:
            del self.headers[auth_header]

        if auth_header not in self.headers:
            response = requests.post(
                self.__prepare_url('/auth/login'),
                data=json.dumps({'username': self.username, 'password': self.password}),
                headers=self.headers,
                timeout=30)
            if response.status_code == 200:
                self.headers[auth_header] = response.json()['token']
            else:
                raise RuntimeError("Authentication failed")


class Users:
    """
    Expose all operations available for users
    """

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def search(self, name=None, group=None, email=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL
        group_re = re.compile(group, re.IGNORECASE) if group else RE_ALL

        def check(user_model):
            name_ok = name_re.search(user_model['login'] + ' ' + user_model['name']) is not None
            group_ok = False
            for g in user_model['groups']:
                group_ok = group_re.search(g['name']) is not None
                if group_ok:
                    break

            email_ok = not email or user_model['email'].lower() == email.lower()

            return name_ok and group_ok and email_ok

        return [u for u in self.ml.get(USER_URL) if check(u)]

    def load(self, user_login):
        try:
            return self.ml.get(USER_URL, 'id', int(user_login))
        except ValueError:
            return self.ml.get(USER_URL, 'login', user_login)
        except MustlinkException as e:
            log.error(e)
            return None

    def update(self, user_model):
        self.ml.put(USER_URL, user_model['login'], body=user_model)

    def change_password(self, login, new_password):
        user = self.load(login)
        if not user:
            raise MustlinkException("User '%s' not found" % login)

        user['password'] = new_password
        self.ml.put(USER_URL, "changepassword", body=user)

    def create(self, user_model):
        self.ml.post(USER_URL, body=user_model)

    def delete(self, user_login):
        self.ml.delete(USER_URL, user_login)

    def whoami(self):
        return self.ml.get(USER_MANAGEMENT_URL, 'userinfo')


class Groups:
    """
    Expose all operations available for groups
    """

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def search(self, name=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL
        return [g for g in self.ml.get(GROUPS_URL) if name_re.search(g['name'] + ' ' + g['description'])]

    def load(self, group):
        try:
            try:
                return self.ml.get(GROUPS_URL, 'id', int(group))
            except ValueError:
                return self.ml.get(GROUPS_URL, 'name', group)
        except MustlinkException as e:
            log.error(e)
            return None

    def update(self, group):
        self.ml.put(GROUP_URL, 'id', group['id'], body=group)

    def add_user(self, group, user):
        group_name = group['name'] if type(group) == dict else group
        user_login = user['login'] if type(user) == dict else user

        self.ml.put(GROUP_URL, 'name', group_name, 'add/user', user_login)

    def remove_user(self, group, user):
        group_name = group['name'] if type(group) == dict else group
        user_login = user['login'] if type(user) == dict else user

        self.ml.put(GROUP_URL, 'name', group_name, 'remove/user', user_login)

    def add_project(self, group, project):
        group_id = group['id'] if type(group) == dict else group
        project_id = project['id'] if type(project) == dict else project

        self.ml.put(GROUP_URL, 'id', group_id, 'add/project', project_id)

    def remove_project(self, group, project):
        group_id = group['id'] if type(group) == dict else group
        project_id = project['id'] if type(project) == dict else project

        self.ml.put(GROUP_URL, 'id', group_id, 'remove/project', project_id)

    def set_projects(self, group, projects):
        group_id = group['id'] if type(group) == dict else group
        projects_id = [p['id'] if type(p) == dict else int(p) for p in projects]
        self.ml.put(GROUP_URL, 'id', group_id, 'set/projects', projects_id)

    def add_role(self, group, role):
        group_id = group['id'] if type(group) == dict else group
        role_id = role['id'] if type(role) == dict else role

        self.ml.put(GROUP_URL, 'id', group_id, 'add/role', role_id)

    def remove_role(self, group, role):
        group_id = group['id'] if type(group) == dict else group
        role_id = role['id'] if type(role) == dict else role

        self.ml.put(GROUP_URL, 'id', group_id, 'remove/role', role_id)

    def set_roles(self, group, roles):
        group_id = group['id'] if type(group) == dict else group
        roles_id = [r['id'] if type(r) == dict else int(r) for r in roles]
        self.ml.put(GROUP_URL, 'id', group_id, 'set/projects', roles_id)


class Projects:
    """
    Expose all operations available for projects
    """

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def search(self, name=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL
        return [p for p in self.ml.get(PROJECTS_URL)
                if name_re.search(p['name'] + ' ' + p['description'])]

    def load(self, project):
        # there is no 'get' endpoint, thus it's simulating a precise load
        try:
            pid = int(project)
            _projects = [p for p in self.ml.get(PROJECTS_URL) if p['id'] == pid]
        except ValueError:
            _projects = [p for p in self.ml.get(PROJECTS_URL) if p['name'].lower() == project]

        if len(_projects) == 1:
            return _projects[0]
        else:
            return None


class Roles:
    """
    Expose all operations available for roles
    """

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def search(self, name=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL
        return [p for p in self.ml.get(ROLES_URL)
                if name_re.search(p['name'] + ' ' + p['description'])]

    def load(self, role):
        # there is no 'get' endpoint, thus it's simulating a precise load
        try:
            pid = int(role)
            _roles = [p for p in self.ml.get(ROLES_URL) if p['id'] == pid]
        except ValueError:
            _roles = [p for p in self.ml.get(ROLES_URL) if p['name'].lower() == role]

        if len(_roles) == 1:
            return _roles[0]
        else:
            return None


class DataProvider:

    def __init__(self, mustlink: Mustlink, data_provider):
        self.ml = mustlink
        self.data_provider = data_provider
        self.dp_segments = (DATA_PROVIDERS_URL, self.data_provider)

    def get_metadata(self):
        return self.ml.get(self.dp_segments)

    def get_tables(self):
        return self.ml.get(self.dp_segments, 'tables')

    def get_table_metadata(self, table):
        return self.ml.get(self.dp_segments, 'table', table, 'metadata')

    def get_table_data(self, table, start, end, mode='BRIEF',
                       max_rows=5000, filter_keys=(' ',), filter_values=(' ',)):
        params = {
            'dateFormat': 'fromTo',
            'from': start,
            'to': end,
            'mode': mode,
            'representation': 'SIMPLE',
            'maxRows': max_rows,
            'filterKeys': ",".join(filter_keys),
            'filterValues': ",".join(filter_values),
        }
        return self.ml.get(self.dp_segments, 'table', table, 'data', parameters=params)

    def get_aggregation_periods(self):
        return self.ml.get(self.dp_segments, 'aggregations')['aggs']

    def get_aggregation_functions(self):
        return self.ml.get(self.dp_segments, 'aggregations')['aggFuncs']

    def get_statistics(self, param_names, start, end,
                       aggregation_function=None,
                       aggregation_type=None,
                       aggregation_value=None):
        params = {
            'key': 'name',
            'values': ",".join(param_names),
            'from': start,
            'to': end,
            'aggregationFunction': aggregation_function,
            'aggregation': aggregation_type,
            'aggregationValue': aggregation_value
        }

        return self.ml.get(self.dp_segments, 'parameters/statistics', parameters=params)


    def get_parameters(self, key=None, value=None, parameter_type=None):
        params = {
            'key' : key,
            'value': value,
            'mode': 'SIMPLE',
            'search': True,
            'parameterType': parameter_type
        }
        return self.ml.get(self.dp_segments, 'parameters', parameters=params)

    def get_parameter_details(self, param_name, nature='SCPN'):
        params = {
            'key': 'name',
            'value': param_name,
            'mode': 'SIMPLE',
            'type': nature,
            'search': False
        }
        return self.ml.get(self.dp_segments, 'parameters', parameters=params)





    def get_parameters_data(self, param_names, start, end,
                            calibrate=None,
                            aggregation_function=None, aggregation_type=None, aggregation_value=None,
                            compression_error=None,
                            chunk_count=None):
        params = {
            'key': 'name',
            'values': ",".join(param_names),
            'from': start,
            'to': end,
            'calibrate': calibrate,
            'aggregationFunction': aggregation_function,
            'aggregation': aggregation_type,
            'aggregationValue': aggregation_value,
            'compressionError': compression_error,
            'chunkCount': chunk_count,
            'mode': 'SIMPLE'
        }

        return self.ml.get(self.dp_segments, 'parameters/data', parameters=params)

    def get_metadata_tree(self, param_id):
        return self.ml.get('/metadata/tree', self.data_provider, parameters={'id': param_id})


class DataProviders:
    """
    Wrapping class for data provider MUSTLINK operations
    """

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def search(self, name=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL
        return [dp for dp in self.ml.get(DATA_PROVIDERS_URL) if name_re.search(dp['name'] + ' ' + dp['description'])]

    def load_data_provider(self, name):
        dp = [dp['name'] for dp in self.ml.get(DATA_PROVIDERS_URL) if name.lower() == dp['name'].lower()]
        if len(dp) == 1:
            return DataProvider(self.ml, dp[0])
        else:
            return None

    def tree_search(self, fields: list, text, data_providers: list):
        # Returns tree nodes which metadata field {field} match the {text} of the given dataprovider(s).
        # Example: field=id,name,description&text=AST*&dataproviders=Cryosat2
        p = {
            'field': ",".join(fields),
            'text': text,
            'dataproviders': ",".join(data_providers)
        }
        return self.ml.get('/metadata/treesearch', parameters=p)


class Exports:
    HEADER_FORMATS = ['Grains', 'Simple', 'tdrs', 'darc']
    DATE_FORMATS = ['Normal', 'DOY', 'CCSDS', 'TDRS', 'S2K DOY']
    SEPARATORS = ['tab', 'comma']

    def __init__(self, mustlink: Mustlink):
        self.ml = mustlink

    def get_parameters_data_chart(self, data_provider, param_names, start, end,
                                  chart_format='png',
                                  aggregation_function=None, aggregation_type=None, aggregation_value=None,
                                  compression_error=None,
                                  chunk_count=None,
                                  detect_gaps=False):
        params = {
            'key': 'name',
            'values': ",".join(param_names),
            'from': start,
            'to': end,
            'aggregationFunction': aggregation_function,
            'aggregation': aggregation_type,
            'aggregationValue': aggregation_value,
            'compressionError': compression_error,
            'chunkCount': chunk_count,
            'detectGaps': detect_gaps,
        }

        return self.ml.get(DATA_PROVIDERS_URL, data_provider, 'chart', chart_format, parameters=params)

    def get_parameters_data_export(self, data_provider, param_names,
                                   start, end, last=None,
                                   # period???, type????
                                   prefix=None, header_format='Grains',
                                   date_format='Normal', use_millis=True,
                                   aggregation_function=None, aggregation_type=None, aggregation_value=None,
                                   separator='tab', combine_files=True,
                                   ):
        assert not header_format or header_format in self.HEADER_FORMATS, \
            "Header format '%s' not valid (must be one of: %s)" % (header_format, self.HEADER_FORMATS)

        assert date_format in self.DATE_FORMATS, \
            "Date format '%s' not valid (must be one of: %s)" % (date_format, self.DATE_FORMATS)

        assert separator in self.SEPARATORS, \
            "Separator '%s' not valid (must be one of: %s)" % (separator, self.SEPARATORS)

        params = {
            'key': 'name',
            'values': ",".join(param_names),
            'type': None,  # TODO ???
            'from': start,
            'to': end,
            'last': int(last) if last else None,
            'period': None,  # TODO ???
            'aggregationFunction': aggregation_function,
            'aggregation': aggregation_type,
            'aggregationValue': aggregation_value,
            'prefix': prefix,
            'headerFormat': header_format,
            'dateFormat': "%s-millis" % date_format if use_millis else date_format,
            'separator': separator,
            'combine': combine_files,
        }

        return self.ml.get(DATA_PROVIDERS_URL, data_provider, 'parameters/data/export', parameters=params)

    def get_table_data_export(self, data_provider, table, start, end,
                              mode='BRIEF', last=None,  # period=???, rows_threshold=???, token=???,
                              filter_keys=(' ',), filter_values=(' ',)):
        params = {
            'dateFormat': 'fromTo',
            'from': start,
            'to': end,
            'last': int(last) if last else None,
            'period': None,  # TODO ???
            'mode': mode,
            'representation': 'SIMPLE',
            'filterKeys': ",".join(filter_keys),
            'filterValues': ",".join(filter_values),
            'rowsThreshold': None,  # TODO ???
            'token': None,  # TODO ???
        }
        return self.ml.get(DATA_PROVIDERS_URL, data_provider, 'table', table, 'export', parameters=params)


class DrMust:
    def __init__(self, ml: Mustlink):
        self.ml = ml

    def get_tasks(self):
        # return self.ml.get(DRMUST_URL)
        raise NotImplemented

    def create_task(self):
        # self.ml.post(DRMUST_URL, body=task)
        raise NotImplemented

    def count_tasks(self):
        # return self.ml.get(DRMUST_URL, 'numresults')["number"]
        raise NotImplemented

    def get_task(self, task_id):
        # return self.ml.get(DRMUST_URL, task_id)
        raise NotImplemented

    def delete_task(self, task_id):
        # self.ml.delete(DRMUST_URL, task_id)
        raise NotImplemented


class Novelty:
    def __init__(self, ml: Mustlink):
        self.ml = ml

    def get_missions(self):
        # return self.ml.get(NOVELTY_URL, 'missions')
        raise NotImplemented

    def get_novelties(self, data_provider, novelty_type, start_date, end_date, last=None, period=None):
        # self.ml.get('analysis', data_provider, 'novelty', parameters= {
        #     'type': novelty_type,
        #     'from': start_date,
        #     'to': end_date,
        #     'last': last,
        #     'period': period,
        # })
        raise NotImplemented


class Scripts:
    """
    Scripts object wrapper
    TODO: complete the implementation when the API URLs have been reorganized
    """
    def __init__(self, ml: Mustlink):
        self.ml = ml

    @staticmethod
    def load_code(script_code):
        _code = None
        try:
            with open(script_code, 'r') as script_file:
                _code = "".join(script_file.readlines())
        except FileNotFoundError:
            _code = script_code

        return _code

    def retrieve_user_id(self, owner):
        try:
            return int(owner)
        except ValueError:
            try:
                return int(owner['id'])
            except:
                return int(self.ml.users.search(name=owner)[0]['userId'])

    def search(self, name=None, user=None, is_public=None):
        name_re = re.compile(name, re.IGNORECASE) if name else RE_ALL

        def is_ok(script):
            name_filter = name_re.search(script['name'] + ' ' + script['description'])
            if user is None:
                user_id_filter = owner_filter = True
            else:
                try:
                    user_id_filter = script['userId'] == int(user)
                    owner_filter = True
                except ValueError:
                    user_id_filter = True
                    owner_filter = re.search(user, script['ownerLogin'], re.IGNORECASE)

            id_public_filter = is_public is None \
                               or (script['publicFlag'] == 0 and is_public is False) \
                               or (script['publicFlag'] == 1 and is_public is True)

            return name_filter and user_id_filter and owner_filter and id_public_filter

        return [s for s in self.ml.get(SCRIPTS_URL) if is_ok(s)]


def get(mustlink=None) -> Mustlink:
    """
    Supply or create Mustlink from the parameter passed or returns the default Mustlink instance
    :param mustlink: configuration dictionary or instance of Mustlink object
    :return:
    """
    global __mustlink
    ml = None

    if mustlink:
        if type(mustlink) == dict:
            ml = Mustlink(**mustlink)
        elif type(mustlink) == Mustlink:
            ml = mustlink
    else:
        ml = __mustlink

    if ml is None:
        raise ValueError("No Mustlink initialized")

    return ml


def init_defaults(mustlink):
    if type(mustlink) == Mustlink:
        ml = mustlink
    else:
        ml = Mustlink(**mustlink)

    global __mustlink
    __mustlink = ml
    log.debug("pynk initialized with instance %s@%s" % (__mustlink.username, __mustlink.base_url))
