import collections

import pynk


def usermod(login, email=None, name=None, password=None, add_groups=None, remove_groups=None, mustlink=None):
    ml = pynk.get(mustlink)
    user = ml.users.load(login)

    if not user:
        raise ValueError("User %s not found" % login)

    if add_groups or remove_groups:
        _all_groups = add_groups if add_groups else [] + remove_groups if remove_groups else []
        for group in [ml.groups.load(g) for g in _all_groups]:
            if not group:
                raise ValueError("Group %s not found", group)

            if add_groups:
                ml.groups.add_user(group, user)
            elif remove_groups:
                ml.groups.remove_user(group, user)

    if email or name:
        if name:
            user['name'] = name
        if email:
            user['email'] = email

        ml.users.update(user)

    if password:
        ml.users.change_password(user['login'], password)


def useradd(login, name, email, password, groups, force_email=False, mustlink=None):
    ml = pynk.get(mustlink)

    if ml.users.load(login):
        raise ValueError("User %s already exists" % login)

    if len(ml.users.search(email=email)) > 0 and force_email is False:
        raise ValueError("Users with email %s already exist" % email)

    user_groups = []
    for g in groups:
        _g = ml.groups.load(g)
        if _g:
            user_groups.append(_g)
        else:
            raise ValueError("Group %s doesn't exist" % g)

    new_user = {
        "login": login,
        "name": name,
        "email": email,
        "password": password,
        "roles": [],
        "groups": user_groups,
        "projects": []
    }
    # the user creation end point doesn't link the user to the groups
    ml.users.create(new_user)

    # thus provide linking groups and users
    for group in user_groups:
        ml.groups.add_user(group, new_user)


def userdel(login, mustlink=None):
    ml = pynk.get(mustlink)
    if ml.users.load(login) is None:
        raise ValueError("User %s doesn't exists" % login)

    ml.users.delete(login)


def groupmod(group_name,
             description=None,
             add_users=None, remove_users=None,
             add_projects=None, remove_projects=None,
             add_roles=None, remove_roles=None,
             mustlink=None):

    ml = pynk.get(mustlink)
    group = ml.groups.load(group_name)
    if not group:
        raise ValueError("Group %s doesn't exist" % group_name)

    if description:
        group['description'] = description
        ml.groups.update(group)

    def do_link(entities, manager, subject, func):
        if entities:
            _list = []
            for _e, _entity in ((e, manager.load(e)) for e in entities):
                if not _entity:
                    raise ValueError("%s %s doesn't exist" % (subject, _e))
                _list.append(_entity)
    
            for _e in _list:
                func(group, _e)

    do_link(add_users, ml.users, 'User', ml.groups.add_user)
    do_link(remove_users, ml.users, 'User', ml.groups.remove_user)
        
    do_link(add_projects, ml.projects, 'Project', ml.groups.add_project)
    do_link(remove_projects, ml.projects, 'Project', ml.groups.remove_project)
    
    do_link(add_roles, ml.roles, 'Role', ml.groups.add_role)
    do_link(remove_roles, ml.roles, 'Role', ml.groups.remove_role)


def list_users(name=None, group=None, email=None, mustlink=None):
    ml = pynk.get(mustlink)

    UserTuple = collections.namedtuple('UserTuple', 'id,login,email,name,groups,projects,roles')

    table = []
    for user in sorted(ml.users.search(name, group, email), key=lambda x: x['login']):
        table.append(UserTuple(
            user['id'],
            user['login'],
            user['email'],
            user['name'],
            [g['name'] for g in user['groups']],
            [p['name'] for p in user['projects']],
            [r['name'] for r in user['roles']]
        ))

    return table


def list_groups(name=None, mustlink=None):
    ml = pynk.get(mustlink)

    GroupTuple = collections.namedtuple('GroupTuple', 'id,name,description,projects,roles')

    table = []
    for group in sorted(ml.groups.search(name), key=lambda g: g['name'].lower()):
        table.append(GroupTuple(
            group['id'],
            group['name'],
            group['description'],
            [p['name'] for p in group['projects']],
            [r['name'] for r in group['roles']],
        ))

    return table


def list_projects(name=None, mustlink=None):
    ml = pynk.get(mustlink)

    ProjectTuple = collections.namedtuple('ProjectTuple', 'id,name,description')
    table = []

    for project in sorted(ml.projects.search(name), key=lambda p: p['name'].lower()):
        table.append(ProjectTuple(
            project['id'],
            project['name'],
            project['description']
        ))

    return table


def check_user(login, mustlink=None):
    user = pynk.get(mustlink).users.load(login)
    if not user:
        raise ValueError("User %s doesn't exists" % login)
