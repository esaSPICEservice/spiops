import shutil


def update_path(files):
    for file in files:
        mk = open(file, 'r')
        new_mk = []
        for line in mk.readlines():
            if 'PATH_VALUES' in line:
                line = line.replace('..', '../..')
            new_mk.append(line)
        mk.close()
        mk = open(file, 'w')
        for line in new_mk:
            mk.write(line)
        mk.close()
    return


def move2formerv(files):
    for file in files:
        shutil.move(file, 'former_versions/' + file)
    return


def update_version(files, version):
    for file in files:
        shutil.copyfile(file, file.replace('.', '_' + version + '.'))
    return

