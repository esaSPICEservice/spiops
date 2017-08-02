import spiceypy as cspice

class Body(object):
    def __init__(self, body, time=object()):

        if isinstance(body, str):
            name = body
            id = cspice.bodn2c(body)
        else:
            id = body
            name = cspice.bodc2n(body)

        self.name = name
        self.id = id
        self.time = time


class Target(Body):
    def __init__(self, body, time=object()):

        super(Target, self).__init__(body, time=object())

        self.frame = 'IAU_{}'.format(self.name)
        self.__getRadii()

    def __getRadii(self):
        self.radii = cspice.bodvar(self.id, 'RADII', 3)



class Observer(Body):
    def __init__(self, body, time=object()):

        super(Observer, self).__init__(body, time=object())

        self.frame = '{}_SPACECRAFT'.format(self.name)
