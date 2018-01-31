# Force field initialization
#
# Written by Konrad Hinsen
#

"""
Force fields

"""

__docformat__ = 'restructuredtext'

import os, string, sys
import six

from pkg_resources import resource_stream

class _ForceFieldLoader(object):

    def __init__(self, module, object):
        self.module = module
        self.object = object
        self.globals = globals()

    __safe_for_unpickling__ = True

    def __call__(self, *args, **kw):
        ffc = getattr(__import__(self.module, self.globals), self.object)
        ffc.description = _description
        return apply(ffc, args, kw)

def _description(self):
    return 'ForceFields.' + self.__class__.__name__ + self.arguments

ff_list = resource_stream(__name__, 'force_fields').readlines()
for line in ff_list:
    if six.PY3:
        line = line.decode()
    line = (line).split()
    six.exec_(line[0] + "=_ForceFieldLoader(line[1], line[2])")

try:
    default_energy_threads = int(os.environ['MMTK_ENERGY_THREADS'])
except KeyError:
    default_energy_threads = 1

del os
del string
del sys
del ff_list
del line

