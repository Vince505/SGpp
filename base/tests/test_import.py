# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

import sys

print "executing pysgpp import test (pysgpp/importTest.py)"
print "trying to import python lib (base) ...",

try:
    import base
except Exception as e:
    print "failed, error:"
    print str(e)
    sys.exit(1)
print "success"

# import os
# print os.environ['LD_LIBRARY_PATH']
# print os.environ['PYTHONPATH']

#content = dict([(name, cls) for name, cls in pysgpp.__dict__.items() if isinstance(cls, type)])