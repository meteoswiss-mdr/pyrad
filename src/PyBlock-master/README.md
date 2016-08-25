PyBlock README
--------------
PyBlock is a Python 2 or 3 module that enables the end user to estimate partial beam blockage using polarimetric radar data. The methodologies it uses depend on the self-consistency of polarimetric radar variables - reflectivity (Zh), differential reflectivity (Zdr), and specific differential phase (Kdp) - in pure rain. There are two methodologies currently available to the end user, both described in Lang et al. (2009): The KDP method, and the Fully Self-Consistent (FSC) method. Briefly, the KDP method will check the behavior of Zh and Zdr for a given range of Kdp both inside and outside of blocked azimuths, and use that to suggest corrections to these measurands. This is effectively a relative calibration of Z and Zdr. The FSC method uses a derived or specified self-consistency relationship to do an absolute calibration of Zh within the blocked regions. PyBlock implements these methodologies within an object-oriented Python framework. 

PyBlock Installation
-------------------
The following dependencies need to be installed first:

- A robust version of Python 2.7 or 3.4 (other versions untested) w/ most standard scientific packages (e.g., numpy, matplotlib, pandas, etc.) - Get one for free [here.](https://store.continuum.io/cshop/anaconda/)
- [The Python Atmospheric Radiation Measurement (ARM) Radar Toolkit (Py-ART)] (https://github.com/ARM-DOE/pyart)
- [CSU_RadarTools](https://github.com/CSU-Radarmet/CSU_RadarTools)
- [SkewT](https://pypi.python.org/pypi/SkewT) - Python 3 version can be found [here.](https://github.com/tjlang/SkewT)
- [DualPol] (https://github.com/nasa/DualPol)

Specific import calls in the PyBlock source code:

```
from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from warnings import warn
import statsmodels.api as sm
import os
import gzip
import pickle
import pyart
import dualpol
from csu_radartools import csu_misc
import six
```

To install PyBlock, in the main directory for the package:

```
python setup.py install
```

Using PyBlock
-------------
To access everything:
```
import pyblock
```

To see PyBlock in action, check out the IPython notebook provided in this distribution.
