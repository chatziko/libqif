# qif

`qif` is a library for Quantitative Information Flow: it computes vulnerability and leakage measures of
channels, constructs optimal privacy mechanisms, and checks refinement relations between channels. It is
written in C++ (libqif) with Python bindings; this documentation covers the Python package.

## Install

```bash
pip install qif
```

- Needs a `haswell` or later CPU (AVX2)
- Python 3.11 to 3.15
- MacOS 15 or later (x86_64 and arm64)

To use libqif from C++ see the [installation instructions](https://github.com/chatziko/libqif/blob/master/INSTALL.md).

## A sample program

```python
from qif import *

def compute_bayes(C):
    pi = probab.uniform(C.shape[0])

    print("Channel:\n", C)
    print("Prior:\n", pi)

    print("Bayes vulnerability", measure.bayes_vuln.posterior(pi, C))
    print("Bayes mult-capacity", measure.bayes_vuln.mult_capacity(C))

compute_bayes(channel.randu(5))

# same, but using rational arithmetic
set_default_type(rat)

C = numpy.array([
    [rat(1,2), rat(1,4), rat(1,4)],
    [rat(1,6), rat(3,6), rat(2,6)],
    [rat(1,2), rat(1,2), rat(0)],
])
compute_bayes(C)
```

```{toctree}
:caption: Modules
:maxdepth: 1

api/qif
api/channel/index
api/probab
api/metric/index
api/measure/index
api/mechanism/index
api/refinement
api/utility
api/lp
```
