[![Build Status](https://github.com/chatziko/libqif/workflows/build/badge.svg)](https://github.com/chatziko/libqif/actions)
[![PyPI version](https://badge.fury.io/py/qif.svg)](https://badge.fury.io/py/qif)

# libqif

## Install

```bash
pip install qif
```

- Needs `python` >= 3.8 and a `sandybridge` or later CPU
- On linux `pip` >= 19 is needed (make sure to `pip install -U pip`)


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
qif.set_default_type(qif.rat)

C = np.array([
	[rat(1,2), rat(1,4), rat(1,4)],
	[rat(1,6), rat(3,6), rat(2,6)],
	[rat(1,2), rat(1,2), rat(0)],
])
compute_bayes(C)
```

## Documentation

A list of methods provided by `qif` is available [here](http://chatziko.github.io/libqif/).

## Use libqif with C++

See the [installation instructions](INSTALL.md).