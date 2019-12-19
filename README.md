# libqif

[![Build Status](https://github.com/chatziko/libqif/workflows/build-pypi/badge.svg)](https://github.com/chatziko/libqif/actions)

## Install

```bash
pip install qif
```

- Needs `python` >= 3.5 and a `sandybridge` or later CPU
- On linux `pip` >= 19 is needed (make sure to `pip install -U pip`)


## A sample program

```python
from qif import *

C = channel.randu(5)
pi = probab.uniform(5)

print("Channel:\n", C)
print("Prior:\n", pi)

print("Bayes vulnerability", measure.bayes_vuln.posterior(pi, C))
print("Bayes mult-capacity", measure.bayes_vuln.mult_capacity(C))
```

## Documentation

A list of available methods is available via `help`:

```python
import qif

help(qif)
help(qif.channel)
...
```

## Use libqif with C++

See the [installation instructions](INSTALL.md).