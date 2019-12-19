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
pi = probab.uniform(3)

print(C)
print(pi)

print("Bayes vulnerability", measure.bayes_vuln.posterior(pi, C))
print("Bayes capacity", measure.bayes_vuln.capacity(C))
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