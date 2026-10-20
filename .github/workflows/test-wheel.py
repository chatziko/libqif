# Sample script testing an installed qif wheel. Run after each wheel build (build.yml) and on the wheels from PyPI (test-pypi.yml)

import os
import sys
import numpy
import qif
from qif import *

print(f"qif {qif.__version__}, python {sys.version.split()[0]}, numpy {numpy.__version__}")

# the sample from the README
def compute_bayes(C):
	pi = probab.uniform(C.shape[0])
	post = measure.bayes_vuln.posterior(pi, C)
	cap = measure.bayes_vuln.mult_capacity(C)
	print("Bayes vulnerability", post, "mult-capacity", cap)
	return post, cap

post, cap = compute_bayes(channel.randu(5))
assert 0 < post <= 1 and 1 <= cap <= 5

set_default_type(rat)
C = numpy.array([
	[rat(1,2), rat(1,4), rat(1,4)],
	[rat(1,6), rat(3,6), rat(2,6)],
	[rat(1,2), rat(1,2), rat(0)],
])
post, cap = compute_bayes(C)
assert post == rat(4,9) and cap == rat(4,3)	# column maxima: 1/2 + 1/2 + 1/3 = 4/3
set_default_type(double)

# factorize solves a linear program, use or-tools' GLOP (unless overridden)
lp.defaults.solver = os.environ.get("QIF_LP_SOLVER", "GLOP")
B = channel.randu(4)
A = B @ channel.randu(4)
X = channel.factorize(A, B)
assert X.size > 0 and numpy.allclose(B @ X, A), "factorize failed"
print(f"factorize with {lp.defaults.solver}: OK")

print("All OK")
