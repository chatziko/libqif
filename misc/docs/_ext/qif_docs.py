"""
Sphinx extension for the qif docs:

- Signatures are taken from the .pyi stubs. The pybind11 signatures are hard to read (C++ types,
  "Overloaded function" listings for the double/rat versions), the stubs are cleaner.
- Warns about qif functions that are not documented in any page.
"""
import ast
import os
import re
import types

from sphinx.util import logging

import qif

logger = logging.getLogger(__name__)

_stub_cache = {}


def _stubs_dir():
	# the stubs are next to the module in an installed wheel, but not in the cmake build dir
	return os.environ.get('QIF_STUBS_DIR') or os.path.dirname(os.path.dirname(qif.__file__))


def _stub_functions(module):
	"""Returns {name: [FunctionDef, ...]} for the stub of the given module (eg 'qif.measure.bayes_vuln')."""
	if module not in _stub_cache:
		base = os.path.join(_stubs_dir(), *module.split('.'))
		_stub_cache[module] = {}
		for path in [base + '.pyi', os.path.join(base, '__init__.pyi')]:
			if os.path.exists(path):
				with open(path) as f:
					tree = ast.parse(f.read())
				for node in tree.body:
					if isinstance(node, ast.FunctionDef):
						_stub_cache[module].setdefault(node.name, []).append(node)
				break
	return _stub_cache[module]


def _unparse(node):
	return re.sub(r'\bt\.', '', ast.unparse(node))		# t.ndarray -> ndarray


def process_signature(app, what, name, obj, options, signature, return_annotation):
	if what != 'function':
		return None
	module, _, func = name.rpartition('.')
	defs = _stub_functions(module).get(func)
	if not defs:
		logger.warning(f'qif_docs: no stub found for {name}')
		return None

	# overloads are returned one per line (the 2nd onwards prefixed with the name), autodoc renders each as a separate signature
	sigs = []
	for i, d in enumerate(defs):
		sig = (func if i else '') + f'({_unparse(d.args)})'
		if d.returns is not None:
			sig += f' -> {_unparse(d.returns)}'
		sigs.append(sig)
	return '\n'.join(sigs), None


def process_docstring(app, what, name, obj, options, lines):
	# remove pybind11's signature listing, keep anything else (the actual docs are in the .md pages)
	func = name.rpartition('.')[2]
	sig_re = re.compile(rf'^(\d+\. )?{re.escape(func)}\(')
	lines[:] = [l for l in lines if l != 'Overloaded function.' and not sig_re.match(l)]
	while lines and not lines[0].strip():
		lines.pop(0)


def _public_functions(module, prefix, seen):
	for attr in dir(module):
		if attr.startswith('_'):
			continue
		obj = getattr(module, attr)
		if isinstance(obj, types.ModuleType):
			if id(obj) not in seen and obj.__name__.lstrip('_').startswith('qif.'):
				seen.add(id(obj))
				yield from _public_functions(obj, f'{prefix}.{attr}', seen)
		elif callable(obj) and type(obj).__name__ == 'builtin_function_or_method':
			yield f'{prefix}.{attr}'


def check_undocumented(app, env):
	documented = env.get_domain('py').objects
	for name in sorted(set(_public_functions(qif, 'qif', set()))):
		if name not in documented and not name.startswith('qif.typing.'):
			logger.warning(f'qif_docs: {name} is not documented')


def setup(app):
	app.connect('autodoc-process-signature', process_signature)
	app.connect('autodoc-process-docstring', process_docstring)
	app.connect('env-check-consistency', check_undocumented)
	return {'parallel_read_safe': True}
