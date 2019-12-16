import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
	def __init__(self, name, sourcedir=''):
		Extension.__init__(self, name, sources=[])
		self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
	def run(self):
		try:
			out = subprocess.check_output(['cmake', '--version'])
		except OSError:
			raise RuntimeError("CMake must be installed to build the following extensions: " +
							   ", ".join(e.name for e in self.extensions))

		if platform.system() == "Windows":
			cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
			if cmake_version < '3.1.0':
				raise RuntimeError("CMake >= 3.1.0 is required on Windows")

		for ext in self.extensions:
			self.build_extension(ext)

	def build_extension(self, ext):
		extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
		cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
					  '-DPYTHON_EXECUTABLE=' + sys.executable]

		cfg = 'Debug' if self.debug else 'Release'
		build_args = ['--config', cfg]

		if platform.system() == "Windows":
			cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
			if sys.maxsize > 2**32:
				cmake_args += ['-A', 'x64']
			build_args += ['--', '/m']
		else:
			cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
			build_args += ['--', '-j2']

		env = os.environ.copy()
		env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
															  self.distribution.get_version())
		if not os.path.exists(self.build_temp):
			os.makedirs(self.build_temp)
		subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
		subprocess.check_call(['cmake', '--build', '.', '-t', os.path.basename(ext.name)] + build_args, cwd=self.build_temp)


# if the current commit matches a vN.N.N tag, use that as the version
tag = subprocess.check_output('git describe --tags --exact-match; exit 0', stderr=subprocess.DEVNULL, shell=True).decode("utf-8")
m = re.match(r"v([\d.]+)", tag)
version = m.group(1) if m else '0.0.1'

# read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, '../README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
	name='qif',
	version=version,
	author='Kostas Chatzikokolakis  ',
	author_email='kostas@chatzi.org',
	url='https://github.com/chatziko/libqif',
	description='Quantitative Information Flow library',
	long_description=long_description,
    long_description_content_type='text/markdown',
	# packages=[''],
	install_requires=[
		'numpy',
	],
	ext_modules=[CMakeExtension('qif_module')],
	cmdclass=dict(build_ext=CMakeBuild),
	zip_safe=False,
	# test_suite='tests/test.py',
)