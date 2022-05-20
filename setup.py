from setuptools import setup
import os


def get_version():
    topdir = os.path.abspath(os.path.join(__file__, '..'))
    with open(os.path.join(topdir, 'vbt3', '__init__.py'), 'r') as f:
        for line in f.readlines():
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
    raise ValueError("Version string not found")

VERSION = get_version()

setup(name='vbt3',
      version=VERSION,
      description='Symbolic calculations for semi-quantitative Valence Bond Theory',
      url='http://github.com/talipovm/vbt3',
      author='Marat Talipov',
      author_email='talipovm@nmsu.edu',
      license='MIT',
      packages=['vbt3'],
      zip_safe=True,
      install_requires=['sympy>=1.8'],
      )