from setuptools import setup

setup(name='MaxwellBloch',
      version='0.3.0.dev',
      description='A Python package for solving the Maxwell-Bloch equations.',
      url='http://github.com/tommyogden/maxwellbloch',
      author='Thomas P Ogden',
      author_email='t@ogden.eu',
      license='MIT',
      packages=['maxwellbloch'],
      install_requires=[
          'qutip',
      ],
      zip_safe=False)
