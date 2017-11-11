from setuptools import setup

setup(name='MaxwellBloch',
      version='0.2.0',
      description='A Python package for solving the Maxwell-Bloch equations.',
      url='http://github.com/tommyogden/maxwellbloch',
      author='Thomas P Ogden',
      author_email='t@ogden.eu',
      license='MIT',
      packages=['maxwellbloch'],
      install_requires=[
          'qutip',
          'tqdm',
      ],
      zip_safe=False)
