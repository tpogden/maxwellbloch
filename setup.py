from setuptools import setup

setup(name='MaxwellBloch',
      version='0.3.0.dev',
      description='A Python package for solving the Maxwell-Bloch equations.',
      url='http://github.com/tommyogden/maxwellbloch',
      author='Thomas P Ogden',
      author_email='t@ogden.eu',
      license='MIT',
      packages=['maxwellbloch'],
      install_requires=['qutip'],
      scripts=['bin/make-mp4-fixed-frame.py',
               'bin/make-gif-ffmpeg.sh'],
      zip_safe=False)
