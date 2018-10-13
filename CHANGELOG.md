# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic
Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

- Added script to make MP4s of systems with two fields

## [0.4.0] 2018-08-28

- Fixed the mp4 and gif filenames so they don't need to be .json.mp4.gif
- Changed t_funcs so that the _1, _2, … suffix is no longer required
- Changed ob_atom to atom
- Added label to atom
- Fixed interpolation on mp4 script
- Changed look of animated plots
- Fixed mp4 script to use real part of rabi freq
- Added automatic update of version number from setup.py
- Changed to allow empty velocity class to be set
- Changed the time function for fields to be constant if none is given.

## [0.3.0] 2018-01-08

- Added Travis CI build
- Fixed bug that tests submodule not in package
- Added git hash to unreleased version number
- Added scripts to make MP4s and gifs to /bin
- Added methods to save field CSV files
- Refactored spectral methods to separate module
- Refactored fixed frame methods to separate module
- Fixed bug where empty decay list would cause exception
- Fixed bug where fixed frame didn't work if no inner z steps
- Added testing.py to run tests
- Added link to video and gif tools
- Fixed multiple coupled levels bug
- Added MBSolve.set_field_rabi_freq_t_func() and set_field_rabi_freq_t_args()
    methods so we can add custom input fields
- Added OBAtom.build() method for resetting operators
- Added FFT methods for spectral analysis of pulse propagation
- Added ability to run `ob_solve` and `mb_solve` from the command line
- Added methods to get results in the fixed frame of reference
- Changed to QuTiP progress bar to remove a dependency

## [0.2.0] 2017-03-27

- Added ability to solve the Maxwell-Bloch equations
