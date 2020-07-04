.. MaxwellBloch documentation master file, created by
   sphinx-quickstart on Thu Aug  8 20:52:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MaxwellBloch Documentation
==========================

MaxwellBloch is a Python package for solving the coupled Maxwell-Bloch equations
describing the nonlinear propagation of near-resonant light through thermal
quantised systems such as atomic vapours.

.. .. image:: example.gif

..    :scale: 100 %

.. figure:: example.gif
   :alt: Propagation of a 4π pulse through a dense atomic vapour.

   Propagation of a 4π pulse through a dense atomic vapour. The pulse
   immediately breaks up on entering the medium and the resultant pulses form
   two optical solitons each with a pulse area of 2π.

MaxwellBloch is used for theoretical research, for modelling experiments and for
undergraduate and graduate teaching.

MaxwellBloch can model two, three or many-level systems with physical
effects including:

* Inhomogeneous broadening due to spontaneous decay,
* Doppler broadening in thermal systems,
* Collision dephasing,
* Sub-level structure.

Modules are also available for:

* Generating hyperfine structure for alkali atoms with the correct channels and 
  angular momentum factors for coupling and decay,
* Specral analysis,
* Plotting and animating solutions.

Some phenomena that can be demonstrated:

* Linear absorption and dispersion,
* Fast light,
* Slow light,
* Electromagnetically Induced Transparency (EIT),
* Storage and Retrieval of Light pulses,
* Self-Induced Transparency (SIT) and Optical Solitons,
* Matched Pulses and Simultons,
* Hyperfine Pumping.

See the Examples section below for details.

.. toctree::
   :maxdepth: 2
   :caption: Installation

   install

.. toctree::
   :maxdepth: 2
   :caption: Usage
   
   usage/two-level
   usage/three-level
   usage/structure
   usage/spectral-analysis
   usage/velocity-classes
   usage/built-in-time-functions
   usage/scripts

.. toctree::
   :maxdepth: 3
   :caption: Examples

   examples/two-level
   examples/three-level
   examples/structure
   
.. toctree::
   :maxdepth: 2
   :caption: Support

   support/troubleshooting

.. Indices and Tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
