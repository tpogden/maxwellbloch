# MaxwellBloch

By [Thomas P Ogden](mailto:t@ogden.eu).

MaxwellBloch is a Python package for solving the coupled Maxwell-Bloch
equations describing the nonlinear propagation of near-resonant light through
thermal atomic vapours.

![](example.gif)

Above is an [example solution][4pi] for the propagation of a 4π pulse through a
dense atomic vapour. The pulse immediately breaks up on entering the medium and
the resultant pulses form two optical solitons each with a pulse
area of 2π.

[4pi]: https://github.com/tommyogden/notebooks-maxwellbloch/blob/master/examples/mb-solve-two-sech-4pi.ipynb

## Install

I recommend using Conda environments. MaxwellBloch requires NumPy and SciPy,
which can be installed via

    conda install numpy scipy

and QuTiP, which can be installed via

    conda install -c jrjohansson qutip=3.1.0

The MaxwellBloch package can then be installed using `pip`

    pip install maxwellbloch

 or manually download and install the [latest release](https://github.com/tommyogden/maxwellbloch/releases).

## Tutorial

A series of tutorial notebooks is in development here:

https://github.com/tommyogden/notebooks-maxwellbloch#tutorial

## Examples

If you prefer to learn by example, there are a large number of example
notebooks available here:

https://github.com/tommyogden/notebooks-maxwellbloch#examples

## Documentation

On the way…

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

MIT License. See [License.txt](LICENSE.txt).