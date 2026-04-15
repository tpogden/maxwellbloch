# Installing MaxwellBloch

## Using uv (recommended)

[uv](https://docs.astral.sh/uv/) is the recommended way to install MaxwellBloch:

```sh
uv pip install maxwellbloch
```

## Using pip

```sh
pip install maxwellbloch
```

## Using Conda

You can create and activate a Conda environment with the required dependencies:

```sh
conda create --name mb -c conda-forge python=3.11 qutip
conda activate mb
pip install maxwellbloch
```

Alternatively you can manually download and install any
[release](https://github.com/tpogden/maxwellbloch/releases).
