# Installing MaxwellBloch

I recommend using Conda environments. 

You can create and activate an environment named `mb` with all the required
dependencies for MaxwellBloch with
```sh
conda create --name mb python=3 numpy=1 scipy=1 qutip=4
conda activate mb
```
The MaxwellBloch package can then be installed from
[PyPI](https://pypi.org/project/MaxwellBloch/) using
```
pip install maxwellbloch
```

Alternatively you can manually download and install any [
release](https://github.com/tpogden/maxwellbloch/releases).

<!-- 
## Testing your Installation

To test MaxwellBloch and its dependencies are working correctly, install `pytest`
```
conda install pytest=3
```
and then call
```
python -c 'from maxwellbloch import testing; testing.run()'
``` -->

