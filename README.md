# Yule-Walker method for IIR filter design given an arbritrary magnitude response

## Description

This implements a modified yule-walker method for estimating the coefficients of an IIR filter that best matches an arbritrary frequency response. This is a port
from sci-lab's version here: [yulewalk.sci](https://github.com/scilab/scilab/blob/master/scilab/modules/signal_processing/macros/yulewalk.sci) and is distributed by
the same license ([CeCILLv2.0](https://opensource.org/licenses/CECILL-2.1)).

## Install

### Using pip

```
$ pip install yulewalker
```

### By cloning this repository

```
$ git clone https://github.com/mmxgn/yulewalker.git
$ cd yulewalker
$ python setup.py install
```

## Testing

Using `doctest`:
```
$ python -m yulewalker.yulewalker -v
```

Using `pytest`:
```
$ pytest
```

## Example

See [Example.ipynb](Example.ipynb)


## Citation

If you are using this code in your research please cite the original paper here:
```
Friedlander, Benjamin, and Boaz Porat. "The modified Yule-Walker method of ARMA spectral stimation." IEEE Transactions on Aerospace and Electronic Systems 2 (1984): 158-173.
```

## License

See [LICENSE](LICENSE)


## Acknowledgements

Thanks to [Nonomo](http://www.nomono.co) for backing this library.
