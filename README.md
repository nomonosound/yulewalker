# Yule-Walker methor for filter design for an arbritrary magnitude response

## Description

This implements a modified yule-walker method for estimating the coefficients of an IIR filter that best matches an arbritrary frequency response. This is a port
from sci-lab's version here: [yulewalk.sci](https://github.com/scilab/scilab/blob/master/scilab/modules/signal_processing/macros/yulewalk.sci) and is distributed by
the same license ([CeCILLv2.0](https://opensource.org/licenses/CECILL-2.1)).

## Install

Clone this repository, then run:

```
python setup.py install
```

## Testing

Using `doctest`:
```
python -m yulewalker.yulewalker -v
```

Using `pytest``:
```
pytest
```

## Example

See (Example.ipynb)[Example.ipynb]

## LICENSE

See (LICENSE)[LICENSE]