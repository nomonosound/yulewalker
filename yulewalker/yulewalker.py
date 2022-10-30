from typing import Tuple
import scipy
import scipy.linalg
import scipy.signal
import numpy as np


def complex_log(z: complex) -> complex:
    """
    Returns the complex log of `z`

    >>> complex_log(1)
    0j

    >>> complex_log(1j)
    1.5707963267948966j

    >>> complex_log(1+1j)
    (0.3465735902799727+0.7853981633974483j)
    """
    return np.log(np.abs(z)) + 1j * np.angle(z)


def mrdivide(B: np.ndarray, A: np.ndarray) -> np.ndarray:
    """
    Equivalent to the '/' operator or `mrdivide` in octave/matlab.
    Solves systems of linear equations of the form xA = B.

    See Nathan Pyle's answer here:
    https://stackoverflow.com/questions/1007442/mrdivide-function-in-matlab-what-is-it-doing-and-how-can-i-do-it-in-python


    >>> A = np.array([[1,1,3],[2,0,4],[-1,6,-1]]); B = np.array([[2,19,8]]); mrdivide(B,A)
    array([[1., 2., 3.]])
    """
    return np.linalg.lstsq(A.conj().T, B.conj().T, rcond=-1)[0].conj().T


def interpolate_freq_mag(
    freq: np.ndarray, mag: np.ndarray, N: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Interpolates freq and mag to N points but keeps the points intact.
    See Jason Orendorf's answer here: https://stackoverflow.com/questions/11884195/python-keep-points-in-spline-interpolation
    """

    assert freq.shape == mag.shape, "`freq` and `mag` must have the same shape"

    interpolated_freq = np.linspace(freq.min(), freq.max(), N - len(freq) + 2)
    new_freq = np.sort(
        np.append(interpolated_freq, freq[1:-1])
    )  # include the original points
    new_mag = scipy.interpolate.interp1d(freq, mag)(new_freq)
    new_freq = np.concatenate((new_freq, 2 - new_freq[-2:0:-1]))
    new_mag = np.concatenate((new_mag, new_mag[-2:0:-1]))

    return new_freq, new_mag


def polystab(p: np.ndarray) -> np.ndarray:
    """
    Reflects the roots of polynomial `p` that are outside the unit circle back inside.

    >>> polystab(np.array([[1.,2.,1.]]))
    array([1., 2., 1.])

    >>> polystab(np.array([[1.,2.,1.01]]))
    array([1.        , 1.98019802, 0.99009901])
    """
    if p.ndim == 2 and (p.shape[0] == 1 or p.shape[1] == 1):
        p = p.flatten()
    if p.size == 1:
        return p

    v = np.roots(p)
    i = np.where(v != 0)[0]
    vs = 0.5 * (np.sign(np.abs(v[i]) - 1) + 1)
    v[i] = (1 - vs) * v[i] + vs / np.conj(v[i])
    b = np.poly(v)
    if not np.imag(p).any():
        b = np.real(b)
    return b


def numf(
    impulse_response: np.ndarray, denominator: np.ndarray, numerator_order: int
) -> np.ndarray:
    """
    Find the numberator of the impulse response `impulse_response` of a frequency
    response H = numerator/denominator.
    """
    nh = np.max(impulse_response.shape)

    impr = scipy.signal.lfilter(
        [1], denominator.flatten(), scipy.signal.unit_impulse(nh)
    )

    toep = scipy.linalg.toeplitz(impr, scipy.signal.unit_impulse(numerator_order + 1))

    b = mrdivide(impulse_response, toep.conj().T)
    return b


def yulewalk(
    filter_order: int, frequencies: np.ndarray, magnitudes: np.ndarray, npt: int = 512
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Recursive filter-design of an arbitrary frequency response using the
    modified Yule-Walker method [1].

    :param filter_order: Estimated filter order.
    :param frequencies: Array of frequencies where the points of magnitude are positioned.
    :param magnitudes: The (linear) magnitudes at points defined by `frequencies`.
    :param N: Number points to use when estimating the frequency response.

    [1] Friedlander, Benjamin, and Boaz Porat. "The modified Yule-Walker method of ARMA
        spectral estimation."
        IEEE Transactions on Aerospace and Electronic Systems 2 (1984): 158-173.
    """

    if frequencies.ndim == 1:
        frequencies = frequencies[:, np.newaxis]
    if magnitudes.ndim == 1:
        magnitudes = magnitudes[:, np.newaxis]

    thelap = np.fix(npt / 25)

    assert len(frequencies) == len(magnitudes)

    nbrk = len(frequencies)

    npt = npt + 1
    Ht = np.zeros((1, npt))

    nint = nbrk - 1
    df = frequencies[1:] - frequencies[:-1]

    assert (df >= 0).all(), "f must be in increasing order"

    nb = 0
    Ht[0, 0] = magnitudes[0]
    for i in range(nint):
        if df[i] == 0:
            nb = nb - thelap / 2
            ne = nb + thelap
        else:
            ne = int(frequencies[i + 1] * npt)

        assert (
            nb >= 0 and ne <= npt
        ), f"nb={nb}, ne={ne} - Too abrupt change near end of frequency range"

        j = np.arange(nb, ne).astype(int)
        if ne == nb:
            inc = 0
        else:
            inc = (j - nb) / (ne - nb)

        Ht[0, j] = inc * magnitudes[i + 1] + (1 - inc) * magnitudes[i]
        nb = ne

    Ht[0, :] = np.interp(
        np.linspace(0, 1, npt), frequencies.flatten(), magnitudes.flatten()
    )

    # Symmetrize the frequency response
    Ht = np.hstack(
        [
            Ht,
            Ht[:, -2:0:-1],
        ]
    )

    n = Ht.shape[1]
    n2 = int((n + 1) / 2)
    nb = filter_order
    nr = 4 * filter_order
    nt = np.arange(0, nr).astype(int)
    R = np.real(np.fft.ifft(Ht * Ht))
    R = R[:, :nr] * (0.54 + 0.46 * np.cos(np.pi * nt / (nr - 1)))
    Rwindow = np.hstack(
        [np.array([[0.5]]), np.ones((1, n2 - 1)), np.zeros((1, n - n2))]
    )
    nr = R.shape[1]
    Rm = scipy.linalg.toeplitz(R[:, filter_order : nr - 1], R[:, filter_order:0:-1])
    Rhs = -R[:, filter_order + 1 : nr]
    denf = np.hstack(
        [
            np.array([[1.0]]),
            np.linalg.lstsq(Rm, Rhs.T, rcond=-1)[0].T,
        ]  # Rhs/Rm' in matlab
    )
    A = polystab(denf)

    hh = np.copy(R)
    hh[0, 0] = R[0, 0] / 2

    Qh = numf(hh, A, filter_order)
    aa = A.reshape(-1, 1)
    bb = Qh.reshape(-1, 1)

    nna = np.max(aa.shape)
    nnb = np.max(bb.shape)

    ssb = np.zeros((1, n))
    ssb[0, :nnb] = Qh

    ssa = np.zeros((1, n))
    ssa[0, :nna] = A

    Ss = 2 * np.real((np.fft.fft(ssb) / np.fft.fft(ssa)))
    hh = np.fft.ifft(np.exp(np.fft.fft(Rwindow * np.fft.ifft(complex_log(Ss)))))

    impr = scipy.signal.lfilter([1], A.flatten(), scipy.signal.unit_impulse(nr))
    B = np.real(
        mrdivide(
            hh[0, :nr],
            scipy.linalg.toeplitz(impr, scipy.signal.unit_impulse(nb + 1)).conj().T,
        )
    )

    B = np.real(numf(hh[0, :nr], A, nb))
    return A, B


if __name__ == "__main__":
    import doctest

    doctest.testmod()
