"""The waveform functions of a separately built liblalsimulation through ctypes.

``scripts/generate_reference_data.py --lalsimulation-library PATH`` uses this
to generate reference data with a lalsimulation that is not the one of the
installed SWIG module, for example lalsimulation 6.2.0 with a patch from
lalsuite master built by ``scripts/build_patched_lalsimulation.sh``. (The
SWIG module ``_lalsimulation.so`` contains its own copy of the library, so
the library cannot be swapped underneath it.)

:class:`Library` provides the functions the reference generator calls, with
the arguments and return values of the SWIG wrappers. LALDicts and mode
arrays are still made with the SWIG modules: they live in liblal, which the
library loaded here shares with them (it must be linked against the same
liblal), and are passed as pointers.
"""
import ctypes
from types import SimpleNamespace

import numpy as np

c_double = ctypes.c_double
c_int = ctypes.c_int
c_uint = ctypes.c_uint
c_void_p = ctypes.c_void_p


class LALUnit(ctypes.Structure):
    _fields_ = [
        ("powerOfTen", ctypes.c_int16),
        ("unitNumerator", ctypes.c_int16 * 7),
        ("unitDenominatorMinusOne", ctypes.c_uint16 * 7),
    ]


class COMPLEX16Sequence(ctypes.Structure):
    _fields_ = [("length", c_uint), ("data", ctypes.POINTER(c_double))]


class REAL8Sequence(ctypes.Structure):
    _fields_ = [("length", c_uint), ("data", ctypes.POINTER(c_double))]


class COMPLEX16FrequencySeries(ctypes.Structure):
    _fields_ = [
        ("name", ctypes.c_char * 64),
        ("gpsSeconds", ctypes.c_int32),
        ("gpsNanoSeconds", ctypes.c_int32),
        ("f0", c_double),
        ("deltaF", c_double),
        ("sampleUnits", LALUnit),
        ("data", ctypes.POINTER(COMPLEX16Sequence)),
    ]


class LALVCSInfo(ctypes.Structure):
    _fields_ = [
        (name, ctypes.c_char_p)
        for name in (
            "name", "version", "vcsId", "vcsDate", "vcsBranch", "vcsTag", "vcsAuthor",
            "vcsCommitter", "vcsClean", "vcsStatus", "configureArgs", "configureDate", "buildDate",
        )
    ]


class SphHarmFrequencySeries(ctypes.Structure):
    pass


SphHarmFrequencySeries._fields_ = [
    ("mode", ctypes.POINTER(COMPLEX16FrequencySeries)),
    ("l", c_uint),
    ("m", c_int),
    ("fdata", ctypes.POINTER(REAL8Sequence)),
    ("next", ctypes.POINTER(SphHarmFrequencySeries)),
]

SeriesPointer = ctypes.POINTER(COMPLEX16FrequencySeries)


def _pointer(swig_object):
    """The address of a SWIG-wrapped LAL object (``None`` for ``NULL``)."""
    return None if swig_object is None else c_void_p(int(swig_object.this))


def _series(pointer):
    """A COMPLEX16FrequencySeries as the SWIG wrapper presents it
    (``series.data.data``, ``series.data.length``)."""
    sequence = pointer.contents.data.contents
    n = sequence.length
    values = np.ctypeslib.as_array(sequence.data, shape=(2 * n,)).copy()
    return SimpleNamespace(data=SimpleNamespace(data=values[0::2] + 1j * values[1::2], length=n))


class Library:
    """The waveform functions of the liblalsimulation at ``path``."""

    def __init__(self, path):
        self.path = str(path)
        lib = ctypes.CDLL(self.path)
        self.lib = lib
        series = ctypes.POINTER(SeriesPointer)
        sequence = ctypes.POINTER(REAL8Sequence)
        lib.XLALSimInspiralChooseFDWaveformSequence.argtypes = (
            [series] * 2 + [c_double] * 12 + [c_void_p, c_int, sequence]
        )
        lib.XLALSimInspiralChooseFDWaveformSequence.restype = c_int
        lib.XLALSimInspiralChooseFDWaveform.argtypes = [series] * 2 + [c_double] * 18 + [c_void_p, c_int]
        lib.XLALSimInspiralChooseFDWaveform.restype = c_int
        lib.XLALSimInspiralChooseFDModes.argtypes = [c_double] * 15 + [c_void_p, c_int]
        lib.XLALSimInspiralChooseFDModes.restype = ctypes.POINTER(SphHarmFrequencySeries)
        lib.XLALSimIMRPhenomXPHMFrequencySequenceOneMode.argtypes = (
            [series] * 2 + [sequence, c_uint, c_int] + [c_double] * 12 + [c_void_p]
        )
        lib.XLALSimIMRPhenomXPHMFrequencySequenceOneMode.restype = c_int
        lib.XLALDestroyCOMPLEX16FrequencySeries.argtypes = [SeriesPointer]
        lib.XLALDestroyCOMPLEX16FrequencySeries.restype = None
        lib.XLALDestroySphHarmFrequencySeries.argtypes = [ctypes.POINTER(SphHarmFrequencySeries)]
        lib.XLALDestroySphHarmFrequencySeries.restype = None
        lib.XLALCreateREAL8Sequence.argtypes = [c_uint]
        lib.XLALCreateREAL8Sequence.restype = sequence
        lib.XLALDestroyREAL8Sequence.argtypes = [sequence]
        lib.XLALDestroyREAL8Sequence.restype = None
        lib.XLALGetBaseErrno.restype = c_int
        lib.XLALClearErrno.restype = None
        info = LALVCSInfo.in_dll(lib, "lalSimulationVCSInfo")
        self.vcs_info = {name: (getattr(info, name) or b"").decode() for name, _ in LALVCSInfo._fields_}
        self.__version__ = self.vcs_info["version"]

    def _sequence(self, values):
        values = np.asarray(values, dtype=float)
        sequence = self.lib.XLALCreateREAL8Sequence(len(values))
        np.ctypeslib.as_array(sequence.contents.data, shape=(len(values),))[:] = values
        return sequence

    def _check(self, status, name):
        if status != 0:
            errno = self.lib.XLALGetBaseErrno()
            self.lib.XLALClearErrno()
            raise RuntimeError(f"{name} failed (XLAL error {errno})")

    def _pair(self, first, second):
        out = _series(first), _series(second)
        self.lib.XLALDestroyCOMPLEX16FrequencySeries(first)
        self.lib.XLALDestroyCOMPLEX16FrequencySeries(second)
        return out

    def SimInspiralChooseFDWaveformSequence(
        self, phi_ref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f_ref, distance, inclination, params,
        approximant, frequencies,
    ):
        hp, hc = SeriesPointer(), SeriesPointer()
        sequence = self._sequence(frequencies.data)
        status = self.lib.XLALSimInspiralChooseFDWaveformSequence(
            ctypes.byref(hp), ctypes.byref(hc), phi_ref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f_ref,
            distance, inclination, _pointer(params), int(approximant), sequence,
        )
        self.lib.XLALDestroyREAL8Sequence(sequence)
        self._check(status, "XLALSimInspiralChooseFDWaveformSequence")
        return self._pair(hp, hc)

    def SimInspiralChooseFDWaveform(
        self, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination, phi_ref, long_asc_nodes,
        eccentricity, mean_per_ano, delta_f, f_min, f_max, f_ref, params, approximant,
    ):
        hp, hc = SeriesPointer(), SeriesPointer()
        status = self.lib.XLALSimInspiralChooseFDWaveform(
            ctypes.byref(hp), ctypes.byref(hc), m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance,
            inclination, phi_ref, long_asc_nodes, eccentricity, mean_per_ano, delta_f, f_min, f_max,
            f_ref, _pointer(params), int(approximant),
        )
        self._check(status, "XLALSimInspiralChooseFDWaveform")
        return self._pair(hp, hc)

    def SimInspiralChooseFDModes(
        self, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, delta_f, f_min, f_max, f_ref, phi_ref, distance,
        inclination, params, approximant,
    ):
        hlms = self.lib.XLALSimInspiralChooseFDModes(
            m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, delta_f, f_min, f_max, f_ref, phi_ref, distance,
            inclination, _pointer(params), int(approximant),
        )
        if not hlms:
            self._check(1, "XLALSimInspiralChooseFDModes")
        nodes = list()
        node = hlms
        while node:
            contents = node.contents
            nodes.append(SimpleNamespace(mode=_series(contents.mode), l=contents.l, m=contents.m, next=None))
            node = contents.next
        for current, following in zip(nodes, nodes[1:]):
            current.next = following
        self.lib.XLALDestroySphHarmFrequencySeries(hlms)
        return nodes[0] if nodes else None

    def SimIMRPhenomXPHMFrequencySequenceOneMode(
        self, frequencies, ell, emm, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination,
        phi_ref, f_ref, params,
    ):
        positive, negative = SeriesPointer(), SeriesPointer()
        sequence = self._sequence(frequencies.data)
        status = self.lib.XLALSimIMRPhenomXPHMFrequencySequenceOneMode(
            ctypes.byref(positive), ctypes.byref(negative), sequence, ell, emm, m1, m2, s1x, s1y, s1z,
            s2x, s2y, s2z, distance, inclination, phi_ref, f_ref, _pointer(params),
        )
        self.lib.XLALDestroyREAL8Sequence(sequence)
        self._check(status, "XLALSimIMRPhenomXPHMFrequencySequenceOneMode")
        return self._pair(positive, negative)


class Overlay:
    """The SWIG module ``swig`` with the waveform functions of ``library``."""

    def __init__(self, swig, library):
        self._swig = swig
        self._library = library
        self.__version__ = library.__version__

    def __getattr__(self, name):
        if name.startswith("Sim") and hasattr(self._library, name):
            return getattr(self._library, name)
        return getattr(self._swig, name)
