import ctypes
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class ThreeBodyPhaseSpace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('gfortran')
        if compiler is None:
            raise unittest.SkipTest('gfortran required')
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        folder = Path(cls.temp.name)
        source = (Path(__file__).resolve().parents[1] / 'physics_Xphasespace.f').read_text()
        start = source.index('      real*8 function xphase_threebody(')
        (folder / 'phase.f').write_text(source[start:])
        library = folder / 'phase.so'
        subprocess.run([compiler, '-shared', '-fPIC', '-O2',
                        '-ffixed-line-length-132', str(folder/'phase.f'),
                        '-o', str(library)], check=True)
        cls.lib = ctypes.CDLL(str(library))
        cls.function = cls.lib.xphase_threebody_
        cls.function.restype = ctypes.c_double
        cls.function.argtypes = [ctypes.POINTER(ctypes.c_double)] * 4

    def phase(self, mass, m1, m2, m3):
        args = [ctypes.c_double(x) for x in (mass, m1, m2, m3)]
        return self.function(*(ctypes.byref(x) for x in args))

    def test_threshold(self):
        for mass in (0, 300, 414.11696):
            self.assertEqual(self.phase(mass, 139.57018, 139.57018, 134.9766), 0)
        self.assertGreater(self.phase(415, 139.57018, 139.57018, 134.9766), 0)

    def test_massless_limit(self):
        for mass in (500, 1000, 2000):
            expected = mass**2/(256*math.pi**3)
            self.assertAlmostEqual(self.phase(mass, 0, 0, 0)/expected, 1, places=6)

    def test_permutation_and_scaling(self):
        masses = (139.57018, 139.57018, 134.9766)
        for mass in (420, 600, 1000, 1800):
            value = self.phase(mass, *masses)
            self.assertTrue(math.isfinite(value))
            self.assertGreater(value, 0)
            swapped = self.phase(mass, masses[2], masses[0], masses[1])
            self.assertAlmostEqual(swapped/value, 1, places=6)
            doubled = self.phase(2*mass, *(2*x for x in masses))
            self.assertAlmostEqual(doubled/value, 4, places=10)


if __name__ == '__main__':
    unittest.main()
