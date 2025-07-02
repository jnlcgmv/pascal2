import unittest
from sun200 import sun_position


class TestSun200(unittest.TestCase):
    def test_sample(self):
        T = 0.244274469541410
        L, B, R, s1, s2, s3 = sun_position(T)
        self.assertAlmostEqual(L, 75.4382471373045, places=9)
        self.assertAlmostEqual(B, 0.000083432781694, places=12)
        self.assertAlmostEqual(R, 1.014712190499918, places=12)
        self.assertAlmostEqual(s1, 38165753.685285926, places=6)
        self.assertAlmostEqual(s2, 134798747.04130605, places=6)
        self.assertAlmostEqual(s3, 58442650.25013706, places=6)


if __name__ == "__main__":
    unittest.main()
