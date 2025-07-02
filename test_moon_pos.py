import unittest
from moon_pos import moon_position

class TestMoonPos(unittest.TestCase):
    def test_sample(self):
        T = 0.244274469541410
        L, B, R = moon_position(T)
        self.assertAlmostEqual(L, 63.73560436541438, places=9)
        self.assertAlmostEqual(B, 3.882806290546285, places=9)
        self.assertAlmostEqual(R, 58.51113474520343, places=9)

if __name__ == "__main__":
    unittest.main()
