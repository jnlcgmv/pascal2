import symforce.symbolic as sf
import math

PI = math.pi

# Conversion of the SUN200 Pascal routine to Python using symforce compatible
# operations. The algorithm computes the ecliptic coordinates (L, B, R) of the
# Sun referred to the mean equinox of date.

def frac(x: float) -> float:
    """Return the fractional part of x."""
    x = x - math.trunc(x)
    return x + 1.0 if x < 0 else x


def add_the(c1: float, s1: float, c2: float, s2: float) -> tuple[float, float]:
    """Return cosine and sine of the sum of two angles given by their cos/sin."""
    c = c1 * c2 - s1 * s2
    s = s1 * c2 + c1 * s2
    return c, s


def sun200(T: float) -> tuple[float, float, float]:
    """Compute L, B, R for a given T using the SUN200 algorithm."""
    C3 = {i: 0.0 for i in range(-1, 8)}
    S3 = {i: 0.0 for i in range(-1, 8)}
    C = {i: 0.0 for i in range(-8, 1)}
    S = {i: 0.0 for i in range(-8, 1)}

    DL = 0.0
    DR = 0.0
    DB = 0.0
    U = 0.0
    V = 0.0

    def term(i1: int, i: int, it: int, dlc: float, dls: float,
             drc: float, drs: float, dbc: float, dbs: float) -> None:
        nonlocal U, V, DL, DR, DB
        if it == 0:
            U, V = add_the(C3[i1], S3[i1], C[i], S[i])
        else:
            U *= T
            V *= T
        DL += dlc * U + dls * V
        DR += drc * U + drs * V
        DB += dbc * U + dbs * V

    def pertven() -> None:
        C[0] = 1.0
        S[0] = 0.0
        C[-1] = math.cos(M2)
        S[-1] = -math.sin(M2)
        for ii in range(-1, -5 - 1, -1):
            C[ii - 1], S[ii - 1] = add_the(C[ii], S[ii], C[-1], S[-1])
        term(1, 0, 0, -0.22, 6892.76, -16707.37, -0.54, 0.0, 0.0)
        term(1, 0, 1, -0.06, -17.35, 42.04, -0.15, 0.0, 0.0)
        term(1, 0, 2, -0.01, -0.05, 0.13, -0.02, 0.0, 0.0)
        term(2, 0, 0, 0.0, 71.98, -139.57, 0.0, 0.0, 0.0)
        term(2, 0, 1, 0.0, -0.36, 0.70, 0.0, 0.0, 0.0)
        term(3, 0, 0, 0.0, 1.04, -1.75, 0.0, 0.0, 0.0)
        term(0, -1, 0, 0.03, -0.07, -0.16, -0.07, 0.02, -0.02)
        term(1, -1, 0, 2.35, -4.23, -4.75, -2.64, 0.0, 0.0)
        term(1, -2, 0, -0.10, 0.06, 0.12, 0.20, 0.02, 0.0)
        term(2, -1, 0, -0.06, -0.03, 0.20, -0.01, 0.01, -0.09)
        term(2, -2, 0, -4.70, 2.90, 8.28, 13.42, 0.01, -0.01)
        term(3, -2, 0, 1.80, -1.74, -1.44, -1.57, 0.04, -0.06)
        term(3, -3, 0, -0.67, 0.03, 0.11, 2.43, 0.01, 0.0)
        term(4, -2, 0, 0.03, -0.03, 0.10, 0.09, 0.01, -0.01)
        term(4, -3, 0, 1.51, -0.40, -0.88, -3.36, 0.18, -0.10)
        term(4, -4, 0, -0.19, -0.09, -0.38, 0.77, 0.0, 0.0)
        term(5, -3, 0, 0.76, -0.68, 0.30, 0.37, 0.01, 0.0)
        term(5, -4, 0, -0.14, -0.04, -0.11, 0.43, -0.03, 0.0)
        term(5, -5, 0, -0.05, -0.07, -0.31, 0.21, 0.0, 0.0)
        term(6, -4, 0, 0.15, -0.04, -0.06, -0.21, 0.01, 0.0)
        term(6, -5, 0, -0.03, -0.03, -0.09, 0.09, -0.01, 0.0)
        term(6, -6, 0, 0.00, -0.04, -0.18, 0.02, 0.0, 0.0)
        term(7, -5, 0, -0.12, -0.03, -0.08, 0.31, -0.02, -0.01)

    def pertmar() -> None:
        C[-1] = math.cos(M4)
        S[-1] = -math.sin(M4)
        for ii in range(-1, -7 - 1, -1):
            C[ii - 1], S[ii - 1] = add_the(C[ii], S[ii], C[-1], S[-1])
        term(1, -1, 0, -0.22, 0.17, -0.21, -0.27, 0.0, 0.0)
        term(1, -2, 0, -1.66, 0.62, 0.16, 0.28, 0.0, 0.0)
        term(2, -2, 0, 1.96, 0.57, -1.32, 4.55, 0.0, 0.01)
        term(2, -3, 0, 0.40, 0.15, -0.17, 0.46, 0.0, 0.0)
        term(2, -4, 0, 0.53, 0.26, 0.09, -0.22, 0.0, 0.0)
        term(3, -3, 0, 0.05, 0.12, -0.35, 0.15, 0.0, 0.0)
        term(3, -4, 0, -0.13, -0.48, 1.06, -0.29, 0.01, 0.0)
        term(3, -5, 0, -0.04, -0.20, 0.20, -0.04, 0.0, 0.0)
        term(4, -4, 0, 0.00, -0.03, 0.10, 0.04, 0.0, 0.0)
        term(4, -5, 0, 0.05, -0.07, 0.20, 0.14, 0.0, 0.0)
        term(4, -6, 0, -0.10, 0.11, -0.23, -0.22, 0.0, 0.0)
        term(5, -7, 0, -0.05, 0.00, 0.01, -0.14, 0.0, 0.0)
        term(5, -8, 0, 0.05, 0.01, -0.02, 0.10, 0.0, 0.0)

    def pertjup() -> None:
        C[-1] = math.cos(M5)
        S[-1] = -math.sin(M5)
        for ii in range(-1, -3 - 1, -1):
            C[ii - 1], S[ii - 1] = add_the(C[ii], S[ii], C[-1], S[-1])
        term(-1, -1, 0, 0.01, 0.07, 0.18, -0.02, 0.0, -0.02)
        term(0, -1, 0, -0.31, 2.58, 0.52, 0.34, 0.02, 0.0)
        term(1, -1, 0, -7.21, -0.06, 0.13, -16.27, 0.0, -0.02)
        term(1, -2, 0, -0.54, -1.52, 3.09, -1.12, 0.01, -0.17)
        term(1, -3, 0, -0.03, -0.21, 0.38, -0.06, 0.0, -0.02)
        term(2, -1, 0, -0.16, 0.05, -0.18, -0.31, 0.01, 0.0)
        term(2, -2, 0, 0.14, -2.73, 9.23, 0.48, 0.0, 0.0)
        term(2, -3, 0, 0.07, -0.55, 1.83, 0.25, 0.01, 0.0)
        term(2, -4, 0, 0.02, -0.08, 0.25, 0.06, 0.0, 0.0)
        term(3, -2, 0, 0.01, -0.07, 0.16, 0.04, 0.0, 0.0)
        term(3, -3, 0, -0.16, -0.03, 0.08, -0.64, 0.0, 0.0)
        term(3, -4, 0, -0.04, -0.01, 0.03, -0.17, 0.0, 0.0)

    def pertsat() -> None:
        C[-1] = math.cos(M6)
        S[-1] = -math.sin(M6)
        C[-2], S[-2] = add_the(C[-1], S[-1], C[-1], S[-1])
        term(0, -1, 0, 0.00, 0.32, 0.01, 0.0, 0.0, 0.0)
        term(1, -1, 0, -0.08, -0.41, 0.97, -0.18, 0.0, -0.01)
        term(1, -2, 0, 0.04, 0.10, -0.23, 0.10, 0.0, 0.0)
        term(2, -2, 0, 0.04, 0.10, -0.35, 0.13, 0.0, 0.0)

    def pertmoo() -> None:
        nonlocal DL, DR, DB
        DL += 6.45 * math.sin(D) - 0.42 * math.sin(D - A) + 0.18 * math.sin(D + A) \
                + 0.17 * math.sin(D - M3) - 0.06 * math.sin(D + M3)
        DR += 30.76 * math.cos(D) - 3.06 * math.cos(D - A) + 0.85 * math.cos(D + A) \
                - 0.58 * math.cos(D + M3) + 0.57 * math.cos(D - M3)
        DB += 0.576 * math.sin(UU)

    # Main computation -----------------------------------------------------
    M2 = 2.0 * PI * frac(0.1387306 + 162.5485917 * T)
    M3 = 2.0 * PI * frac(0.9931266 + 99.9973604 * T)
    M4 = 2.0 * PI * frac(0.0543250 + 53.1666028 * T)
    M5 = 2.0 * PI * frac(0.0551750 + 8.4293972 * T)
    M6 = 2.0 * PI * frac(0.8816500 + 3.3938722 * T)
    D = 2.0 * PI * frac(0.8274 + 1236.8531 * T)
    A = 2.0 * PI * frac(0.3749 + 1325.5524 * T)
    UU = 2.0 * PI * frac(0.2591 + 1342.2278 * T)

    C3[0] = 1.0
    S3[0] = 0.0
    C3[1] = math.cos(M3)
    S3[1] = math.sin(M3)
    C3[-1] = C3[1]
    S3[-1] = -S3[1]
    for ii in range(2, 8):
        C3[ii], S3[ii] = add_the(C3[ii - 1], S3[ii - 1], C3[1], S3[1])

    pertven()
    pertmar()
    pertjup()
    pertsat()
    pertmoo()

    DL += 6.40 * math.sin(2.0 * PI * (0.6983 + 0.0561 * T)) + \
          1.87 * math.sin(2.0 * PI * (0.5764 + 0.4174 * T)) + \
          0.27 * math.sin(2.0 * PI * (0.4189 + 0.3306 * T)) + \
          0.20 * math.sin(2.0 * PI * (0.3581 + 2.4814 * T))

    L = 360.0 * frac(0.7859453 + M3 / (2.0 * PI) + ((6191.2 + 1.1 * T) * T + DL) / 1296.0e3)
    R = 1.0001398 - 0.0000007 * T + DR * 1e-6
    B = DB / 3600.0
    return L, B, R


def sun_position(T: float) -> tuple[float, float, float, float, float, float]:
    """Return L, B, R and Cartesian Sun position."""
    L, B, R = sun200(T)
    rad = 2 * PI / 360.0
    sun_pos_1 = R * 149597870.0 * (math.cos(L * rad) * math.cos(B * rad))
    sun_pos_2 = R * 149597870.0 * (
        math.cos(23.43929111 * rad) * math.sin(L * rad) * math.cos(B * rad)
        - math.sin(23.43929111 * rad) * math.sin(B * rad)
    )
    sun_pos_3 = R * 149597870.0 * (
        math.sin(23.43929111 * rad) * math.sin(L * rad) * math.cos(B * rad)
        + math.cos(23.43929111 * rad) * math.sin(B * rad)
    )
    return L, B, R, sun_pos_1, sun_pos_2, sun_pos_3


if __name__ == "__main__":
    T = 0.244274469541410
    L, B, R, s1, s2, s3 = sun_position(T)
    print(f"Para T={T:.15f}:")
    print(f"L (longitud ecl\xEDptica) = {L:.15f} grados")
    print(f"B (latitud ecl\xEDptica)   = {B:.15f} grados")
    print(f"R (distancia)           = {R:.15f} AU")
    print(f"sun_pos_1 (distancia)           = {s1:.15f} km")
    print(f"sun_pos_2 (distancia)           = {s2:.15f} km")
    print(f"sun_pos_3 (distancia)           = {s3:.15f} km")
