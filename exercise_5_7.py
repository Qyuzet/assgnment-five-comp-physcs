# exercise_5_7.py
# Newman, Computational Physics (2012), Exercise 5.7 (pp. 162-163)
#
# Evaluate  I = integral_0^1 sin^2(sqrt(100*x)) dx
#
# (a) Adaptive Trapezoidal Rule  (Section 5.3, Eq. 5.34)
# (b) Romberg Integration        (Section 5.3, Eq. 5.49)
#
# Target accuracy: eps = 1e-6
#
# Analytical result (substitution u = sqrt(100x), then integration by parts):
#   I = (1/100) * [50 - 5*sin(20) - cos(20)/4 + 1/4]  ~= 0.45583

from math import sin, sqrt, cos


def f(x):
    """Integrand: sin^2(sqrt(100 x))."""
    return sin(sqrt(100 * x)) ** 2


a, b = 0.0, 1.0
eps = 1e-6

# Analytical reference value
I_exact = (1 / 100) * (50 - 5 * sin(20) - cos(20) / 4 + 0.25)

# ─── Part (a): Adaptive Trapezoidal Rule ─────────────────────────────────────

print("=" * 66)
print("Exercise 5.7(a): Adaptive Trapezoidal Rule")
print("=" * 66)
print(f"Integrand : f(x) = sin^2(sqrt(100x))")
print(f"Interval  : [{a}, {b}]")
print(f"Tolerance : eps = {eps}")
print(f"Exact     : I = {I_exact:.10f}\n")

# --- Initialise with N = 1 slice ---
N = 1
h = b - a
I = 0.5 * h * (f(a) + f(b))

print(f"{'N':>10}   {'Integral':>15}   {'Error est.':>12}")
print(f"{'-'*10}   {'-'*15}   {'-'*12}")
print(f"{N:>10}   {I:>15.10f}   {'---':>12}")

while True:
    # Eq. (5.34): new estimate = 0.5 * old + (h/2) * sum of f at new midpoints
    new_sum = sum(f(a + (k + 0.5) * h) for k in range(N))
    I_new = 0.5 * I + 0.5 * h * new_sum

    # Error estimate: (1/3) |I_new - I|   [Richardson, trapezoidal is O(h^2)]
    error = abs(I_new - I) / 3

    N *= 2
    h /= 2
    I = I_new

    print(f"{N:>10}   {I:>15.10f}   {error:>12.3e}")

    if error < eps:
        break

print(f"\nResult  : I = {I:.10f}")
print(f"Exact   : I = {I_exact:.10f}")
print(f"Slices  : N = {N}")

# ─── Part (b): Romberg Integration ───────────────────────────────────────────

print("\n" + "=" * 66)
print("Exercise 5.7(b): Romberg Integration")
print("=" * 66)
print(f"Integrand : f(x) = sin^2(sqrt(100x))")
print(f"Interval  : [{a}, {b}]")
print(f"Tolerance : eps = {eps}")
print(f"Exact     : I = {I_exact:.10f}\n")

# R[i][j]:
#   j = 0  →  trapezoidal with N = 2^i slices
#   j >= 1 →  Richardson extrapolation:
#             R[i][j] = R[i][j-1] + (R[i][j-1] - R[i-1][j-1]) / (4^j - 1)
R = []

# Level 0: trapezoidal with N = 1
h = b - a
R.append([0.5 * h * (f(a) + f(b))])

print("Romberg table  R[i][j]  (rows = level i, columns = extrapolation order j):")
print(f"  {'  '.join(f'{v:.8f}' for v in R[0])}")

for i in range(1, 30):
    h /= 2
    n_mid = 2 ** (i - 1)  # number of new midpoints at level i

    # New trapezoidal estimate: reuse previous via Eq. (5.34)
    new_sum = sum(f(a + (2 * k + 1) * h) for k in range(n_mid))
    Ri0 = 0.5 * R[i - 1][0] + h * new_sum

    row = [Ri0]
    for j in range(1, i + 1):
        # Eq. (5.49)
        Rij = row[j - 1] + (row[j - 1] - R[i - 1][j - 1]) / (4 ** j - 1)
        row.append(Rij)

    R.append(row)
    print(f"  {'  '.join(f'{v:.8f}' for v in row)}")

    # Convergence check: compare best estimates of consecutive levels
    error = abs(R[i][i] - R[i - 1][i - 1])
    if error < eps:
        print(f"\nResult : I = {R[i][i]:.10f}")
        print(f"Exact  : I = {I_exact:.10f}")
        print(f"Level  : i = {i}  (N = {2**i} slices,  error ~ {error:.3e})")
        break
