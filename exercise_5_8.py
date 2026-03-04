# exercise_5_8.py
# Newman, Computational Physics (2012), Exercise 5.8 (p. 163)
#
# Evaluate  I = integral_0^1 sin^2(sqrt(100*x)) dx
#
# Method: Adaptive Simpson's Rule  (Section 5.3, Eqs. 5.35–5.39)
#
# The key recurrence (Eq. 5.36):
#   S_{2N} = (4 * T_{2N} - T_N) / 3
# where T_N is the trapezoidal estimate with N slices.
# Halving h from T_N → T_{2N} is done with Eq. (5.34).
#
# Error estimate (Richardson, Simpson's is O(h^4)):
#   eps_i ~ |S_{new} - S_{old}| / (2^4 - 1) = |S_{new} - S_{old}| / 15
#
# Start with N_S = 2 slices; double until error < 1e-6.
#
# Analytical result:
#   I = (1/100) * [50 - 5*sin(20) - cos(20)/4 + 1/4]  ~= 0.45583

from math import sin, sqrt, cos


def f(x):
    """Integrand: sin^2(sqrt(100 x))."""
    return sin(sqrt(100 * x)) ** 2


a, b = 0.0, 1.0
eps = 1e-6

# Analytical reference value
I_exact = (1 / 100) * (50 - 5 * sin(20) - cos(20) / 4 + 0.25)

print("=" * 66)
print("Exercise 5.8: Adaptive Simpson's Rule")
print("=" * 66)
print(f"Integrand : f(x) = sin^2(sqrt(100x))")
print(f"Interval  : [{a}, {b}]")
print(f"Tolerance : eps = {eps}")
print(f"Exact     : I = {I_exact:.10f}\n")

# ─── Seed two trapezoidal estimates ──────────────────────────────────────────
# T1: N = 1 slice, step h = 1
# T2: N = 2 slices, step h/2 = 0.5
h = b - a
T1 = 0.5 * h * (f(a) + f(b))
T2 = 0.5 * T1 + 0.5 * h * f(a + 0.5 * h)

# Initial Simpson's estimate: S2 = (4*T2 - T1)/3  (N_S = 2)
S = (4 * T2 - T1) / 3
T_old = T2          # the "older" trapezoidal estimate for the next iteration
N = 2               # current number of Simpson's slices
h = h / 2           # step corresponding to T_old (= (b-a)/N)

print(f"{'N_S':>10}   {'Integral':>18}   {'Error est.':>12}")
print(f"{'-'*10}   {'-'*18}   {'-'*12}")
print(f"{N:>10}   {S:>18.10f}   {'---':>12}")

while True:
    # Double Simpson's slices: N → 2N
    # Compute T_{2N} from T_old (= T_N) via Eq. (5.34):
    #   T_{2N} = 0.5 * T_N + (h/2) * sum_k f(a + (k+0.5)*h)
    new_sum = sum(f(a + (k + 0.5) * h) for k in range(N))
    T_new = 0.5 * T_old + 0.5 * h * new_sum

    # New Simpson's estimate: S_{2N} = (4*T_{2N} - T_N) / 3
    S_new = (4 * T_new - T_old) / 3

    # Error estimate: |S_{2N} - S_N| / 15   (Simpson's is O(h^4) → factor 15)
    error = abs(S_new - S) / 15

    N *= 2
    h /= 2
    T_old = T_new
    S = S_new

    print(f"{N:>10}   {S:>18.10f}   {error:>12.3e}")

    if error < eps:
        break

print(f"\nResult  : I = {S:.10f}")
print(f"Exact   : I = {I_exact:.10f}")
print(f"Slices  : N = {N}")
