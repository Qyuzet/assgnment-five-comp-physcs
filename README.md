# Computational Physics – Week 5 Assignment
**Course:** SCIE6063001 – Computational Physics
**Textbook:** Newman, M. (2012). *Computational Physics*. CreateSpace. ISBN: 1480145513

---

## Problem

Evaluate the integral:

$$I = \int_0^1 \sin^2\!\left(\sqrt{100x}\right)\, dx \approx 0.45583253$$

Analytical result (via substitution $u = \sqrt{100x}$ and integration by parts):

$$I = \frac{1}{100}\left[50 - 5\sin(20) - \frac{\cos(20)}{4} + \frac{1}{4}\right]$$

---

## Exercise 5.7 — [`exercise_5_7.py`](exercise_5_7.py)

### (a) Adaptive Trapezoidal Rule *(Eq. 5.34)*

Starts with $N = 1$ slice and doubles until the error estimate falls below $\varepsilon = 10^{-6}$.

Recurrence relation:
$$I_{2N} = \frac{1}{2}I_N + \frac{h}{2}\sum_{k=0}^{N-1} f\!\left(a + \left(k+\tfrac{1}{2}\right)h\right)$$

Error estimate: $\varepsilon \approx \tfrac{1}{3}|I_{2N} - I_N|$

**Result:** converged at **N = 4096 slices**

### (b) Romberg Integration *(Eq. 5.49)*

Builds a triangular table by applying Richardson extrapolation to successive trapezoidal estimates:

$$R_{i,j} = R_{i,j-1} + \frac{R_{i,j-1} - R_{i-1,j-1}}{4^j - 1}$$

**Result:** converged at **level i = 7 (N = 128 slices)**

---

## Exercise 5.8 — [`exercise_5_8.py`](exercise_5_8.py)

### Adaptive Simpson's Rule *(Eqs. 5.35–5.39)*

Derives Simpson's estimates from pairs of trapezoidal estimates:

$$S_{2N} = \frac{4T_{2N} - T_N}{3}$$

Error estimate: $\varepsilon \approx \tfrac{1}{15}|S_{2N} - S_N|$
(factor 15 = $2^4 - 1$, since Simpson's rule is $O(h^4)$)

**Result:** converged at **N = 256 slices**

---

## Comparison

| Method               | Slices to reach $\varepsilon = 10^{-6}$ |
|----------------------|----------------------------------------:|
| Adaptive Trapezoidal | 4096                                    |
| Adaptive Simpson's   | 256                                     |
| Romberg              | 128                                     |

Romberg integration reaches the target accuracy considerably faster than the trapezoidal rule alone, while adaptive Simpson's rule lies between the two — consistent with Newman's commentary on p. 163.

---

## Usage

```bash
python exercise_5_7.py
python exercise_5_8.py
```

Requires Python 3 (standard library only — no external dependencies).
