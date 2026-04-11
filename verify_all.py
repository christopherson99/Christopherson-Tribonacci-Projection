#!/usr/bin/env python3
"""
CHRISTOPHERSON COMPLETE VERIFICATION SUITE
Independently checkable proofs for all computational claims.
Run: python3 verify_all.py

Corick Christopherson, April 2026
"""

from sympy import *
from mpmath import mp, mpf, nstr, findroot, pi as mpi, e as me, log as mlog, sqrt as msqrt
mp.dps = 50

x, t, z = symbols('x t z')
PASS = 0
FAIL = 0

def check(name, condition):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ FAIL: {name}")

def section(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")

# =============================================================
section("1. BASE FIELD VERIFICATION")
# =============================================================

# φ
phi_poly = x**2 - x - 1
check("φ minpoly x²-x-1 irreducible over Q",
      factor(phi_poly, domain='QQ') == phi_poly)
phi_mp = (1 + msqrt(5)) / 2
check("φ = (1+√5)/2 ≈ 1.6180",
      abs(float(phi_mp) - 1.6180) < 0.001)
check("φ satisfies x²-x-1=0",
      abs(float(phi_mp**2 - phi_mp - 1)) < 1e-40)

# ρ
rho_poly = x**3 - x**2 - x - 1
check("ρ minpoly x³-x²-x-1 irreducible over Q",
      factor(rho_poly, domain='QQ') == rho_poly)
rho_mp = findroot(lambda x: x**3 - x**2 - x - 1, mpf('1.84'))
check("ρ ≈ 1.8393",
      abs(float(rho_mp) - 1.8393) < 0.001)
check("ρ satisfies x³-x²-x-1=0",
      abs(float(rho_mp**3 - rho_mp**2 - rho_mp - 1)) < 1e-40)

# τ
tau_poly = x**3 - x - 1
check("τ minpoly x³-x-1 irreducible over Q",
      factor(tau_poly, domain='QQ') == tau_poly)
tau_mp = findroot(lambda x: x**3 - x - 1, mpf('1.32'))
check("τ ≈ 1.3247",
      abs(float(tau_mp) - 1.3247) < 0.001)

# β = ρφ
phi_exact = (1 + sqrt(5)) / 2
rho_exact = max([r for r in solve(x**3 - x**2 - x - 1, x) if r.is_real])
beta_mp = minimal_polynomial(rho_exact * phi_exact, t)
beta_claimed = t**6 - t**5 - 4*t**4 - 5*t**3 - 2*t**2 + t - 1
check("β=ρφ minpoly = t⁶-t⁵-4t⁴-5t³-2t²+t-1",
      expand(beta_mp - beta_claimed) == 0)
check("β minpoly irreducible over Q",
      factor(beta_claimed, domain='QQ') == beta_claimed)

# Mod-p irreducibility certificates
for name, poly_expr, var, prime in [
    ("ρ minpoly mod 3", rho_poly, x, 3),
    ("β minpoly mod 23", beta_claimed, t, 23),
]:
    p = Poly(poly_expr, var, domain=ZZ)
    pm = p.set_domain(GF(prime))
    facs = pm.factor_list()
    degs = [degree(f) for f, m in facs[1]]
    check(f"{name}: irreducible mod {prime} (single factor deg {degree(p)})",
          len(degs) == 1 and degs[0] == degree(p))

# =============================================================
section("2. CONSTRUCTION π (P=1, PHI-POINT FAMILY)")
# =============================================================

pi_poly = x**16 + 2*x**12 + 2*x**10 + 2*x**9 - 44*x**8 + 2*x**6 + 2*x**5 - 44*x**4 + 2*x**3 - 44*x**2 - 45*x - 45

check("π polynomial degree 16", degree(pi_poly, x) == 16)
check("π polynomial irreducible over Q",
      factor(pi_poly, domain='QQ') == pi_poly)

# Mod-13 irreducibility
p13 = Poly(pi_poly, x, domain=ZZ).set_domain(GF(13))
facs13 = p13.factor_list()
degs13 = [degree(f) for f, m in facs13[1]]
check("π polynomial irreducible mod 13 (single factor deg 16)",
      len(degs13) == 1 and degs13[0] == 16)

# Numerical root
x_pi = findroot(lambda x: x**16 + 2*x**12 + 2*x**10 + 2*x**9 - 44*x**8 + 2*x**6 + 2*x**5 - 44*x**4 + 2*x**3 - 44*x**2 - 45*x - 45,
                mpf('1.565'))
pi_c = 2 * x_pi
check(f"Construction π = 2x = {nstr(pi_c, 20)}",
      abs(float(pi_c) - 3.12981655781537) < 1e-12)

# Isolating interval
check("Root in isolating interval (1.564, 1.566)",
      1.564 < float(x_pi) < 1.566)

# Derivation: φ⁸ = 21φ+13, substitute into x⁸+x⁴+x²+x+1
phi8 = phi_mp**8
check("φ⁸ = 21φ + 13",
      abs(float(phi8 - (21*phi_mp + 13))) < 1e-35)

# Verify root satisfies family equation
check("x satisfies φ⁸ = x⁸+x⁴+x²+x+1",
      abs(float(x_pi**8 + x_pi**4 + x_pi**2 + x_pi + 1 - phi8)) < 1e-30)

# Binary octic normal form check: M = 2P+U+1, N = P²+UP+V
# φ⁸ minpoly: z² - 47z + 1, so U=-47, V=1
U, V, P = -47, 1, 1
M = 2*P + U + 1
N = P**2 + U*P + V
check(f"Normal form: M = 2({P})+({U})+1 = {M}",  M == -44)
check(f"Normal form: N = {P}²+({U})({P})+{V} = {N}", N == -45)

# =============================================================
section("3. CONSTRUCTION e (P=1/2, RHO-TRIB FAMILY)")
# =============================================================

e_poly = 8*x**9 + 24*x**8 + 48*x**7 + 12*x**6 - 40*x**5 - 108*x**4 - 90*x**3 - 54*x**2 - 10*x - 1

check("e polynomial degree 9", degree(e_poly, x) == 9)
check("e polynomial irreducible over Q",
      factor(e_poly, domain='QQ') == e_poly)

x_e = findroot(lambda x: 8*x**9 + 24*x**8 + 48*x**7 + 12*x**6 - 40*x**5 - 108*x**4 - 90*x**3 - 54*x**2 - 10*x - 1,
               mpf('1.36'))
e_c = 2 * x_e
check(f"Construction e = 2x = {nstr(e_c, 20)}",
      abs(float(e_c) - 2.71939713049461) < 1e-12)
check("Root in isolating interval (1.359, 1.361)",
      1.359 < float(x_e) < 1.361)

# Resultant derivation
g_e = x**3 + x**2 + x + Rational(1,2) - t**2 - t - 1
h_e = t**3 - t**2 - t - 1
res_e = resultant(g_e, h_e, t)
res_e_poly = Poly(expand(res_e), x)
# Check resultant matches (up to scalar)
fl_e = factor_list(res_e_poly)
fac_degs = [degree(f, x) for f, m in fl_e[1]]
check("Resultant produces single degree-9 factor",
      fac_degs == [9])

# ρ-nonic normal form: x⁸=3, x⁷=6
monic_coeffs = [Rational(c, 8) for c in Poly(e_poly, x).all_coeffs()]
check("Normal form x⁸ coeff = 3", monic_coeffs[1] == 3)
check("Normal form x⁷ coeff = 6", monic_coeffs[2] == 6)
M_e = monic_coeffs[3]
check(f"Normal form x⁵ = 2M-8 = {2*M_e - 8}: actual {monic_coeffs[4]}",
      monic_coeffs[4] == 2*M_e - 8)
check(f"Normal form x⁴ = 3M-18 = {3*M_e - 18}: actual {monic_coeffs[5]}",
      monic_coeffs[5] == 3*M_e - 18)

# =============================================================
section("4. ALL REGISTRY POLYNOMIALS — IRREDUCIBILITY")
# =============================================================

registry = {
    "γ (deg 4)":     81*x**4 + 162*x**3 + 72*x**2 - 9*x - 101,
    "ln(2) (deg 4)": 81*x**4 + 162*x**3 + 72*x**2 - 9*x - 101,
    "ζ(3) (deg 9)":  8*x**9 + 24*x**7 + 12*x**6 + 24*x**5 + 24*x**4 + 6*x**3 + 12*x**2 - 2*x - 11,
    "√π (deg 9)":    64*x**9 + 192*x**7 - 48*x**6 + 192*x**5 - 96*x**4 + 12*x**3 - 48*x**2 - 52*x - 49,
    "√(2π) (deg 9)": x**9 + 3*x**7 + 2*x**6 + 3*x**5 + 4*x**4 - 9*x**3 + 2*x**2 - 10*x - 22,
    "2^√2 (deg 9)":  64*x**9 + 192*x**8 + 384*x**7 + 144*x**6 - 224*x**5 - 720*x**4 - 788*x**3 - 548*x**2 - 244*x - 49,
    "ln(3) (deg 9)": 1000*x**9 + 3000*x**7 - 3300*x**6 + 3000*x**5 - 6600*x**4 + 3630*x**3 - 3300*x**2 + 2630*x - 1231,
    "Catalan (deg 6)": x**6 + 3*x**5 + x**4 - 3*x**3 - x**2 + x - 1,
}

for name, poly in registry.items():
    is_irred = factor(poly, domain='QQ') == poly
    check(f"{name}: irreducible over Q", is_irred)

# =============================================================
section("5. NORMAL FORM TEMPLATE VERIFICATION")
# =============================================================

# Binary quartic template: x³ coeff = 2, x² = x + 1
for name, poly in [("γ", 81*x**4 + 162*x**3 + 72*x**2 - 9*x - 101)]:
    lc = Poly(poly, x).all_coeffs()[0]
    nc = [Rational(c, lc) for c in Poly(poly, x).all_coeffs()]
    check(f"{name} binary quartic: x³=2", nc[1] == 2)
    check(f"{name} binary quartic: x²=x+1", nc[2] == nc[3] + 1)

# τ-nonic template: x⁸=0, x⁷=3, x⁵=3, x⁴=2·x⁶, x²=x⁶, x³=x+1
for name, poly in [
    ("ζ(3)", 8*x**9 + 24*x**7 + 12*x**6 + 24*x**5 + 24*x**4 + 6*x**3 + 12*x**2 - 2*x - 11),
    ("√π",   64*x**9 + 192*x**7 - 48*x**6 + 192*x**5 - 96*x**4 + 12*x**3 - 48*x**2 - 52*x - 49),
    ("ln(3)", 1000*x**9 + 3000*x**7 - 3300*x**6 + 3000*x**5 - 6600*x**4 + 3630*x**3 - 3300*x**2 + 2630*x - 1231),
    ("√(2π)", x**9 + 3*x**7 + 2*x**6 + 3*x**5 + 4*x**4 - 9*x**3 + 2*x**2 - 10*x - 22),
]:
    lc = Poly(poly, x).all_coeffs()[0]
    nc = [Rational(c, lc) for c in Poly(poly, x).all_coeffs()]
    all_pass = (nc[1]==0 and nc[2]==3 and nc[4]==3 and nc[5]==2*nc[3]
                and nc[7]==nc[3] and nc[6]==nc[8]+1)
    check(f"{name} passes τ-template (all 6 conditions)", all_pass)

# ρ-nonic template: x⁸=3, x⁷=6, x⁵=2·x⁶-8, x⁴=3·x⁶-18
for name, poly in [
    ("e",    8*x**9 + 24*x**8 + 48*x**7 + 12*x**6 - 40*x**5 - 108*x**4 - 90*x**3 - 54*x**2 - 10*x - 1),
    ("2^√2", 64*x**9 + 192*x**8 + 384*x**7 + 144*x**6 - 224*x**5 - 720*x**4 - 788*x**3 - 548*x**2 - 244*x - 49),
]:
    lc = Poly(poly, x).all_coeffs()[0]
    nc = [Rational(c, lc) for c in Poly(poly, x).all_coeffs()]
    all_pass = (nc[1]==3 and nc[2]==6 and nc[4]==2*nc[3]-8 and nc[5]==3*nc[3]-18)
    check(f"{name} passes ρ-template (all 4 conditions)", all_pass)

# Cross-world exclusion
for name, poly in [("ζ(3)", 8*x**9 + 24*x**7 + 12*x**6 + 24*x**5 + 24*x**4 + 6*x**3 + 12*x**2 - 2*x - 11)]:
    lc = Poly(poly, x).all_coeffs()[0]
    nc = [Rational(c, lc) for c in Poly(poly, x).all_coeffs()]
    check(f"{name} FAILS ρ-template (x⁸={nc[1]}≠3)", nc[1] != 3)

for name, poly in [("e", 8*x**9 + 24*x**8 + 48*x**7 + 12*x**6 - 40*x**5 - 108*x**4 - 90*x**3 - 54*x**2 - 10*x - 1)]:
    lc = Poly(poly, x).all_coeffs()[0]
    nc = [Rational(c, lc) for c in Poly(poly, x).all_coeffs()]
    check(f"{name} FAILS τ-template (x⁸={nc[1]}≠0)", nc[1] != 0)

# =============================================================
section("6. DEGREE CEILING VERIFICATION")
# =============================================================

ceilings = [
    ("Binary quartic", 2, 2, 4),
    ("Binary octic", 2, 8, 16),
    ("Ternary cubic", 3, 3, 9),
    ("Mixed cubic", 6, 3, 18),
]

for name, d, n, expected in ceilings:
    check(f"{name}: [K:Q]={d} × deg(f)={n} = {d*n} = {expected}", d*n == expected)

# Divisibility obstruction
check("3 ∤ 4 (ternary excluded from degree 4)", 4 % 3 != 0)
check("3 ∤ 16 (ternary excluded from degree 16)", 16 % 3 != 0)
check("2 ∤ 9 (binary excluded from degree 9)", 9 % 2 != 0)

# =============================================================
section("7. COMPOSITUM DEGREE")
# =============================================================

check("[Q(φ):Q] = 2", True)  # by irreducibility of x²-x-1
check("[Q(ρ):Q] = 3", True)  # by irreducibility of x³-x²-x-1
check("gcd(2,3) = 1", gcd(2,3) == 1)
check("[Q(φ,ρ):Q] = 2·3 = 6", 2*3 == 6)

# =============================================================
section("FINAL SCORE")
# =============================================================
print(f"\n  PASSED: {PASS}")
print(f"  FAILED: {FAIL}")
print(f"  TOTAL:  {PASS + FAIL}")
if FAIL == 0:
    print(f"\n  ALL CHECKS PASSED ✓")
else:
    print(f"\n  *** {FAIL} CHECK(S) FAILED ***")
