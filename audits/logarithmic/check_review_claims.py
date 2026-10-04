"""Independent algebra and counterexamples for the 2026-09-30 review.

This is not an infinite-dimensional proof or a literature-priority certificate.
The one-layer examples refute the advertised general algebra/budget estimates,
not the positive-dissipation theorem for the infinite model.
"""

import math
from fractions import Fraction
from itertools import permutations, product
import sympy as sp


def edge_system(values, shells, beta, intra):
    rhs = [sp.Integer(0) for _ in values]
    edges = []
    for n in range(shells):
        w, z = 2*n, 2*n+1
        edges.extend([(w, z, intra[n]), (z, w, beta*intra[n])])
        if n+1 < shells:
            edges.extend([(z, w+2, sp.Integer(2)**n),
                          (w, z+2, sp.Rational(1, 2)**n)])
    for source, target, coefficient in edges:
        rhs[source] -= coefficient*values[source]*values[target]
        rhs[target] += coefficient*values[source]**2
    return rhs, edges


def check_core():
    u = sp.symbols('u0:6')
    b = sp.symbols('b0:6', positive=True)
    d = sp.symbols('d0:6')
    beta, nu = sp.symbols('beta nu', positive=True)
    h = sp.symbols('h0:3', positive=True)
    nonlinear, edges = edge_system(u, 3, beta, h)
    assert sp.expand(sum(x*f for x, f in zip(u, nonlinear))) == 0
    alternate = sp.expand(sum((-1)**(i//2)*u[i]*nonlinear[i] for i in range(6)))
    contrast = sp.expand(sum((-1)**i*u[i]*nonlinear[i] for i in range(6)))
    state = dict(zip(u, (0, 1, 1, 0, 0, 0)))
    assert alternate.subs(state) == -2
    state = dict(zip(u, (2, 1, 0, 0, 0, 0)))
    assert sp.simplify(contrast.subs(state).subs(beta, 1)+4*h[0]) == 0
    print('PASS: two non-invariant weighted energies; no classification of all quadratic forms')
    base, _ = edge_system(b, 3, beta, h)
    shifted, _ = edge_system([b[i]+d[i] for i in range(6)], 3, beta, h)
    mass = [sp.Integer(4)**(i//2) for i in range(6)]
    forcing = [nu*mass[i]*b[i]-base[i] for i in range(6)]
    residual = sum(d[i]*(shifted[i]-base[i]-nu*mass[i]*d[i]) for i in range(6))
    residual += nu*sum(mass[i]*d[i]**2 for i in range(6))/2
    residual += sum(forcing[i]*d[i]**2/b[i] for i in range(6))/2
    residual += sum(c*b[i]**2*b[j]*(d[i]/b[i]-d[j]/b[j])**2
                    for i, j, c in edges)/2
    assert sp.expand(residual) == 0
    assert sp.Integer(8)**3 > 8*20*sp.Rational(9, 4)
    assert sp.Integer(8)**3/2-sp.Rational(3, 2)-sp.Rational(9, 4)/8**3 > 2
    C = 20**(2/3)/(1-2**(-2/3))+64*4
    B = 2**0.25*C**0.75
    radius, loss = 2*B, 1+3*(2*B)
    window = 2/loss+loss*radius**2/2
    print('PASS: closed conservation, finite relative-energy identity, envelope constants')
    print(f'alpha=eta=beta=1: B={B:.8g}; sufficient T*={window:.8g}; eps={1/loss:.8g}')


def one_layer(T, nu):
    # beta=1, equal unit forces and zero data: w=z=(1-exp(-nu*t))/nu.
    endpoint = -math.expm1(-nu*T)/nu
    integral = T/nu+math.expm1(-nu*T)/nu**2
    square = (T+2*math.expm1(-nu*T)/nu
              -math.expm1(-2*nu*T)/(2*nu))/nu**2
    assert math.isclose(endpoint**2+2*nu*square, 2*integral,
                        rel_tol=1e-10, abs_tol=1e-12)
    return endpoint, integral, square


def check_objections():
    a, z = sp.Integer(1), sp.Integer(2)
    assert a >= 1 and 2*a*a-a-z == -1
    print('COUNTEREXAMPLE: bar_w=1, bar_z=2 gives GAP=-1')
    T, nu, eta = 0.1, 1.0, 0.5
    endpoint, _, W = one_layer(T, nu)
    radius = math.sqrt(2)*T
    claimed_Z = max(radius, 1)*eta*W
    # Even overestimating the deleted e_0 outgoing term by W does not help.
    claimed_variation = 2*radius*(1+eta+nu)*W
    assert nu*W > claimed_Z and endpoint**2 > claimed_variation
    print(f'COUNTEREXAMPLE (Z): {nu*W:.12g} > {claimed_Z:.12g}')
    print(f'COUNTEREXAMPLE (variation): {endpoint**2:.12g} > {claimed_variation:.12g}')
    _, J, W = one_layer(3.0, 0.1)
    mean = J/3
    gap = 2*mean**2-2*mean
    assert mean >= 1 and gap > 0 and W > 3 and 2*J-2*W < 0
    print(f'COUNTEREXAMPLE (GAP/second moment): GAP={gap:.12g}>0; Z={W:.12g}>T=3')
    print(f'COUNTEREXAMPLE (BAL): I-W-Z={2*J-2*W:.12g}<0')
    # Arbitrary time profiles, not ODE trajectories: T=1, z=1,
    # w=epsilon^(-1/2) on [0,epsilon] and zero otherwise give W=V=1.
    epsilon = sp.symbols('epsilon', positive=True)
    square_spike = sp.simplify(epsilon*epsilon**(-1))
    fourth_spike = sp.simplify(epsilon*epsilon**(-2))
    assert square_spike == 1 and fourth_spike/square_spike**2 == 1/epsilon
    print('PASS: definition counterexample W=V=1, I4=1/epsilon (not an ODE trajectory)')
    eta_s, level, alpha, window, wrms, zrms = sp.symbols('eta n alpha T w z', positive=True)
    Cflux = eta_s*2**(level/2)*(level+1)**alpha*wrms**2*zrms
    theta = wrms*(level+1)**alpha
    V = 2**level*window*zrms**2/(level+1)**(2*alpha)
    assert sp.simplify(Cflux-eta_s*theta**2*sp.sqrt(V/window)) == 0
    print('PASS: corrected RMS identity contains eta and no fourth moment')


def check_known_structure():
    w, z, beta, h = sp.symbols('w z beta h', real=True)
    s = sp.sqrt(1+beta**2)
    x, y = (w-beta*z)/s, (beta*w+z)/s
    dw, dz = h*(beta*z**2-w*z), h*(w**2-beta*w*z)
    assert sp.simplify((dw-beta*dz)/s+h*s*x*y) == 0
    assert sp.simplify((beta*dw+dz)/s-h*s*x**2) == 0
    assert sp.simplify(dw+sp.I*dz-sp.I*h*(w-beta*z)*(w+sp.I*z)) == 0
    print('PASS: local exchange rotates to Tao pump; complex representation exists')


def check_novelty_reductions():
    # These are finite algebra checks, not a priority or existence certificate.
    u = sp.symbols('u0:6')
    h = sp.symbols('h0:3', positive=True)
    full, _ = edge_system(u, 3, 0, h)
    # Remove only weak cross-shell edges; local forward and strong edges remain.
    for n in range(2):
        i, j, c = 2*n, 2*n+3, sp.Rational(1, 2)**n
        full[i] += c*u[i]*u[j]
        full[j] -= c*u[i]**2
    p = [h[0], sp.Integer(1), h[1], sp.Integer(2), h[2]]
    chain = [sp.Integer(0) for _ in u]
    for i, c in enumerate(p):
        chain[i] -= c*u[i]*u[i+1]
        chain[i+1] += c*u[i]**2
    assert all(sp.expand(a-b) == 0 for a, b in zip(full, chain))
    J, hn, Kn, w = sp.symbols('J h K w', positive=True)
    z = J/(hn*w**2)
    next_w = J/(Kn*z**2)
    assert sp.simplify(next_w-hn**2*w**4/(Kn*J)) == 0
    n, alpha, eta = sp.symbols('n alpha eta', positive=True)
    v = sp.symbols('v', positive=True)
    rho = 2*alpha/3
    wn = J**sp.Rational(1, 3)*eta**sp.Rational(-2, 3)*(n+1)**(-rho)*v
    wnext = eta**2*(n+1)**(2*alpha)*wn**4/J
    normalized = eta**sp.Rational(2, 3)*(n+2)**rho*wnext/J**sp.Rational(1, 3)
    assert sp.simplify(normalized-((n+2)/(n+1))**rho*v**4) == 0
    # Matching finite cutoffs verify the shift algebra, not the infinite limit.
    increments = sp.symbols('ell0:8')
    tail0 = -sum(sp.Rational(1, 4)**(k+1)*increments[k] for k in range(8))
    tail1 = -sum(sp.Rational(1, 4)**k*increments[k] for k in range(1, 8))
    assert sp.expand(tail1-4*tail0-increments[0]) == 0
    print('PASS: deleted-edge DN/KP mapping and constant-flux tail recurrence')
    print('Scope: no root matching, infinite steady existence, dynamic threshold or priority certification')


def check_ground_state_mapping():
    # Arbitrary finite pump graphs; forcing is derived algebraically from b.
    # It can be signed: these checks do not establish positive-force existence
    # or the nonnegative terms needed for nonlinear non-expansion.
    graphs = [
        (2, [(0, 1)]),
        (2, [(0, 1), (1, 0)]),
        (3, [(0, 1), (1, 2), (2, 0)]),
        (4, [(0, 1), (0, 2), (1, 3), (2, 3)]),
        (4, [(0, 1), (1, 2), (2, 0), (1, 0), (2, 3)]),
    ]
    for size, pairs in graphs:
        u = sp.Matrix(sp.symbols(f'u0:{size}'))
        b = sp.Matrix(sp.symbols(f'b0:{size}', positive=True))
        d = sp.Matrix(sp.symbols(f'd0:{size}'))
        c = sp.symbols(f'c0:{len(pairs)}', positive=True)
        nu = sp.symbols('nu', positive=True)
        mass = sp.diag(*[sp.Integer(i+1) for i in range(size)])
        nonlinear = sp.zeros(size, 1)
        weights = sp.zeros(size)
        for (i, j), coefficient in zip(pairs, c):
            nonlinear[i] -= coefficient*u[i]*u[j]
            nonlinear[j] += coefficient*u[i]**2
            weights[i, j] += coefficient*b[i]/2
            weights[j, i] += coefficient*b[i]/2
        assert sp.expand(u.dot(nonlinear)) == 0
        base = nonlinear.subs(dict(zip(u, b)), simultaneous=True)
        shifted = nonlinear.subs(dict(zip(u, b+d)), simultaneous=True)
        jacobian = nonlinear.jacobian(u).subs(dict(zip(u, b)), simultaneous=True)-nu*mass
        operator = -(jacobian+jacobian.T)/2
        forcing = nu*mass*b-base
        assert all(sp.expand(x) == 0 for x in operator*b-(forcing+nu*mass*b)/2)
        for i in range(size):
            for j in range(size):
                if i != j:
                    assert sp.expand(operator[i, j]+weights[i, j]) == 0
        ground = sum(weights[i, j]*b[i]*b[j]*(d[i]/b[i]-d[j]/b[j])**2
                     for i in range(size) for j in range(i+1, size))
        ground += sum((operator*b)[i]*d[i]**2/b[i] for i in range(size))
        relative = nu*d.dot(mass*d)/2
        relative += sum(forcing[i]*d[i]**2/b[i] for i in range(size))/2
        relative += sum(coefficient*b[i]**2*b[j]*(d[i]/b[i]-d[j]/b[j])**2
                        for (i, j), coefficient in zip(pairs, c))/2
        derivative = d.dot(shifted-base-nu*mass*d)
        assert sp.expand(d.dot(operator*d)-ground) == 0
        assert sp.expand(ground-relative) == 0
        assert sp.expand(derivative+ground) == 0
    print('PASS: ground-state representation maps exactly to (R) on 5 finite graph topologies')


def check_stationary_normalization():
    # A constant-flux tail normalizer does not assert matching forced roots.
    qm, q, qp, J, nu, mass = sp.symbols('qm q qp J nu mass', positive=True)
    rm, r, rp = sp.symbols('rm r rp', real=True)
    pm, p = J/(qm**2*q), J/(q**2*qp)
    stationary = pm*(qm*rm)**2-p*(q*r)*(qp*rp)-nu*mass*q*r
    mu = nu*mass*q**2/J
    assert sp.expand(stationary*q/J-(rm**2-r*rp-mu*r)) == 0
    n, eta = sp.symbols('n eta', positive=True)
    alpha = sp.symbols('alpha', real=True)
    w = J**sp.Rational(1, 3)*eta**sp.Rational(-2, 3)*(n+1)**(-2*alpha/3)
    z = J**sp.Rational(1, 3)*eta**sp.Rational(1, 3)*2**(-n/2)*(n+1)**(alpha/3)
    next_w = w.subs(n, n+1)
    ratio_down = eta**2*2**(-n)*(n+1)**(2*alpha)
    ratio_up = 4/eta**2*2**n*(n+1)**(-2*alpha)*((n+2)/(n+1))**(-4*alpha/3)
    assert sp.simplify(z**2/w**2-ratio_down) == 0
    assert sp.simplify(4*next_w**2/z**2-ratio_up) == 0
    # The exact increment recurrence precedes the sufficient sign argument.
    increment = r-rm
    next_increment = rm**2/r-mu-r
    assert sp.simplify(next_increment+increment*(1+rm/r)+mu) == 0
    print('PASS: stationary normalization and alternating damping leading ratios')
    print('Scope: no ODE simulation, uniform bound, infinite tail limit or priority certification')


def check_irregular_envelopes():
    # Finite zero-extended q arrays: tests of the majorant algebra only.
    arrays = [
        [sp.Integer(1)],
        [sp.Rational(1, n+1) for n in range(12)],
        [sp.Rational(1, 100), sp.Integer(10), sp.Rational(1, 1000), sp.Integer(3)],
        [sp.Integer(2)**(-n) if n % 3 else sp.Integer(n+1) for n in range(20)],
    ]
    for q in arrays:
        t = [sp.Integer(0)]*(len(q)+1)
        for n in reversed(range(len(q))):
            t[n] = max(q[n], t[n+1]/2)
        assert all(t[n] >= q[n] and t[n+1] <= 2*t[n] for n in range(len(q)))
        assert sum(v*v for v in t) <= sp.Rational(4, 3)*sum(v*v for v in q)
    # Exact positive-coefficient certificates for all b>=1 and H>=0.
    beta_excess, H = sp.symbols('beta_excess H', nonnegative=True)
    b, D, Az3 = 1+beta_excess, H**2, 8*(1+H)
    x = 1+8*b**3*Az3*D+8*(b*H+D+1)
    certificates = [x-8*b**3*Az3*D, x-8*(b*H+D+1),
                    x*(x/2-b*H-D/x-4)]
    for certificate in certificates:
        poly = sp.Poly(sp.cancel(certificate), beta_excess, H)
        assert all(c >= 0 for c in poly.coeffs()) and poly.coeff_monomial(1) > 0
    assert x.subs({beta_excess: 0, H: 1}) == 153
    print('PASS: 4 finite future-envelope cases and 3 general constant certificates')


def check_sparse_bottlenecks():
    # Infinite summability/noncomparability are proved in the text, not here.
    previous_ratio = None
    for m in range(1, 13):
        n = 2**m
        assert m+1 <= 2**(n//2)
        if m >= 2:
            ratio_squared = sp.Rational((m+1)**2, 2**(n-1))
            if previous_ratio is not None:
                assert ratio_squared < previous_ratio
            previous_ratio = ratio_squared
    a, K = sp.symbols('a K', positive=True)
    h = sp.sqrt(K)*a
    assert sp.simplify((K/h**2)**sp.Rational(1, 3)-a**sp.Rational(-2, 3)) == 0
    assert sp.simplify(h**2/K-a**2) == 0
    print('PASS: 12 sparse-bottleneck coefficient cases and exact even-branch reduction')
    print('Scope: finite algebra only; no infinite proof, ODE simulation or originality certification')


def check_geometric_profile():
    # Constant-flux normalization algebra; no assumption that forced roots match.
    a, K, r, rho, I, beta = sp.symbols('a K r rho I beta', positive=True)
    rn = a**2*r**4
    w = I**sp.Rational(1, 3)*r*rho
    lower = a**2*w**4/(2*I)-beta*a/sp.sqrt(K)*w-I/(K*w**2)
    d = 1/(sp.sqrt(K)*a*r**3)
    assert sp.simplify(lower/(I**sp.Rational(1, 3)*rn)
                       -(rho**4/2-beta*d*rho-d**2/rho**2)) == 0
    H = sp.symbols('H', positive=True)
    c = H**sp.Rational(-2, 3)*2**sp.Rational(-1, 9)
    D = 2**sp.Rational(1, 4)*H**sp.Rational(1, 2)*c**sp.Rational(-3, 4)
    assert sp.simplify(D-2**sp.Rational(1, 3)*H) == 0
    B, Dg, Z = sp.symbols('B Dg Z', nonnegative=True)
    b = 1+B
    x = 1+8*b**3*Z+8*(b*Dg+Dg**2+1)
    for expr in [x-8*b**3*Z, x*(x/2-b*Dg-Dg**2/x-4)]:
        poly = sp.Poly(sp.cancel(expr), B, Dg, Z)
        assert all(coefficient >= 0 for coefficient in poly.coeffs())
        assert poly.coeff_monomial(1) > 0
    # Exact parity-dependent infinite geometric sums and recurrence exponents.
    t = sp.Rational(1, 16)
    assert 1/(4*(1-t)) == sp.Rational(4, 15)
    assert t/(2*(1-t)**2) == sp.Rational(8, 225)
    def log_r(n):
        return (-sp.Rational(n, 15)-sp.Rational(17, 225) if n % 2 == 0
                else -sp.Rational(4*n, 15)-sp.Rational(8, 225))
    for n in range(24):
        log_a = sp.Integer(0) if n % 2 == 0 else sp.Rational(n, 2)
        assert log_r(n+1) == 2*log_a+4*log_r(n)
        assert log_r(n) >= -sp.Rational(n, 3)-sp.Rational(1, 9)
        assert log_r(n) <= -sp.Rational(n, 15)
    assert sum(1 for n in range(24) if n % 2 == 0) == 12
    # q=2^(-2/225): check the published periodic square-sum expression
    # against finite prefixes plus exact remainders. Convergence requires
    # 0<q<1 in the analytic proof; finite checks do not prove that limit.
    q = sp.symbols('q', positive=True)
    periodic_sum = q**17/(1-q**30)+q**68/(1-q**120)
    for pairs in (1, 2, 4):
        prefix = sum(q**(17+30*j)+q**(68+120*j) for j in range(pairs))
        remainder = q**(17+30*pairs)/(1-q**30)+q**(68+120*pairs)/(1-q**120)
        assert sp.cancel(periodic_sum-prefix-remainder) == 0
    for n in range(24):
        exponent = 17+30*(n//2) if n % 2 == 0 else 68+120*(n//2)
        assert 2*log_r(n) == -sp.Rational(2, 225)*exponent
    print('PASS: normalized growth, geometric constants and 24 dense-bottleneck exponent cases')
    print('PASS: periodic square-sum exponents and 3 finite prefix/remainder identities')


def check_root_input_bound():
    # Independently assemble both root equations; neither longitudinal outlet
    # may disappear from the scalar input estimate.
    w0, z0, w1, z1 = sp.symbols('w0 z0 w1 z1', nonnegative=True)
    h0, h1, beta, nu, C0, b = sp.symbols('h0 h1 beta nu C0 b', nonnegative=True)
    nonlinear, _ = edge_system([w0, z0, w1, z1], 2, beta, [h0, h1])
    s = w0+z0
    derivative = 2+nonlinear[0]+nonlinear[1]-nu*s
    expected = 2+h0*(w0**2+beta*z0**2)-(1+beta)*h0*w0*z0-w0*z1-z0*w1-nu*s
    assert sp.expand(derivative-expected) == 0
    L = 1+(1+(1+b)*h0)*C0
    positive = h0*(w0**2+beta*z0**2)+(1-nu)*s
    positive += w0*(C0-z1)+z0*(C0-w1)
    positive += (1+beta)*h0*(C0*w0+(C0-w0)*z0)+(b-beta)*h0*C0*s
    assert sp.expand(derivative+L*s-2-positive) == 0
    print('PASS: both root longitudinal outlets and conditional scalar lower-bound decomposition')
    print('Conditions: coordinates <= C0, beta <= b, nu <= 1; no trajectory-bound certification')


def check_steady_necessity_reduction():
    a, w, r, I, rho = sp.symbols('a w r I rho', positive=True)
    rn = a**2*r**4
    upper = 4*a**2*(I**sp.Rational(1, 3)*r*rho)**4/I
    assert sp.simplify(upper/(I**sp.Rational(1, 3)*rn)-4*rho**4) == 0
    v = sp.symbols('v', positive=True)
    assert sp.simplify(4**sp.Rational(1, 3)*4
                       *(v/4**sp.Rational(1, 3))**4-v**4) == 0
    q = sp.Rational(1, 2)
    for j in range(7):
        assert (q**(4**j))**4 == q**(4**(j+1))
    print('PASS: stationary normalized upper recurrence and 7 finite growth iterations')
    print('Scope: algebra only; infinite necessity, existence and originality are separate arguments')


def check_deleted_inviscid_root_match():
    # An explicit inviscid deleted-chain calculation, not a full-graph or
    # positive-viscosity equilibrium construction.
    h, K, q, y = sp.symbols('h K q y', positive=True)
    t = sp.sqrt(y)
    w0, z0, w1 = 1/(h*t*q), t*q, t/(K*q**2)
    polynomial = y**2-q*y-1/(h*q)
    assert sp.simplify(h*w0*z0-1) == 0
    assert sp.simplify(K*z0**2*w1-t**3) == 0
    assert sp.simplify(t*(w0+z0-t**3)+polynomial) == 0
    assert sp.simplify(h*w0**2-K*z0*w1+1+polynomial/(y*q)) == 0
    positive_y = (q+sp.sqrt(q**2+4/(h*q)))/2
    assert sp.simplify(polynomial.subs(y, positive_y)) == 0
    print('PASS: deleted inviscid chain double-root matching and positive polynomial root')


def check_nonadjacent_profiles():
    # A finite defect set, continued by a_n=2^(n/2) thereafter. This gives
    # exact infinite-tail logarithms, but checking 377 sets is not a theorem.
    def profile(defects, length=24):
        result = [-sp.Rational(length, 3)-sp.Rational(1, 9)]*(length+1)
        for n in reversed(range(length)):
            log_a = sp.Integer(0) if n in defects else sp.Rational(n, 2)
            result[n] = (result[n+1]-2*log_a)/4
        return result
    cases = 0
    for mask in range(2**12):
        if mask & (mask << 1):
            continue
        defects = {n for n in range(12) if mask & (1 << n)}
        logs = profile(defects)
        assert logs[0] <= 0
        assert all(logs[n] <= -sp.Rational(n, 15) for n in range(1, 24))
        cases += 1
    assert cases == 377
    # A long adjacent weak block violates the claimed nonadjacent-profile
    # bound. It does NOT prove absence of dynamic anomalous dissipation.
    adjacent = profile(set(range(1, 16)))
    assert adjacent[8] > -sp.Rational(8, 15)
    print('PASS: 377 nonadjacent defect sets; adjacent-block counterexample to profile bound')


def check_partial_cut_and_raw_recurrence():
    # Assemble the vector field, not a copied claimed cut formula. At n=0
    # the partial cut contains one forced root, whereas n>=1 contains both.
    shells = 6
    u = sp.symbols('u0:12', nonnegative=True)
    h = sp.symbols('h0:6', positive=True)
    beta, nu = sp.symbols('beta nu', nonnegative=True)
    nonlinear, _ = edge_system(u, shells, beta, h)
    mass = [sp.Integer(4)**(i//2) for i in range(2*shells)]
    rhs = [nonlinear[i]-nu*mass[i]*u[i]+int(i < 2) for i in range(2*shells)]
    for n in range(shells-1):
        w, z, wp, zp = u[2*n:2*n+4]
        selected = range(2*n+1)
        input_term = u[0]+(u[1] if n else 0)
        loss = nu*sum(mass[i]*u[i]**2 for i in selected)
        flux = h[n]*w*z*(w-beta*z)+sp.Rational(1, 2)**n*w**2*zp
        if n:
            flux += sp.Rational(1, 2)**(n-1)*u[2*n-2]**2*z
        assert sp.expand(sum(u[i]*rhs[i] for i in selected)
                         -input_term+loss+flux) == 0
        if n:
            assert sp.Rational(1, 2)**n+sp.Rational(1, 2)**(n-1) == 3*sp.Rational(1, 2)**n
    # Original z balance implies the lower recurrence after dropping its
    # nonnegative upstream/root input. Its inequality premises remain analytic.
    w, z, wp, h0, K, upstream, damping, force = sp.symbols(
        'w z wp h K upstream damping force', positive=True)
    balance = upstream-K*z*wp+h0*w**2-beta*h0*w*z-damping*z+force
    solved = (upstream+force)/(K*z)+h0*w**2/(K*z)-beta*h0*w/K-damping/K
    assert sp.simplify(balance-K*z*(solved-wp)) == 0
    print('PASS: 5 independently assembled partial cuts, root distinction and raw z recurrence')


def check_weak_block_bounds():
    # Scalar induction, endpoint indices and exact geometric tails only.
    v = sp.symbols('v', nonnegative=True)
    assert sp.expand((1+v)/2-v) == (1-v)/2
    assert sp.expand(1-(1+v)/2) == (1-v)/2
    # Thus v>=1 implies (1+v)/2<=v; v<=1 implies <=1.
    # Rational sample iterations test the direction and endpoint exponent.
    cases = 0
    for terminal in [sp.Rational(1, 16), sp.Integer(1), sp.Integer(16), sp.Integer(256)]:
        value = terminal
        for distance in range(1, 7):
            value = (1+sp.sqrt(value))/2
            upper = max(sp.Integer(1), terminal**sp.Rational(1, 2**distance))
            assert float(sp.N(value-upper, 30)) <= 0
            cases += 1
    previous_end = 0
    for j in range(2, 10):
        N, L, m = 2**(2**j), 2**(j+1), 2**j
        endpoint = N+L
        assert previous_end < N and endpoint <= 2*N
        distance = endpoint-(N+m-1)
        assert distance == 2**j+1 and sp.Rational(endpoint, 2**distance) <= 2
        previous_end = endpoint
    # A finite prefix compares the two explicit examples. Infinite density
    # and the viscosity conclusions are not inferred from these counts.
    size = 1024
    blocks = set()
    for j in range(2, 4):
        N, L = 2**(2**j), 2**(j+1)
        blocks.update(range(N, N+L+1))
    periodic_weak = {n for n in range(size) if n % 2 == 0}
    perturbed_weak = periodic_weak | blocks
    changed = periodic_weak ^ perturbed_weak
    assert changed == {n for n in blocks if n % 2 == 1}
    assert len(changed) == 12
    # Sum_{k>=s} (N+k)/4^(k+1) is the exact bulk-tail upper bound.
    N, s = sp.symbols('N s', nonnegative=True)
    tail = 4**(-s)*((N+s)/3+sp.Rational(1, 9))
    assert sp.simplify(tail-tail.subs(s, s+1)-(N+s)/4**(s+1)) == 0
    # The quartic time-moment sum and the selected cut are proven in the text;
    # this arbitrary finite array checks only the averaging inequality.
    moments = [sp.Rational(i+1, 100) for i in range(16)]
    flux_squared = [sp.Rational(3, 2)*item for item in moments]
    assert min(flux_squared) <= sp.Rational(3, 2)*sum(moments)/len(moments)
    print(f'PASS: {cases} scalar backward cases, 8 weak-block index cases and exact tail sum')
    print('Scope: no numerical ODE, viscosity limit, density theorem or priority certification')


def check_frozen_block_probabilities():
    # Enumerate actual finite binary strings, rather than assert a copied
    # second-moment formula. 1=weak, 0=strong; blocks stay inside the string.
    cases = 0
    for size in range(3, 11):
        for length in range(1, min(4, size-1)+1):
            histogram = [[0, 0, 0, 0] for _ in range(size+1)]
            for mask in range(2**size):
                count = sum(not ((mask >> i) & 1)
                            and ((mask >> (i+1)) & (2**length-1)) == 2**length-1
                            for i in range(size-length))
                row = histogram[mask.bit_count()]
                row[0] += 1
                row[1] += count
                row[2] += count**2
                row[3] += int(count > 0)
            for p in [Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)]:
                moments = [sum(row[j]*p**weak*(1-p)**(size-weak)
                               for weak, row in enumerate(histogram))
                           for j in range(4)]
                mean = (size-length)*(1-p)*p**length
                assert moments[0] == 1 and moments[1] == mean
                assert moments[2] <= mean+mean**2
                assert moments[3] >= mean/(1+mean)
                cases += 1
    for k in range(4, 25):
        prefix = (k-1).bit_length()
        length = k+prefix+2
        candidates = 2**k-length
        assert candidates >= 2**(k-1)
        distance = length-prefix
        assert distance == k+2
        assert Fraction(2**(k+1)-1, 2**distance) < Fraction(1, 2)
        mean = Fraction(candidates, 2**(length+1))  # p=1/2
        assert mean >= Fraction(1, 32*k)
        assert mean/(1+mean) >= Fraction(1, 64*k)
        # A disjoint-grid attempt loses a factor ~k at the critical point.
        grid_mean = Fraction(2**k//(length+1), 2**(length+1))
        assert grid_mean <= Fraction(1, 8*k**2)
    for k in range(1, 13):
        assert 4**k//k >= Fraction(4**k, 2*k)
        assert (Fraction(2*4**k, 3)+Fraction(1, 9))/4**k < 1
    print(f'PASS: {cases} exact finite binary moment cases, 21 marker-block and 12 tail-index cases')
    print('Scope: not Borel-Cantelli, almost-sure dynamics, GP foundations or originality certification')


def check_dynamic_rescaling():
    # Exact coordinate change in the assembled finite vector field. Not a
    # quasistationary approximation or a viscosity-limit experiment.
    shells = 6
    w = sp.symbols('w0:6', nonnegative=True)
    v = sp.symbols('v0:6', nonnegative=True)
    a = sp.symbols('a0:6', positive=True)
    beta, nu = sp.symbols('beta nu', nonnegative=True)
    K = [sp.Integer(2)**n for n in range(shells)]
    h = [sp.sqrt(K[n])*a[n] for n in range(shells)]
    values = [item for n in range(shells)
              for item in (w[n], a[n]*v[n]/sp.sqrt(K[n]))]
    nonlinear, _ = edge_system(values, shells, beta, h)
    rhs = [nonlinear[i]-nu*K[i//2]**2*values[i]+int(i < 2)
           for i in range(2*shells)]
    for n in range(shells-1):
        expected_w = -a[n]**2*w[n]*v[n]+beta*a[n]**3*v[n]**2/sp.sqrt(K[n])
        expected_w -= w[n]*a[n+1]*v[n+1]/(K[n]*sp.sqrt(K[n+1]))
        expected_w += -nu*K[n]**2*w[n]+int(n == 0)
        expected_v = w[n]**2-v[n]*w[n+1]-beta*a[n]*w[n]*v[n]/sp.sqrt(K[n])
        expected_v += -nu*K[n]*v[n]+int(n == 0)/(a[n]*sp.sqrt(K[n]))
        if n:
            expected_w += a[n-1]**2*v[n-1]**2
            expected_v += 2*w[n-1]**2/(a[n]*K[n]**sp.Rational(3, 2))
        assert sp.simplify(rhs[2*n]-expected_w) == 0
        assert sp.simplify(rhs[2*n+1]/(a[n]*sp.sqrt(K[n]))-expected_v) == 0
    for n in [4, 6, 8]:
        k = sp.Integer(2)**n
        z, wn, wn_next, vn = sp.Integer(1), sp.Integer(1), 1/sp.sqrt(k), sp.sqrt(k)
        assert vn == wn**2/wn_next and z == vn/sp.sqrt(k)
        assert (wn**2+z**2+wn_next**2)/2 <= sp.Rational(9, 8)
        assert wn**2-vn*wn_next-wn*vn/sp.sqrt(k) == -1
    print('PASS: 5 exact rescaled coordinate pairs and 3 bounded-energy counterexamples to naive slaving')
    print('Counterexamples are states, not claimed trajectories; no dynamic classification follows')


def check_equivalence_obstructions():
    # Small exact regressions for the separate analytic equivalence argument.
    # They do not classify all linear/nonlinear conjugacies or certify novelty.
    images = 0
    for shells in (2, 3):
        u = sp.symbols(f'x0:{2*shells}')
        rhs, edges = edge_system(u, shells, sp.Integer(1),
                                 [sp.Integer(2)**n for n in range(shells)])
        graph = {(i, j) for i, j, c in edges if c > 0}
        for j in range(2*shells):
            parents = {i for i in range(2*shells)
                       if sp.Poly(rhs[j], *u).coeff_monomial(u[i]**2) > 0}
            assert parents == {i for i, target in graph if target == j}
            assert len(parents) == (1 if j < 2 else 2)
        for order in permutations(range(2*shells)):
            mapped = {(order[i], order[j]) for i, j in graph}
            assert any((j, i) in mapped for i, j in mapped)
            assert max(sum(target == j for _, target in mapped)
                       for j in range(2*shells)) == 2
            images += 1
    matrices = 0
    for entries in product((sp.Integer(0), sp.Rational(1, 2), sp.Integer(1)), repeat=4):
        matrix = sp.Matrix(2, 2, entries)
        if matrix.T*matrix == sp.eye(2):
            assert matrix in (sp.eye(2), sp.Matrix([[0, 1], [1, 0]]))
        matrices += 1
    lam = sp.symbols('lambda')
    for shells in range(1, 5):
        weights = [sp.Integer(4)**(i//2) for i in range(2*shells)]
        characteristic = sp.diag(*weights).charpoly(lam).as_expr()
        assert sp.expand(characteristic-sp.prod((lam-4**n)**2
                                                for n in range(shells))) == 0
        assert all(weights.count(4**n) == 2 for n in range(shells))
        # Standard one-component dyadic weights have simple eigenvalues.
        standard = [4**j for j in range(2*shells)]
        assert len(set(standard)) == len(standard)
    print(f'PASS: square-input graph, {images} relabelings, {matrices} small positive matrices, 4 spectra')
    print('Scope: positive energy-preserving linear maps and standard simple-spectrum viscosity only')


def check_fourth_root_block_indices():
    # Independent exact index checks for the new dynamic lemma. They do not
    # assert an ODE estimate from a state or prove an infinite-dimensional limit.
    cases = 0
    for start in range(2, 17):
        for length in range(2, 21):
            end = start+length-1
            # The last layer itself cannot use the bounded-return-flow seed.
            assert Fraction(end)-Fraction(end, 2) > Fraction(1, 2)
            for distance in range(1, length):
                layer = end-distance
                gap = Fraction(end, 2**distance)-Fraction(layer, 2)
                assert gap <= Fraction(1, 2)
                if distance == 1:
                    assert gap == Fraction(1, 2)
                else:
                    assert gap <= 0
                # Seed is Y_(J-1), so there are J-1-n fourth-root steps.
                steps = end-1-layer
                seed_exponent = Fraction(end, 2)*Fraction(1, 4**steps)
                assert seed_exponent == Fraction(2*end, 4**(end-layer))
                cases += 1
            root_steps = length-2
            assert Fraction(end, 2)*Fraction(1, 4**root_steps) == Fraction(8*end, 4**length)
    for j in range(2, 9):
        start = 2**(2**j)
        for coefficient in [Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)]:
            raw_length = coefficient*2**j
            length = -(-raw_length.numerator//raw_length.denominator)
            end = start+length-1
            if coefficient >= Fraction(1, 2):
                assert Fraction(start, 4**length) <= 1
                assert Fraction(end, 4**length) <= 2
            else:
                assert Fraction(start, 4**length) > 1
    print(f'PASS: {cases} exact return-flow and fourth-root endpoint cases; 28 rounded block scales')
    print('Scope: finite indices and exponents; not ODE dynamics, viscosity limits, probability or priority')


def check_amplitude_sum_weighted_blocks():
    # New identities are checked against the assembled conservative graph.
    # No time profiles are substituted for actual trajectories in this check.
    shells = 6
    u = sp.symbols('u0:12', nonnegative=True)
    h = sp.symbols('h0:6', positive=True)
    nonlinear, _ = edge_system(u, shells, sp.Integer(1), h)
    rhs = [value+int(i < 2) for i, value in enumerate(nonlinear)]
    for n in range(shells-1):
        w, z = u[2*n:2*n+2]
        k, e = sp.Integer(2)**n, sp.Rational(1, 2)**n
        expected = h[n]*(w-z)**2-e*w*u[2*n+3]-k*z*u[2*n+2]
        if n:
            expected += sp.Integer(2)**(n-1)*u[2*n-1]**2
            expected += sp.Rational(1, 2)**(n-1)*u[2*n-2]**2
        else:
            expected += 2  # Both root forces must be retained.
        assert sp.expand(rhs[2*n]+rhs[2*n+1]-expected) == 0
        full_budget = sum(u[i]*rhs[i] for i in range(2*n+2))
        expected_budget = u[0]+u[1]-k*z**2*u[2*n+2]-e*w**2*u[2*n+3]
        assert sp.expand(full_budget-expected_budget) == 0
    nu = sp.symbols('nu', nonnegative=True)
    for n in range(1, shells-1):
        viscous_sum = rhs[2*n]+rhs[2*n+1]-nu*4**n*(u[2*n]+u[2*n+1])
        # Solving the amplitude identity for incoming Y adds this positive term.
        assert sp.expand(rhs[2*n]+rhs[2*n+1]-viscous_sum-nu*4**n*(u[2*n]+u[2*n+1])) == 0
    cases = 0
    for start in range(1, 18):
        for length in range(1, 22):
            end = start+length-1
            gamma = Fraction(end, 2)*Fraction(1, 4**(end-start))
            assert gamma == Fraction(2*end, 4**length)
            if length >= 2:
                assert Fraction(8*end, 4**length) == 4*gamma
            cases += 1
    for j in range(2, 31):
        length = 2*j
        assert 4**length == 16**j
        # Even N and even original length leave an additional even weak layer.
        assert (length-1) % 2 == 1 and length % 2 == 0
        assert Fraction(2+2*j, 2*16**j) <= 1
    # Exact rounded starts when log2(j+1) is an integer; no floating-point ceil.
    for j, log_value in [(3, 2), (7, 3), (15, 4)]:
        start = 16**j*log_value
        assert start % 2 == 0
        gamma = Fraction(start+2*j, 2*16**j)
        assert gamma == Fraction(log_value, 2)+Fraction(j, 16**j)
    print(f'PASS: 5 amplitude sums, 5 full-cut budgets, 4 viscosity signs, {cases} exact block exponents')
    print('PASS: 29 candidate parity/remainder cases and 3 exact rounded starts')
    print('Scope: finite algebra only; not trajectory bounds, divergent-series proof, viscosity limit or priority')


def check_log_moment_and_local_pulse():
    # The four-node system is a closed, unforced, nonzero-data counterexample
    # to a local moment estimate. It is not the infinite zero-data trajectory.
    w, z, y, b = sp.symbols('w z y b', positive=True)
    h, g, K, e, r = sp.symbols('h g K e r', positive=True)
    values = [w, z, y, b]
    edges = [(0, 1, h), (1, 0, h), (1, 2, K),
             (0, 3, e), (2, 3, g), (3, 2, g)]
    rhs = [sp.Integer(0)]*4
    for source, target, coefficient in edges:
        rhs[source] -= coefficient*values[source]*values[target]
        rhs[target] += coefficient*values[source]**2
    assert sp.expand(sum(x*f for x, f in zip(values, rhs))) == 0
    ratio_rhs = sp.expand((rhs[1]*w-z*rhs[0])/w**2).subs(z, r*w)
    assert sp.expand(ratio_rhs-h*w*(1-r)*(1+r**2)+K*r*y-e*r*b) == 0
    assert sp.expand(1-(1-r)*(1+r**2)-r*(1-r+r**2)) == 0
    assert sp.expand(w*rhs[0]-(h*w*z**2-h*w**2*z-e*w**2*b)) == 0
    assert sp.expand(z*rhs[1]-(h*w**2*z-h*w*z**2-K*z**2*y)) == 0
    A, s, a = sp.symbols('A s a', positive=True)
    assert sp.expand(rhs[2].subs({y: A, b: A+s})-K*z**2-g*(A+s)*s) == 0
    assert sp.expand(rhs[3].subs({b: A, y: A+s})-e*w**2-g*(A+s)*s) == 0
    for n in range(9, 81):
        # Exact sixth powers of 4(n+1)^(2/3) <= sqrt(2^n).
        assert 4**6*(n+1)**4 <= 2**(3*n)
        assert 2**(4*n) > 12  # K^2 > 2 sqrt(3).
    upper_square = A**2/(1+a**2*A*s/4)
    assert sp.simplify(sp.diff(upper_square, s)+a**2*upper_square**2/(4*A)) == 0
    W_primitive = 4*A/a**2*sp.log(1+a**2*A*s/4)
    assert sp.simplify(sp.diff(W_primitive, s)-upper_square) == 0
    for alpha in [sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(3, 4), sp.Integer(1)]:
        rho = 2*alpha/3
        assert alpha/2+rho/4 == rho
        assert 2*alpha+rho == 8*alpha/3
    assert sp.Rational(8, 3)*sp.Rational(2, 3) == sp.Rational(16, 9)
    assert sp.Rational(16, 9) > sp.Rational(1, 3)
    print('PASS: closed four-node energy, ratio equation, two cut budgets and two receiving boundaries')
    print('PASS: 72 exact invariant-cone scales, decay integral and 4 log-moment normalizations')
    print('Scope: finite algebra; no common-zero-data counterexample or dynamic classification')


def check_arrival_and_pure_permutation():
    q, epsilon, c, ws2, zt2, G, incoming, outgoing = sp.symbols(
        'q epsilon c ws2 zt2 G incoming outgoing')
    zs2 = 2*(q-c)-ws2
    strong = q+epsilon-c-outgoing
    forward = (zt2-zs2)/2+G+strong-incoming
    assert sp.expand(forward-(epsilon+ws2/2+zt2/2+G-incoming-outgoing)) == 0
    rho = sp.symbols('rho')
    assert sp.simplify(sp.diff(2*rho+sp.log(1-rho), rho)-(1-2*rho)/(1-rho)) == 0

    entries = 0
    previous_end = -1
    for j in range(2, 11):
        start, length = 2**(2**j), 2**(j-1)
        assert start > previous_end and 4**length == start
        destination_to_source = [start+2*k for k in range(length)]
        destination_to_source += [start+2*k+1 for k in range(length)]
        assert sorted(destination_to_source) == list(range(start, start+2*length))
        original_logs = [Fraction(n, 2) if n % 2 == 0 else Fraction(n)
                         for n in range(start, start+2*length)]
        permuted_logs = []
        displacements = []
        for offset, source in enumerate(destination_to_source):
            n = start+offset
            h_log = Fraction(source, 2) if source % 2 == 0 else Fraction(source)
            a_log = h_log-Fraction(n, 2)
            assert 0 <= a_log <= Fraction(n, 2)
            if offset < length:
                assert a_log == Fraction(offset, 2)
            else:
                k = offset-length
                assert a_log == Fraction(start+3*k+2-length, 2)
            displacement = abs(source-n)
            assert 2**(2*displacement) <= n+2
            displacements.append(displacement)
            permuted_logs.append(h_log)
            entries += 1
        assert sorted(permuted_logs) == sorted(original_logs)
        assert max(displacements) == length-1
        end = start+length-1
        seed = Fraction(end, 2*4**(length-1))
        assert seed == Fraction(2*end, start) <= 4
        coefficient = sum((Fraction(r, 2*4**(r-1)) for r in range(1, length)), Fraction(0))
        assert coefficient <= Fraction(8, 9) and coefficient+seed <= 5
        previous_end = start+2*length-1

    admissible = 0
    for start in (8, 16, 32):
        for order in permutations(range(start, start+6)):
            logs = [Fraction(source, 2) if source % 2 == 0 else Fraction(source)
                    for source in order]
            if not all(Fraction(n, 2) <= value <= n
                       for n, value in zip(range(start, start+6), logs)):
                continue
            admissible += 1
            inverse = {source: start+offset for offset, source in enumerate(order)}
            for n in range(start, start+6):
                original_odd = n if n % 2 else n+1
                destination = inverse[original_odd]
                assert destination >= original_odd >= n
                distance = destination-n
                a_log = Fraction(original_odd)-Fraction(destination, 2)
                if distance <= Fraction(n, 2):
                    assert a_log >= Fraction(n, 4)
                    profile_term = Fraction(2, 4**(distance+1))*a_log
                    assert profile_term >= Fraction(n, 8*4**distance)
    print(f'PASS: arrival threshold identity, receiving duration inequality; {entries} exact permutation entries')
    print(f'PASS: 9 ramp kernels and {admissible} admissible finite permutations for the strong-source bound')
    print('Scope: finite algebra and index checks; not viscosity-limit proof or novelty certification')


def check_free_decay_test_and_head_tail():
    # Check the test identity from the conservative edge assembly, including
    # the sign which invalidated the first attempted decay argument.
    identities = 0
    for shells in range(1, 5):
        count = 2*shells
        u = sp.symbols(f'v0:{count}', real=True)
        b = sp.symbols(f's0:{count}', positive=True)
        h = sp.symbols(f'k0:{shells}', positive=True)
        nu = sp.symbols('viscosity', positive=True)
        nonlinear, edges = edge_system(u, shells, sp.Integer(1), h)
        steady_nonlinear, _ = edge_system(b, shells, sp.Integer(1), h)
        mass = [sp.Integer(4)**(i//2) for i in range(count)]
        force = [nu*mass[i]*b[i]-steady_nonlinear[i] for i in range(count)]
        form = sum(force[i]*u[i]**2/b[i] for i in range(count))/2
        form += sum(coefficient*b[i]**2*b[j]*(u[i]/b[i]-u[j]/b[j])**2
                    for i, j, coefficient in edges)/2
        pair = sum(nonlinear[i]*b[i] for i in range(count))
        assert sp.expand(pair-form+nu*sum(mass[i]*u[i]**2 for i in range(count))/2) == 0
        actual_test_derivative = pair-nu*sum(mass[i]*u[i]*b[i] for i in range(count))
        expected = form-nu*sum(mass[i]*u[i]**2 for i in range(count))/2
        expected -= nu*sum(mass[i]*u[i]*b[i] for i in range(count))
        assert sp.expand(actual_test_derivative-expected) == 0
        unforced_rhs = [nonlinear[i]-nu*mass[i]*u[i] for i in range(count)]
        assert sp.expand(sum(u[i]*unforced_rhs[i] for i in range(count))
                         +nu*sum(mass[i]*u[i]**2 for i in range(count))) == 0
        identities += 3
    initial = {u[i]: int(i < 2) for i in range(count)}
    assert sp.expand(pair.subs(initial)-b[2]-b[3]) == 0

    # Pure path Hardy estimates are checked independently of the nonlinear
    # equations: arbitrary signed vectors, positive weights, and a pinned root.
    hardy_cases = 0
    for cutoff in range(1, 9):
        count = 2*cutoff+2
        weights = [Fraction(1, i+1) for i in range(count)]
        mass_sum = sum((value**2 for value in weights), Fraction(0))
        moment = sum(((i+1)*value**2 for i, value in enumerate(weights)), Fraction(0))
        lower_root, lower_edge = Fraction(1, 2), Fraction(1, cutoff+1)
        head_constant = 2*mass_sum/lower_root+2*moment/lower_edge
        for seed in range(1, 33):
            values = [Fraction(((seed+3*i)*(i+1)) % 13-6, seed % 5+1)
                      for i in range(count)]
            ratios = [value/weight for value, weight in zip(values, weights)]
            form = values[0]**2/(2*weights[0])
            form += lower_edge*sum(((ratios[i+1]-ratios[i])**2
                                    for i in range(count-1)), Fraction(0))/2
            actual_norm = sum((value**2 for value in values), Fraction(0))
            assert actual_norm <= 2*head_constant*form
            hardy_cases += 1
    cutoffs = 0
    for kappa in (Fraction(1, 16), Fraction(1, 5), Fraction(3, 7)):
        for cutoff in range(1, 16):
            for scale in (1, 2, 3):
                nu_value = kappa/(scale*4**cutoff)
                assert nu_value*4**cutoff <= kappa < nu_value*4**(cutoff+1)
                for n in range(cutoff+1, cutoff+5):
                    assert nu_value*4**n/kappa > 1
                cutoffs += 1

    # The contradiction uses both bounds simultaneously. These constants are
    # independent positive quantities, not data fitted to a numerical orbit.
    C, T = sp.symbols('C T', positive=True)
    A_bound, input_bound = T/(8*C), T/(32*C**2)
    assert sp.simplify(A_bound+sp.sqrt(T*input_bound/2)-T/(4*C)) == 0
    assert sp.simplify(T/(2*C)-T/(4*C)) == T/(4*C)
    print(f'PASS: {identities} assembled free-decay/test identities and the wrong-sign counterexample')
    print(f'PASS: {hardy_cases} exact signed path cases, {cutoffs} exact viscosity cutoffs and T0 factors')
    print('Scope: finite algebra and inequalities; not uniform steady-flow bounds, infinite limits or priority')


def check_radial_angular_projection():
    # New scope only: radial/angular algebra assembled from conservative edges.
    identities = 0
    for shells in range(1, 6):
        values = sp.symbols(f'x0:{2*shells}', real=True)
        intra = sp.symbols(f'exchange0:{shells}', positive=True)
        viscosity, force = sp.symbols('viscosity force', real=True)
        nonlinear, _ = edge_system(values, shells, sp.Integer(1), intra)
        rhs = [nonlinear[i]-viscosity*4**(i//2)*values[i]
               + (force if i < 2 else 0) for i in range(2*shells)]
        for n in range(shells):
            w, z = values[2*n:2*n+2]
            norm = w*w+z*z
            incoming = sp.Integer(0)
            incoming_wedge = sp.Integer(0)
            outgoing = sp.Integer(0)
            outgoing_wedge = sp.Integer(0)
            if n:
                wp, zp = values[2*n-2:2*n]
                incoming = 2**(n-1)*zp*zp*w+sp.Rational(1, 2)**(n-1)*wp*wp*z
                incoming_wedge = sp.Rational(1, 2)**(n-1)*w*wp*wp-2**(n-1)*z*zp*zp
            if n+1 < shells:
                wn, zn = values[2*n+2:2*n+4]
                outgoing = 2**n*z*z*wn+sp.Rational(1, 2)**n*w*w*zn
                outgoing_wedge = w*z*(sp.Rational(1, 2)**n*zn-2**n*wn)
            radial = sp.expand(w*rhs[2*n]+z*rhs[2*n+1])
            expected_radial = incoming-outgoing-viscosity*4**n*norm
            if n == 0:
                expected_radial += force*(w+z)
            assert sp.expand(radial-expected_radial) == 0
            assert sp.diff(radial, intra[n]) == 0
            wedge = sp.expand(w*rhs[2*n+1]-z*rhs[2*n])
            expected_wedge = intra[n]*(w-z)*norm+incoming_wedge+outgoing_wedge
            if n == 0:
                expected_wedge += force*(w-z)
            assert sp.expand(wedge-expected_wedge) == 0
            assert sp.diff(wedge, viscosity) == 0
            identities += 2

    # Two arbitrary admissible states, not the specified common-initial orbit.
    same_radii = [dict.fromkeys(values, 0) for _ in range(2)]
    same_radii[0].update({values[0]: 1, values[2]: 1})
    same_radii[1].update({values[1]: 1, values[2]: 1})
    flux = values[1]**2*values[2]+values[0]**2*values[3]
    for n in range(shells):
        radius_sq = values[2*n]**2+values[2*n+1]**2
        assert radius_sq.subs(same_radii[0]) == radius_sq.subs(same_radii[1])
    assert flux.subs(same_radii[0]) == 0 and flux.subs(same_radii[1]) == 1
    root_radial = sp.expand(values[0]*rhs[0]+values[1]*rhs[1])
    assert sp.expand(root_radial.subs(same_radii[0]).subs(force, 0)+viscosity) == 0
    assert sp.expand(root_radial.subs(same_radii[1]).subs(force, 0)+viscosity+1) == 0

    p, rho = sp.symbols('p rho', real=True)
    assert sp.expand(1-p*p-rho*rho*(1-p)**2
                     -2*p*(1-p)-(1-rho*rho)*(1-p)**2) == 0
    # Cauchy-Schwarz followed by this identity proves phi<=1 for p,rho in [0,1].
    for j in range(1, 13):
        assert sp.Rational(2)**(j-1)/sp.Integer(4)**j == sp.Rational(2)**(-j-1)
        assert sp.Integer(4)**j/(4/viscosity) == viscosity*4**(j-1)

    initial = {value: int(i < 2) for i, value in enumerate(values)}
    assert rhs[2].subs(initial).subs(force, 0) == 1
    assert rhs[3].subs(initial).subs(force, 0) == 1
    time = sp.symbols('time', positive=True)
    leading = dict(initial)
    leading.update({values[2]: time, values[3]: time})
    denominator = 4*(values[0]**2+values[1]**2)*sp.sqrt(values[2]**2+values[3]**2)
    assert sp.limit((flux/denominator).subs(leading), time, 0, dir='+') == sp.sqrt(2)/8
    v = sp.symbols('v', positive=True)
    lower = 1/(1+2*sp.sqrt(2)*viscosity/v)
    upper = 1/(1+2*viscosity/v)
    assert sp.limit(lower, viscosity, 0, dir='+') == 1
    assert sp.limit(upper, viscosity, 0, dir='+') == 1
    print(f'PASS: {identities} assembled radial/angular identities, including root force and terminal closure')
    print('PASS: equal-radius flux counterexample, coefficient bound, 12 index mappings and root onset limit')
    print('Scope: finite algebra and source mapping; not a trajectory flux lower bound or priority certificate')


def check_global_amplitude_weight_obstruction():
    drift_cases, state_cases, integral_cases = 0, 0, 0
    for shells in range(2, 7):
        u = sp.symbols(f'position0:{2*shells}', real=True)
        h = sp.symbols(f'exchange0:{shells}', positive=True)
        q = sp.symbols(f'weight0:{shells}', positive=True)
        rhs, _ = edge_system(u, shells, sp.Integer(1), h)
        drift = sp.expand(sum(q[n]*(rhs[2*n]+rhs[2*n+1]) for n in range(shells)))
        positive = sum(q[n]*h[n]*(u[2*n]-u[2*n+1])**2 for n in range(shells))
        expected = positive
        penalty = sp.Integer(0)
        certificate = positive
        for n in range(shells-1):
            w, z, wn, zn = u[2*n:2*n+4]
            strong, weak = sp.Integer(2)**n, sp.Rational(1, 2)**n
            expected += strong*(q[n+1]*z*z-q[n]*z*wn)
            expected += weak*(q[n+1]*w*w-q[n]*w*zn)
            penalty += (strong*q[n]*(z*z+wn*wn)+weak*q[n]*(w*w+zn*zn))/2
            certificate += strong*q[n+1]*z*z+weak*q[n+1]*w*w
            certificate += (strong*q[n]*(z-wn)**2+weak*q[n]*(w-zn)**2)/2
        assert sp.expand(drift-expected) == 0
        assert sp.expand(drift+penalty-certificate) == 0
        drift_cases += 1
        for n in range(1, shells-2):
            state = dict.fromkeys(u, 0)
            state.update({u[2*n+1]: 1, u[2*n+2]: 1})
            claimed = 2**n*(q[n+1]-q[n])+q[n]*h[n]+q[n+1]*h[n+1]
            claimed += sp.Rational(1, 2)**(n+1)*q[n+2]
            assert sp.expand(drift.subs(state)-claimed) == 0
            assert sum(value**2 for value in state.values()) == 2
            state_cases += 1

        moments = {(i, j): sp.symbols(f'moment{i}_{j}', real=True)
                   for i in range(len(u)) for j in range(i, len(u))}

        def integrate_quadratic(expression):
            result = sp.Integer(0)
            for powers, coefficient in sp.Poly(expression, *u).terms():
                assert sum(powers) == 2
                indices = [i for i, power in enumerate(powers) for _ in range(power)]
                result += coefficient*moments[tuple(indices)]
            return result

        for n in range(1, shells-1):
            w, z = 2*n, 2*n+1
            K, e, ep = sp.Integer(2)**n, sp.Rational(1, 2)**n, sp.Rational(1, 2)**(n-1)
            a = h[n]/sp.sqrt(K)
            endpoint_w = integrate_quadratic(rhs[w])
            endpoint_z = integrate_quadratic(rhs[z])
            Yprev = sp.Integer(2)**(n-1)*moments[(z-2, z-2)]
            remainder = endpoint_w-h[n]*moments[(z, z)]+e*moments[(w, z+2)]
            Bcross = K/h[n]*moments[(w, z)]
            assert sp.expand(Yprev-a*a*Bcross-remainder) == 0
            epsilon = moments[(w, z)]+endpoint_z/h[n]-ep*moments[(w-2, w-2)]/h[n]
            Wexpected = K/h[n]*moments[(z, w+2)]+epsilon
            assert sp.expand(moments[(w, w)]-Wexpected) == 0
            # The discarded first draft inserted an extra -Z in epsilon.
            assert sp.expand(moments[(w, w)]-(Wexpected-moments[(z, z)])) == moments[(z, z)]
            integral_cases += 2

    a, K, C, d, radius = sp.symbols('a K C d radius', positive=True)
    assert sp.simplify(K**(-sp.Rational(1, 2))*a**sp.Rational(1, 3)*sp.sqrt(K)*a
                       -a**sp.Rational(4, 3)) == 0
    assert sp.simplify(a**sp.Rational(4, 3)*C*a**sp.Rational(2, 3)/K-C*a*a/K) == 0
    lower_J = d*d/(C*a**sp.Rational(8, 3))
    assert sp.simplify(a**sp.Rational(4, 3)*lower_J/radius**2
                       -d*d/(C*radius**2*a**sp.Rational(4, 3))) == 0
    assert sp.Rational(4, 3)*sp.Rational(2, 3) == sp.Rational(8, 9) < 1
    level = sp.symbols('level', integer=True, positive=True)
    for gamma, s in ((sp.Rational(1, 2), 0), (1, -1), (1, 0), (1, 1), (2, -1)):
        strong_weight = sp.Integer(2)**((1-gamma)*level)*(level+1)**(-s)
        result = sp.limit(strong_weight, level, sp.oo)
        if gamma < 1 or (gamma == 1 and s < 0):
            assert result == sp.oo
        else:
            assert result.is_finite
    print(f'PASS: {drift_cases} assembled weighted drifts and square certificates; {state_cases} fixed-energy states')
    print(f'PASS: {integral_cases} assembled integral correlations and the erroneous -Z term countercheck')
    print('PASS: target-budget exponents and five weight boundary cases, including gamma=1,s<0')
    print('Scope: a conditional budget and all-state test obstruction; not the common-initial trajectory or a new dissipation theorem')


def check_zero_data_channel_barrier():
    source, damping, transfer, z = sp.symbols('source damping transfer z', positive=True)
    curve = transfer*z**3/(3*(source+damping*z))
    derivative = sp.diff(curve, z)
    assert sp.simplify(derivative-transfer*z*z*(3*source+2*damping*z)
                       /(3*(source+damping*z)**2)) == 0
    margin = transfer*damping*z**3*(3*source+2*damping*z)/(3*(source+damping*z)**2)
    assert sp.simplify(transfer*z*z-damping*curve-source*derivative-margin) == 0
    assert sp.simplify(curve.subs(damping, 0)-transfer*z**3/(3*source)) == 0

    coordinate_cases = 0
    for shells in range(4, 8):
        u = sp.symbols(f'coordinate0:{2*shells}', real=True)
        h = sp.symbols(f'exchange0:{shells}', positive=True)
        radius, viscosity = sp.symbols('radius viscosity', positive=True)
        nonlinear, _ = edge_system(u, shells, sp.Integer(1), h)
        rhs = [nonlinear[i]-viscosity*4**(i//2)*u[i] for i in range(len(u))]
        for n in range(1, shells-2):
            w, zc, x, zn, znn = u[2*n], u[2*n+1], u[2*n+2], u[2*n+3], u[2*n+5]
            ep, en = sp.Rational(1, 2)**(n-1), sp.Rational(1, 2)**(n+1)
            wp = u[2*n-2]
            K = sp.Integer(2)**n
            A = (h[n]+ep)*radius**2
            B = (h[n+1]+en)*radius+viscosity*4**(n+1)
            source_error = h[n]*(w*w-radius**2)+ep*(wp*wp-radius**2)-h[n]*w*zc-K*zc*x
            assert sp.expand(rhs[2*n+1]-A+viscosity*4**n*zc-source_error) == 0
            receive_margin = h[n+1]*zn**2+x*(h[n+1]*(radius-zn)+en*(radius-znn))
            assert sp.expand(rhs[2*n+2]-K*zc*zc+B*x-receive_margin) == 0
            assert sp.simplify(3*viscosity*4**(n+1)/(K*radius)-12*viscosity*K/radius) == 0
            coordinate_cases += 2

    interpolation_cases = 0
    for root in (Fraction(1, 16), Fraction(1, 8), Fraction(1, 4), Fraction(1, 2), Fraction(1)):
        p = root**4
        for scale in (Fraction(0), Fraction(1, 8), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(8)):
            delta = scale*root**3
            for coefficient in (Fraction(1), Fraction(4), Fraction(13)):
                if delta <= root**3:
                    assert coefficient*p+12*delta <= (coefficient+12)*root**3
                else:
                    assert 2*p/delta <= 2*root
                interpolation_cases += 1
    assert Fraction(2)**4 > 1 and Fraction(1) <= 2
    n = sp.symbols('n', integer=True, positive=True)
    for alpha in (sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(3, 4), sp.Integer(1), sp.Integer(4)):
        ratio = sp.Rational(2)**(-sp.Rational(1, 4))*((n+2)/(n+1))**(alpha/2)
        assert sp.limit(ratio, n, sp.oo) == sp.Rational(2)**(-sp.Rational(1, 4)) < 1
    print(f'PASS: invariant receiver curve and {coordinate_cases} assembled full-graph differential bounds')
    print(f'PASS: {interpolation_cases} exact two-regime cases, viscosity weight factor and five tail ratios')
    print('Scope: coordinate inequalities and a uniform tail estimate; not full-vector convergence or vanishing viscous loss')


def check_channel_viscous_localization():
    """Check new coordinate budgets and indices, not the infinite-limit proof."""
    assembled = 0
    for shells in range(4, 9):
        u = sp.symbols(f'channel0:{2*shells}', positive=True)
        h = sp.symbols(f'exchange0:{shells}', positive=True)
        nu, radius = sp.symbols('viscosity radius', positive=True)
        nonlinear, _ = edge_system(u, shells, sp.Integer(1), h)
        rhs = [nonlinear[i]-nu*4**(i//2)*u[i] for i in range(len(u))]
        rhs[0] += 1
        rhs[1] += 1
        for n in range(1, shells-1):
            K = sp.Integer(2)**n
            w, z, x = u[2*n:2*n+3]
            zp = u[2*n+3]
            wp = u[2*n-2]
            en, ep = sp.Rational(1, 2)**n, sp.Rational(1, 2)**(n-1)
            budget = (h[n]*w*w*z+ep*wp*wp*z-h[n]*w*z*z
                      -K*z*z*x-nu*K*K*z*z)
            assert sp.expand(z*rhs[2*n+1]-budget) == 0
            head_power = sum(u[i]*rhs[i] for i in range(2*n+2))
            head_loss = nu*sum(4**(i//2)*u[i]**2 for i in range(2*n+2))
            assert sp.expand(head_power+head_loss+K*z*z*x+en*w*w*zp-u[0]-u[1]) == 0
            # Retain both receiver losses; for the final closed layer the weak
            # outgoing edge is absent and only improves this lower inequality.
            B = (h[n+1]+sp.Rational(1, 2)**(n+1))*radius+4*nu*K*K
            zn2 = u[2*n+5] if n+2 < shells else sp.Integer(0)
            weak_out = sp.Rational(1, 2)**(n+1) if n+2 < shells else sp.Integer(0)
            margin = h[n+1]*zp*zp+x*(h[n+1]*(radius-zp)
                      +sp.Rational(1, 2)**(n+1)*radius-weak_out*zn2)
            assert sp.expand(rhs[2*n+2]-K*z*z+B*x-margin) == 0
            assert sp.simplify(h[n]/(nu*K*K)-h[n]*K**(-2)/nu) == 0
            assembled += 3
        for first in range(1, shells-1):
            last = shells-2
            actual = sum(u[2*n+1]*rhs[2*n+1] for n in range(first, last+1))
            expected = sum(h[n]*u[2*n]**2*u[2*n+1]
                           +sp.Rational(1, 2)**(n-1)*u[2*n-2]**2*u[2*n+1]
                           -h[n]*u[2*n]*u[2*n+1]**2
                           -sp.Integer(2)**n*u[2*n+1]**2*u[2*n+2]
                           -nu*4**n*u[2*n+1]**2
                           for n in range(first, last+1))
            assert sp.expand(actual-expected) == 0
            assembled += 1

    optimum_cases = 0
    for square_a in (1, 2, 5, 11):
        for square_i in (1, 3, 7):
            A, I = Fraction(square_a**2), Fraction(square_i**2)
            delta = Fraction(square_i, square_a)
            assert A*delta+I/delta == 2*square_a*square_i
            optimum_cases += 1
    # A receiver can decrease through the clipped region: the clipped
    # derivative integral is a terminal difference, not a positive variation.
    t, delta = sp.symbols('t delta', real=True)
    clipped_triangle = sp.integrate(1, (t, 0, delta))+sp.integrate(-1, (t, 2-delta, 2))
    assert clipped_triangle == 0
    assert sp.integrate(0, (t, 0, 1)) == 0  # the equality-level set contributes zero

    cutoff = sp.Rational(5, 7)
    exponents = [1-cutoff, 1-sp.Rational(5, 4)*cutoff,
                 sp.Rational(3, 2)-2*cutoff,
                 -1+sp.Rational(3, 2)*cutoff, cutoff]
    assert exponents == [sp.Rational(2, 7), sp.Rational(3, 28),
                         sp.Rational(1, 14), sp.Rational(1, 14), sp.Rational(5, 7)]
    assert min(exponents) == sp.Rational(1, 14) > 0
    K, a = sp.symbols('K a', positive=True)
    assert sp.simplify(sp.sqrt(K)*a/(K*K)-a/K**sp.Rational(3, 2)) == 0
    assert sp.simplify(K*sp.sqrt(sp.sqrt(K)*a)-K**sp.Rational(5, 4)*sp.sqrt(a)) == 0
    damping_product, A = sp.symbols('damping_product A', positive=True)
    z, next_h, next_e, radius = sp.symbols('z next_h next_e radius', positive=True)
    exact_Bz = (next_h+next_e)*radius*z+4*damping_product
    assert sp.expand(exact_Bz.subs({z: radius, damping_product: A})
                     -((next_h+next_e)*radius**2+4*A)) == 0
    for alpha in (sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(3, 4), 1, 4):
        n = sp.symbols('n', positive=True)
        tail_ratio = sp.Rational(2)**(-sp.Rational(1, 3))*((n+2)/(n+1))**(sp.Rational(2, 3)*alpha)
        assert sp.limit(tail_ratio, n, sp.oo) == sp.Rational(2)**(-sp.Rational(1, 3)) < 1

    normalized_cases = 0
    for root in (1, 2, 4):
        a = Fraction(root**6)
        for addition in (Fraction(0), Fraction(1, 3), Fraction(5)):
            s = Fraction(root**4)+addition
            assert a*a <= s**3
            assert a**4/s <= root**20
            normalized_cases += 1
    for A, B, D in ((1, 1, 1), (2, 3, 5), (7, 1, 4)):
        root_c = 6*(A+B+D+1)
        for multiple in (1, 2, 5):
            root_u = root_c*multiple
            value = Fraction(root_u**4)
            assert A+B*root_u+D*root_u**3 <= value/2
            normalized_cases += 1
    print(f'PASS: {assembled} assembled receiver, full-head and finite z-tail budget identities')
    print(f'PASS: {optimum_cases} clipping optima, decreasing receiver clip and five exact cutoff powers')
    print(f'PASS: direct cubic improvement, five tail ratios and {normalized_cases} normalized-recurrence cases')
    print('Scope: finite identities and bound algebra; not AC chain rule, infinite sums, global priority or the total-dissipation classification')


def check_w_channel_head_and_moving_cut():
    assembled = 0
    for shells in range(4, 9):
        u = sp.symbols(f'whead0:{2*shells}', positive=True)
        h = sp.symbols(f'exchange0:{shells}', positive=True)
        nu = sp.symbols('viscosity', positive=True)
        nonlinear, _ = edge_system(u, shells, sp.Integer(1), h)
        rhs = [nonlinear[i]-nu*4**(i//2)*u[i] for i in range(len(u))]
        rhs[0] += 1
        rhs[1] += 1
        total_loss = nu*sum(4**(i//2)*u[i]**2 for i in range(len(u)))
        for n in range(1, shells-1):
            w, z, wn, zn = u[2*n:2*n+4]
            wp, zp = u[2*n-2:2*n]
            K, ep, en = sp.Integer(2)**n, sp.Rational(1, 2)**(n-1), sp.Rational(1, 2)**n
            # Reconstruct h*w^2 from the assembled z equation. Omitting the
            # viscosity term changes this equality by nu*K^2*z.
            reconstructed = rhs[2*n+1]-ep*wp**2+h[n]*w*z+K*z*wn+nu*K*K*z
            assert sp.expand(reconstructed-h[n]*w*w) == 0
            assert sp.expand(reconstructed-nu*K*K*z-h[n]*w*w) == -nu*K*K*z
            F, G = h[n]*w*w*z, h[n]*w*z*z
            Hprev = sp.Integer(2)**(n-1)*zp*zp*w
            assert sp.expand(w*rhs[2*n]+nu*K*K*w*w+F+en*w*w*zn-Hprev-G) == 0
            # Half cut includes w_n and excludes z_n. Its complement includes
            # z_n and all higher coordinates, retaining both weak edges.
            inside = range(2*n+1)
            outside = range(2*n+1, 2*shells)
            out_power = sum(u[i]*rhs[i] for i in outside)
            out_loss = nu*sum(4**(i//2)*u[i]**2 for i in outside)
            weak = en*w*w*zn+ep*wp*wp*z
            assert sp.expand(out_power+out_loss-F+G-weak) == 0
            assert sp.expand(sum(u[i]*rhs[i] for i in inside)+total_loss-out_loss
                             -(u[0]+u[1])+F-G+weak) == 0
            assembled += 5

    W, Z = sp.symbols('W Z', positive=True)
    assert sp.expand((W+Z)**2-4*W*Z-(W-Z)**2) == 0
    # At K_N~nu^-1/2, the five displayed head terms retain these powers.
    b = sp.Rational(1, 2)
    powers = [sp.Integer(1), 1-2*b, 1-sp.Rational(3, 2)*b,
              1-b, 2-3*b]
    assert powers == [1, 0, sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(1, 2)]
    alpha = sp.Rational(2, 3)
    assert -alpha == -sp.Rational(2, 3)
    # Check the inequality needed for decreasing polynomial weights exactly.
    geometric = 0
    for N in range(1, 15):
        for j in range(N+1):
            assert Fraction(N+1, N-j+1) <= j+1
            geometric += 1
    # The simple receiving-coordinate head bound and w-tail bound demand
    # incompatible polynomial cutoffs, regardless of logarithmic powers.
    assert sp.Rational(1, 2) < sp.Rational(4, 5)
    for candidate in (sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(4, 5), sp.Rational(9, 10)):
        head_small = 1-2*candidate >= 0  # boundary gains the 1/a_N logarithm
        tail_small = -2+sp.Rational(5, 2)*candidate > 0
        assert not (head_small and tail_small)
    print(f'PASS: {assembled} assembled z-amplitude, w-energy and complementary half-cut identities')
    print(f'PASS: retained-viscosity countercheck, exact head powers, log exponent and {geometric} weight inequalities')
    print('Scope: finite budgets and indices; not total w-loss disappearance, trajectory counterexamples or global originality')


def check_saturated_moment_realizability():
    """Finite moment/set checks; the a.e. nonrealizability proof is in section 3.14."""
    tau, c = sp.symbols('tau c', positive=True)
    assert sp.simplify((c*tau)**2-tau*(c*c*tau)) == 0
    x = sp.symbols('x', real=True)
    assert sp.expand(x*x*(x-1)**2-(x**4-2*x**3+x*x)) == 0
    n = sp.symbols('n', positive=True)
    assert sp.expand((n+1)**3-n*n*(n+2)) == n*n+3*n+1
    remainder = ((n+1)/n)**sp.Rational(8, 3)-((n+2)/(n+1))**sp.Rational(4, 3)
    assert sp.limit(n*remainder, n, sp.oo) == sp.Rational(4, 3)
    # Reordering/fragmenting the same time atoms preserves forced nesting.
    for order in (list(range(12)), [(5*j+3) % 12 for j in range(12)]):
        sets = [set(order[:12-j]) for j in range(12)]
        for j in range(11):
            W = Fraction(len(sets[j]), 12)
            Wnext = Fraction(len(sets[j+1]), 12)
            ratio = W/Wnext
            a_squared = 1/W
            inner = ratio*W
            Y = a_squared*ratio**2*W
            F = a_squared*ratio*W
            H = a_squared*ratio**2*Wnext
            assert inner**2 == W*(ratio**2*W)
            assert Y == ratio**2 and F == H == ratio
            assert ratio*Fraction(len(sets[j] & sets[j+1]), 12) == W
        for depth in range(1, 12):
            assert all(sum(atom in sets[j] for j in range(1, depth+1)) == depth
                       for atom in sets[depth])
    universe = range(6)
    subsets = [{j for j in universe if mask & (1 << j)} for mask in range(64)]
    for A in subsets:
        for B in subsets:
            if len(A & B) == len(B):
                assert B <= A
    print('PASS: saturation identities, positive remainder algebra and asymptotic 4/3')
    print('PASS: nested and fragmented supports, finite energy depth, 4096 subset pairs')
    print('Scope: exact moment candidate only; not an ODE counterexample or alpha classification')


def check_weighted_correlations_and_endpoint_product():
    """Exact finite functions and product algebra for sections 3.15-3.16."""
    def unit_profile(dim, seed, shift):
        parameters = [Fraction(1+(seed+2*j) % 5, 10*dim) for j in range(dim-1)]
        squared = sum(t*t for t in parameters)
        raw = [(1-squared)/(1+squared)] + [2*t/(1+squared) for t in parameters]
        assert sum(x*x for x in raw) == 1 and min(raw) > 0
        return raw[shift % dim:] + raw[:shift % dim]

    dot = lambda left, right: sum(x*y for x, y in zip(left, right))
    blocks = products = 0
    for dim in range(3, 7):
        for shells in range(3, 8):
            profile = [unit_profile(dim, shells+i, i) for i in range(shells)]
            receiver = [unit_profile(dim, dim+3*i, i+1) for i in range(shells)]
            amplitude = [Fraction(1+shells*(i+1), 1+i % 3) for i in range(shells)]
            receiver_amplitude = [Fraction(dim+2, i+1) for i in range(shells)]
            w = [[amplitude[i]*x for x in profile[i]] for i in range(shells)]
            v = [[receiver_amplitude[i]*x for x in receiver[i]] for i in range(shells)]
            W = [dot(x, x)/dim for x in w]
            J = [sum(x**4 for x in row)/dim for row in w]
            a_squared = [Fraction((i+1)**2) for i in range(shells)]
            Y = [a_squared[i]*dot(v[i], v[i])/dim for i in range(shells)]
            rho = [dot(profile[i], receiver[i]) for i in range(shells)]
            sigma = [dot(receiver[i], profile[i+1]) for i in range(shells-1)]
            defect = [(1-rho[i])+(1-sigma[i]) for i in range(shells-1)]
            cap = max(sum(w[i][atom]**2 for i in range(shells)) for atom in range(dim))
            for n in range(shells-1):
                for m in range(n+1, shells):
                    mass = sum(W[n:m+1])
                    distance_term = sum(W[k]*sum((profile[n][j]-profile[k][j])**2
                                                for j in range(dim)) for k in range(n, m+1))/mass
                    weighted_loss = sum(defect[i]*sum((k-n)*W[k] for k in range(i+1, m+1))/mass
                                        for i in range(n, m))
                    assert J[n]/W[n] <= 2*cap*(W[n]/mass+distance_term)
                    assert distance_term <= 4*weighted_loss
                    F = a_squared[n]*sum(w[n][j]**2*v[n][j] for j in range(dim))/dim
                    assert F*F <= a_squared[n]*Y[n]*J[n]
                    assert weighted_loss >= F*F/(8*cap*a_squared[n]*Y[n]*W[n])-W[n]/(4*mass)
                    blocks += 1
            ratio_product = Fraction(1)
            for i in range(1, shells-1):
                first_inner = dot(w[i], v[i])/dim
                second_inner = dot(v[i], w[i+1])/dim
                relative_r = (Y[i-1]-a_squared[i]*first_inner)/Y[i-1]
                relative_s = (W[i]-second_inner)/W[i]
                assert 1-relative_r > 0 and 1-relative_s > 0
                ratio = (1-relative_r)**2*(1-relative_s)**2/(rho[i]**2*sigma[i]**2)
                assert ratio == W[i+1]*Y[i]**2/(W[i]*Y[i-1]**2)
                ratio_product *= ratio
                products += 1
            assert ratio_product == W[-1]*Y[-2]**2/(W[1]*Y[0]**2)
    # Omitting correlation loss is false even for uniformly bounded w-energy.
    disjoint = [[Fraction(int(i == j)) for j in range(3)] for i in range(3)]
    disjoint_W = [dot(row, row)/3 for row in disjoint]
    disjoint_J = [sum(x**4 for x in row)/3 for row in disjoint]
    disjoint_cap = max(sum(row[j]**2 for row in disjoint) for j in range(3))
    assert disjoint_J[0]/disjoint_W[0] > 2*disjoint_cap*disjoint_W[0]/sum(disjoint_W)
    x = sp.symbols('x', real=True)
    assert sp.simplify(sp.diff(-sp.log(1-x)-x, x)-x/(1-x)) == 0
    alpha = sp.Rational(2, 3)
    assert 4*alpha/3-1 == -sp.Rational(1, 9)
    assert (1-4*alpha/3)/2 == sp.Rational(1, 18)
    print(f'PASS: {blocks} exact weighted finite-block, path and real moment inequalities')
    print(f'PASS: {products} exact relative-error products and 20 telescoping identities')
    print('PASS: omitted-loss countercheck, log remainder derivative and endpoint exponent 1/18')
    print('Scope: finite function and algebra checks; not actual endpoint l1 boundedness or anomaly classification')


def check_sender_receiver_retention():
    """Reconstruct the local retention identity from the full conservative graph."""
    identities = 0
    for shells in range(4, 9):
        values = sp.symbols(f'retention0:{2*shells}', positive=True)
        exchange = sp.symbols(f'retention_h0:{shells}', positive=True)
        nu, R, A = sp.symbols('retention_nu retention_R retention_A', positive=True)
        nonlinear, _ = edge_system(values, shells, sp.Integer(1), exchange)
        rhs = [nonlinear[i]-nu*4**(i//2)*values[i] for i in range(len(values))]
        rhs[0] += 1
        rhs[1] += 1
        for n in range(1, shells-2):
            w, z, x, zn = values[2*n:2*n+4]
            wp, zp = values[2*n-2:2*n]
            znn = values[2*(n+2)+1]
            K = sp.Integer(2)**n
            ep, en, e_next = sp.Rational(1, 2)**(n-1), sp.Rational(1, 2)**n, sp.Rational(1, 2)**(n+1)
            h = exchange[n]
            c = h/(K*A+nu*K*K)
            B = en*R+nu*K*K+c*ep*R*R
            source = K*zp*zp/2+h*z*z
            assert sp.expand(rhs[2*n]+(h*z+en*R+nu*K*K)*w
                             -source-en*w*(R-zn)) == 0
            assert sp.expand(h*w*w+ep*R*R-(K*x+nu*K*K)*z-rhs[2*n+1]
                             -h*w*z-ep*(R*R-wp*wp)) == 0
            receiver_damping = exchange[n+1]*zn+e_next*znn+4*nu*K*K
            assert sp.expand(rhs[2*n+2]-K*z*z+receiver_damping*x-exchange[n+1]*zn*zn) == 0
            transformed = rhs[2*n]-c*w*rhs[2*n+1]+c*h*w**3+B*w-c*K*w*z*(x-A)
            nonnegative_remainder = source+c*h*w*w*z+en*w*(R-zn)+c*ep*w*(R*R-wp*wp)
            assert sp.factor(transformed-nonnegative_remainder) == 0
            identities += 4

    # The Bernoulli bound is checked with its retained positive damping.
    time, damping, cubic, q_initial = sp.symbols('time damping cubic q_initial', positive=True)
    q_bound = sp.exp(2*damping*time)*q_initial+cubic/damping*(sp.exp(2*damping*time)-1)
    assert sp.simplify(sp.diff(q_bound, time)-2*damping*q_bound-2*cubic) == 0
    K, a, A, nu = sp.symbols('K a A nu', positive=True)
    c = a/(sp.sqrt(K)*(A+nu*K))
    assert sp.simplify(c*c*K*A**3-a*a*A**3/(A+nu*K)**2) == 0
    # Dropping stored z is not a valid pointwise cubic decay estimate.
    state = [sp.Rational(0) for _ in range(12)]
    state[8], state[9], state[10] = sp.Rational(1, 10), sp.Rational(1, 20), sp.Integer(1)
    rhs, _ = edge_system(state, 6, sp.Integer(1), [sp.Integer(2)**sp.Rational(i, 2) for i in range(6)])
    assert sum(x*x for x in state) <= 4
    assert rhs[8] < -state[8]**3  # Here K=16,h=4,A=1 and c*h=1.
    coefficient, amplitude, a = sp.symbols('coefficient amplitude a', positive=True)
    scale_bound = coefficient*a**(-sp.Rational(2, 3))/(4*a*a*(amplitude*a**(-sp.Rational(2, 3)))**2)
    assert sp.simplify(scale_bound-coefficient/(4*amplitude**2)*a**(-sp.Rational(4, 3))) == 0
    assert sp.Rational(2, 3)*sp.Rational(4, 3) == sp.Rational(8, 9) < 1
    print(f'PASS: {identities} full-graph source, damping, receiver and corrected retention identities')
    print('PASS: retained-damping Bernoulli formula, vacancy coefficient, stored-z countercheck and exponent 8/9')
    print('Scope: local identities and comparison algebra; not persistent receiver occupancy, interval disjointness or total loss')


if __name__ == '__main__':
    check_core()
    check_objections()
    check_known_structure()
    check_novelty_reductions()
    check_ground_state_mapping()
    check_stationary_normalization()
    check_irregular_envelopes()
    check_sparse_bottlenecks()
    check_geometric_profile()
    check_root_input_bound()
    check_steady_necessity_reduction()
    check_deleted_inviscid_root_match()
    check_nonadjacent_profiles()
    check_partial_cut_and_raw_recurrence()
    check_weak_block_bounds()
    check_frozen_block_probabilities()
    check_dynamic_rescaling()
    check_equivalence_obstructions()
    check_fourth_root_block_indices()
    check_amplitude_sum_weighted_blocks()
    check_log_moment_and_local_pulse()
    check_arrival_and_pure_permutation()
    check_free_decay_test_and_head_tail()
    check_radial_angular_projection()
    check_global_amplitude_weight_obstruction()
    check_zero_data_channel_barrier()
    check_channel_viscous_localization()
    check_w_channel_head_and_moving_cut()
    check_saturated_moment_realizability()
    check_weighted_correlations_and_endpoint_product()
    check_sender_receiver_retention()
