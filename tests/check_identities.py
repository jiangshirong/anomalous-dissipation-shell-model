"""Exact finite algebra checks. NOT a proof of a viscosity-limit anomaly."""

import sympy as sp


def verify_closed_truncations():
    checked = 0
    for last in range(9):
        a = sp.symbols(f"a0:{last + 1}", real=True)
        p = sp.symbols(f"p0:{last + 1}", real=True)
        q = sp.symbols(f"q0:{last + 1}", real=True)
        nu, force = sp.symbols("nu F", real=True)
        rhs = [-nu * 4**i * a[i] + (force if i == 0 else 0)
               for i in range(last + 1)]
        for i in range(last + 1):
            for step, coefficient in ((1, p[i]), (2, q[i])):
                j = i + step
                if j <= last:
                    rhs[i] -= coefficient * a[i] * a[j]
                    rhs[j] += coefficient * a[i] ** 2
        energy_rate = sum(a[i] * rhs[i] for i in range(last + 1))
        target = force * a[0] - nu * sum(4**i * a[i] ** 2
                                        for i in range(last + 1))
        assert sp.expand(energy_rate - target) == 0, last
        checked += 1
    return checked


def verify_internal_cuts():
    checked = 0
    for cut in range(9):
        a = sp.symbols(f"a0:{cut + 3}", real=True)
        p = sp.symbols(f"p0:{cut + 3}", real=True)
        q = sp.symbols(f"q0:{cut + 3}", real=True)
        rate = sp.Integer(0)
        for n in range(cut + 1):
            nonlinear = -a[n] * (p[n] * a[n + 1] + q[n] * a[n + 2])
            if n >= 1:
                nonlinear += p[n - 1] * a[n - 1] ** 2
            if n >= 2:
                nonlinear += q[n - 2] * a[n - 2] ** 2
            rate += a[n] * nonlinear
        flux = p[cut] * a[cut] ** 2 * a[cut + 1]
        flux += q[cut] * a[cut] ** 2 * a[cut + 2]
        if cut >= 1:
            flux += q[cut - 1] * a[cut - 1] ** 2 * a[cut + 1]
        assert sp.expand(rate + flux) == 0, cut
        checked += 1
    return checked


def verify_rescaling():
    i = sp.symbols("i", integer=True, nonnegative=True)
    I, nu, A, B, s = sp.symbols("I nu A B s", positive=True)
    ai = s * I ** sp.Rational(1, 3) * 2 ** (-sp.Rational(2, 3) * i)
    for step in (3, 4):
        k = i + step
        raw = A * 4**i * ai**4 / I - B * nu
        scaled = 2 ** (sp.Rational(2, 3) * k) * raw / I ** sp.Rational(1, 3)
        t = nu * I ** (-sp.Rational(1, 3)) * 2 ** (sp.Rational(2, 3) * i)
        target = 2 ** (sp.Rational(2, 3) * step) * (A * s**4 - B * t)
        assert sp.simplify(sp.expand_power_base(scaled - target)) == 0, step
    return 2


def verify_coverage_constants():
    """Check both strong-edge lengths; the infinite iteration is not tested."""
    i = sp.symbols("i", integer=True, nonnegative=True)
    kappa = sp.symbols("kappa", positive=True)
    for step in (1, 2):
        normalized = sp.powsimp((kappa * 4**i)**2 /
                                (2 * 4**(i + step)) / (kappa**2 * 4**i))
        assert sp.simplify(normalized - sp.Rational(1, 2 * 4**step)) == 0
        assert normalized >= sp.Rational(1, 32)
        for outgoing_step in (1, 2):
            assert 2 <= step + outgoing_step <= 4
    # H=6 with kappa=1 ensures main term >= 4*a_i and >= nu.
    assert sp.Rational(1, 32) * 6**3 >= 4
    assert sp.Rational(1, 32) * 6**4 >= 1
    assert 12**3 < 64**2
    assert 72 * (32 + 16) < 4096
    assert 1 + sp.sqrt(4096) == 65
    return 2


def verify_sparse_cut_support():
    """Enumerate all crossing edges, independently of the three-term formula."""
    checked = 0
    eta = sp.symbols("eta", positive=True)
    defects = {j*j for j in range(2, 8)}
    repairs = defects - {9, 25}  # Both repaired and unrepaired cuts.
    a = sp.symbols("a0:52", real=True)
    for cut in sorted(defects):
        flux = sp.Integer(0)
        for source in range(cut + 1):
            primary = sp.Rational(1, 2**source) if source in defects else 4**source
            bypass = eta * 4**source if source in repairs else 0
            for step, coefficient in ((1, primary), (2, bypass)):
                target = source + step
                if target > cut:
                    flux += coefficient * a[source]**2 * a[target]
        expected = sp.Rational(1, 2**cut) * a[cut]**2 * a[cut + 1]
        if cut in repairs:
            expected += eta * 4**cut * a[cut]**2 * a[cut + 2]
        assert sp.expand(flux - expected) == 0, cut
        checked += 1
    return checked


def verify_difference_energy():
    """Algebra behind the evolution-uniqueness estimate, not Gronwall itself."""
    checked = 0
    for last in range(9):
        a = sp.symbols(f"a0:{last + 1}", real=True)
        b = sp.symbols(f"b0:{last + 1}", real=True)
        p = sp.symbols(f"p0:{last + 1}", real=True)
        q = sp.symbols(f"q0:{last + 1}", real=True)
        d = [a[i] - b[i] for i in range(last + 1)]
        rhs_a = [sp.Integer(0) for _ in a]
        rhs_b = [sp.Integer(0) for _ in b]
        target = sp.Integer(0)
        for i in range(last + 1):
            for step, coefficient in ((1, p[i]), (2, q[i])):
                j = i + step
                if j <= last:
                    rhs_a[i] -= coefficient * a[i] * a[j]
                    rhs_a[j] += coefficient * a[i]**2
                    rhs_b[i] -= coefficient * b[i] * b[j]
                    rhs_b[j] += coefficient * b[i]**2
                    target += coefficient * (b[i] * d[i] * d[j] - b[j] * d[i]**2)
        pairing = sum(d[i] * (rhs_a[i] - rhs_b[i]) for i in range(last + 1))
        assert sp.expand(pairing - target) == 0, last
        checked += 1
    return checked


def verify_stationary_elimination():
    x, z, p, r, q, damping = sp.symbols("x z p r q damping", positive=True)
    for viscous_term in (sp.Integer(0), damping):
        y = p * x**2 / (viscous_term + r * z)
        assert sp.simplify(p * x**2 - (viscous_term + r*z) * y) == 0
        flux = p * x**2 * y + q * x**2 * z
        expected = q * x**2 * z + p**2 * x**4 / (viscous_term + r*z)
        assert sp.simplify(flux - expected) == 0
        # Not a constant effective cubic coefficient on the retained variables.
        ratio_derivative = sp.factor(sp.diff(flux / (x**2 * z), x))
        assert sp.simplify(ratio_derivative -
                           2*p**2*x / (z*(viscous_term + r*z))) == 0
        assert ratio_derivative.is_positive
    return 2


def verify_relative_lyapunov():
    """Exact graph identity with forcing chosen as the equilibrium residual."""
    checked = 0
    for last in range(9):
        b = sp.symbols(f"b0:{last + 1}", positive=True)
        d = sp.symbols(f"d0:{last + 1}", real=True)
        a = [b[i] + d[i] for i in range(last + 1)]
        p = sp.symbols(f"p0:{last + 1}", positive=True)
        q = sp.symbols(f"q0:{last + 1}", nonnegative=True)
        nu = sp.symbols("nu", positive=True)
        nonlinear_a = [sp.Integer(0) for _ in a]
        nonlinear_b = [sp.Integer(0) for _ in b]
        edge_squares = sp.Integer(0)
        for i in range(last + 1):
            for step, c in ((1, p[i]), (2, q[i])):
                j = i + step
                if j <= last:
                    nonlinear_a[i] -= c*a[i]*a[j]
                    nonlinear_a[j] += c*a[i]**2
                    nonlinear_b[i] -= c*b[i]*b[j]
                    nonlinear_b[j] += c*b[i]**2
                    edge_squares += c*b[i]**2*b[j]*(d[i]/b[i]-d[j]/b[j])**2
        forcing = [nu*4**i*b[i] - nonlinear_b[i] for i in range(last + 1)]
        energy_derivative = sum(
            d[i]*(nonlinear_a[i] - nonlinear_b[i] - nu*4**i*d[i])
            for i in range(last + 1))
        target = (-sp.Rational(1, 2)*edge_squares
                  -nu*sp.Rational(1, 2)*sum(4**i*d[i]**2 for i in range(last + 1))
                  -sp.Rational(1, 2)*sum(forcing[i]*d[i]**2/b[i]
                                         for i in range(last + 1)))
        assert sp.expand(energy_derivative-target) == 0, last
        checked += 1
    return checked


def verify_positive_coordinate_changes():
    """Square-input support under weighted permutations, not all conjugacies."""
    checked = 0
    for last in range(6, 11):
        a = sp.symbols(f"a0:{last + 1}", real=True)
        y = sp.symbols(f"y0:{last + 1}", real=True)
        rho = sp.symbols(f"rho0:{last + 1}", positive=True)
        p = sp.symbols(f"p0:{last + 1}", positive=True)
        q4 = sp.symbols("q4", positive=True)
        rhs = [sp.Integer(0) for _ in a]
        for i in range(last):
            rhs[i] -= p[i]*a[i]*a[i+1]
            rhs[i+1] += p[i]*a[i]**2
        rhs[4] -= q4*a[4]*a[6]
        rhs[6] += q4*a[4]**2
        permutations = (list(reversed(range(last + 1))),
                        [(i + 2) % (last + 1) for i in range(last + 1)])
        for sigma in permutations:
            substitution = {a[i]: y[sigma[i]]/rho[i] for i in range(last+1)}
            transformed = [sp.Integer(0) for _ in a]
            for i in range(last+1):
                transformed[sigma[i]] = rho[i]*rhs[i].subs(substitution,
                                                           simultaneous=True)
            for i in range(last+1):
                polynomial = sp.Poly(transformed[sigma[i]], *y)
                squares = {j for j in range(last+1)
                           if polynomial.coeff_monomial(y[j]**2) != 0}
                expected = set() if i == 0 else {sigma[i-1]}
                if i == 6:
                    expected.add(sigma[4])
                    assert sp.simplify(polynomial.coeff_monomial(
                        y[sigma[4]]**2) - rho[6]*q4/rho[4]**2) == 0
                assert squares == expected, (last, sigma, i)
            checked += 1
    return checked


def verify_n1_face_and_support():
    """Read-back checks of PSF arXiv v2 Eq. (233) and Appendix C support."""
    checked = 0
    # At a non-root zero coordinate, the only nonzero amplitude is U_{n+1}=1.
    # All omitted m terms vanish on this finite-support state.
    u = lambda offset: sp.Integer(1 if offset == 1 else 0)
    for lam in (sp.Integer(2), sp.Integer(3), sp.Integer(4)):
        kn = lam**3
        rate = sp.Integer(0)
        for m in range(-1, 6):
            coefficient = (sp.Rational(387, 1000) if m == -1 else
                           sp.Integer(0) if m == 0 else sp.Rational(219, 100))
            rate += kn*coefficient*(lam**(-sp.Rational(5, 2)*m)*u(-m)*u(0)
                                    - lam**(-sp.Rational(3, 2)*m)*u(m)**2)
        expected = -kn*sp.Rational(219, 100)*lam**(-sp.Rational(3, 2))
        assert sp.simplify(rate-expected) == 0
        assert expected.is_negative
        checked += 1
    # Conjugation changes coefficients, not offsets, on a real restriction.
    # Finite cases check transcription; the universal support argument is in prose.
    for m in range(1, 13):
        offsets = {tuple(sorted(pair)) for pair in (
            (-m, -1), (-m+1, 1), (-1, -m), (m-1, m),
            (0, -m), (m, m), (-m+2, 0))}
        assert (-2, -2) not in offsets
        assert (-1, -1) in offsets if m == 1 else (-1, -1) not in offsets
        checked += 1
    return checked


def verify_fixed_window_comparison():
    """Scalar identities only; uniform dynamical norm control is a theorem."""
    F, L, C, T, t = sp.symbols("F L C T t", positive=True)
    initial = sp.symbols("initial", nonnegative=True)
    checked = 0
    for value in (sp.Integer(0), initial):
        g = value*sp.exp(-L*t) + F/L*(1-sp.exp(-L*t))
        assert sp.simplify(sp.diff(g, t) - (F-L*g)) == 0
        assert sp.simplify(g.subs(t, 0)-value) == 0
        input_integral = sp.integrate(F*g, (t, 0, T))
        expected = F**2*T/L + (F*value/L-F**2/L**2)*(1-sp.exp(-L*T))
        assert sp.simplify(input_integral-expected) == 0
        checked += 1
    bound = F**2/L-(F**2/L**2+C**2/2)/T
    threshold = 2/L+L*C**2/F**2
    assert sp.simplify(bound.subs(T, threshold)-F**2/(2*L)) == 0
    assert sp.diff(bound, T).is_positive
    checked += 1
    explicit_threshold = threshold.subs({F: 1, L: 129, C: 128})
    assert explicit_threshold == 2113536+sp.Rational(2, 129)
    assert bound.subs({F: 1, L: 129, C: 128, T: 2113537}) >= sp.Rational(1, 258)
    checked += 1
    return checked


def verify_closed_window_bound():
    """Exact scalar integration and a necessary, not sufficient, cutoff test."""
    A, F, T, nu, t = sp.symbols("A F T nu t", positive=True)
    last = sp.symbols("last", integer=True, nonnegative=True)
    integrated = sp.integrate(nu*4**last*(A+F*t)**2, (t, 0, T))/T
    target = nu*4**last*(A**2+A*F*T+F**2*T**2/3)
    assert sp.simplify(integrated-target) == 0
    assert sp.limit(target, nu, 0, dir="+") == 0
    checked = 1
    # Exact arithmetic on coupled refinements. These concern an upper bound,
    # not computed solutions and not sufficiency of a selected shell count.
    for m in (1, 2, 4, 8):
        viscosity = sp.Rational(1, 16**m)
        assert viscosity*4**m == sp.Rational(1, 4**m)
        assert viscosity*4**(2*m) == 1
        checked += 1
    return checked


def verify_strong_path_obstruction():
    """Check local elimination, all crossing edges, and scalar bound constants."""
    checked = 0
    nu = sp.symbols("nu", nonnegative=True)
    for cut in range(0, 16, 2):
        a = sp.symbols(f"a0:{cut + 3}", nonnegative=True)
        p = lambda i: sp.Integer(4)**i if i % 2 == 0 else sp.Rational(1, 2**i)
        q = lambda i: sp.Integer(4)**i if i % 2 == 0 else sp.Integer(0)
        # Enumerate crossing edges rather than inserting the proposed flux.
        flux = sp.Integer(0)
        for i in range(cut + 1):
            for step, coefficient in ((1, p(i)), (2, q(i))):
                if i + step > cut:
                    flux += coefficient*a[i]**2*a[i + step]
        j = cut + 1
        odd_rhs = p(j - 1)*a[j - 1]**2 - a[j]*(p(j)*a[j + 1]+nu*4**j)
        if j >= 2:
            odd_rhs += q(j - 2)*a[j - 2]**2
        y, z = a[j], a[j + 1]
        eliminated = (nu*4**j+p(j)*z)*y*(y+z)
        # The difference is exactly the stationary residual times y+z.
        assert sp.expand(flux-eliminated-odd_rhs*(y+z)) == 0
        for n in (cut, cut + 1):
            total = p(n)+q(n)+(q(n-1) if n else 0)
            assert total >= sp.Rational(1, 4)*4**n
        assert q(cut) == 4**cut
        checked += 1
    m = sp.symbols("m", integer=True, nonnegative=True)
    delta = sp.Integer(4)**(-m)
    assert sp.simplify(sp.Integer(64)**(-m)*sp.Integer(4)**(2*m)-delta) == 0
    assert sp.simplify(2*sp.Integer(2)**(-(2*m+1))-delta) == 0
    # At the lower endpoint nu=64^(-(m+1)), delta=4*nu^(1/3).
    assert sp.simplify(4*sp.Integer(64)**(-(m+1)/sp.Integer(3))-delta) == 0
    checked += 1
    R = sp.symbols("R", positive=True)
    assert sp.expand(19*R**3-(9*R**2+R**3)-9*R**2*(2*R-1)) == 0
    assert sp.simplify(1/(1+2*R)-1/(4*R)-(2*R-1)/(4*R*(1+2*R))) == 0
    assert 4*19*4 == 304
    checked += 1
    return checked


def verify_weakly_fed_outlets():
    """Finite algebra supporting the separate infinite stationary estimates."""
    checked = 0
    for m in range(8):
        r = sp.Rational(1, 2**(6*m+7))
        s = sp.Rational(1, 2**(2*m+1))
        c = sp.Integer(16)**m
        r_next = sp.Rational(1, 2**(6*m+13))
        assert r/s == sp.Rational(1, 2**(4*m+6))
        assert 16*r/s == sp.Rational(1, 4*16**m)
        assert r*r_next/(s*c) == sp.Rational(1, 2**(14*m+19))
        assert r_next/(16*c) == sp.Rational(1, 2**(10*m+17))
        assert s/(4**(2*m+1)) == sp.Rational(1, 2**(3*(2*m+1)))
        assert r*r_next/(s*c) <= sp.Rational(1, 2**19*16**m)
        assert r_next/(16*c) <= sp.Rational(1, 16**m)
        checked += 1
    r, s, c, c_next, r_next, nu = sp.symbols(
        "r s c c_next r_next nu", positive=True)
    x, x_next, x_after, y, y_next = sp.symbols(
        "x x_next x_after y y_next", positive=True)
    odd_residual = r*x**2-(s*x_next+4*nu*c)*y
    outgoing = c_next*x_after+r_next*y_next+nu*c_next
    even_residual = c*x**2+s*y**2-x_next*outgoing
    gap = r/(s*c)*outgoing-y
    positive_terms = 4*nu*c*y/(s*x_next)+r*y**2/(c*x_next)
    # Zero residuals make the gap nonnegative; no infinite claim is tested.
    assert sp.simplify(gap-positive_terms
                       -odd_residual/(s*x_next)
                       +r*even_residual/(s*c*x_next)) == 0
    checked += 1
    m = sp.symbols("m", integer=True, nonnegative=True)
    I, R = sp.symbols("I R", positive=True)
    thresholds = [
        4*I**sp.Rational(1, 3)*16**(-m/sp.Integer(3)),
        4*(nu*I)**sp.Rational(1, 4)*16**(-m/sp.Integer(4)),
        4*(I*R)**sp.Rational(1, 4)*16**(-m/sp.Integer(2))]
    targets = [4*thresholds[0], 16*nu, 16*16**(-m)*R]
    for threshold, target in zip(thresholds, targets):
        main = 16**m*threshold**4/(16*I)
        assert sp.simplify(main-target) == 0
        checked += 1
    assert sp.Rational(1, 4)+sp.Rational(1, 2**19) < sp.Rational(1, 3)
    assert 1/(9*(1-sp.Rational(1, 16**2))) < sp.Rational(1, 8)
    # First geometric ratio is <=1/4 because its cube is 1/256<=1/64.
    assert sp.Rational(1, 16**2) <= sp.Rational(1, 4**3)
    assert 48*sp.Rational(4, 3) == 64
    assert sp.expand(4*(R**2-(128+192*R+R**2/4+sp.Rational(1, 4)))
                     -(3*R**2-768*R-513)) == 0
    assert 3*257**2-768*257-513 > 0
    assert 6*257-768 > 0
    assert 1/(1+sp.Rational(129, 128)*257) > sp.Rational(1, 261)
    checked += 1
    return checked


def verify_feeding_perturbation_bounds():
    """Finite algebra for the stated band, not a general stability test."""
    checked = 0
    rho, rho_next = sp.symbols("rho rho_next", positive=True)
    for m in range(6):
        c = sp.Integer(16)**m
        s = sp.Rational(1, 2**(2*m+1))
        r = rho*sp.Rational(1, 2**(6*m+7))
        rn = rho_next*sp.Rational(1, 2**(6*m+13))
        K = r/s
        assert sp.simplify((16*K).subs(rho, 2)
                           -sp.Rational(1, 2*16**m)) == 0
        assert sp.simplify((K*rn/c).subs({rho: 2, rho_next: 2})
                           -sp.Rational(1, 2**(14*m+17))) == 0
        assert sp.Rational(1, 2**(14*m+17)) <= sp.Rational(1, 2**17*16**m)
        assert (rn/(16*c)).subs(rho_next, 2) <= sp.Rational(1, 16**m)
        checked += 1
    R = sp.symbols("R", positive=True)
    assert sp.Rational(1, 2)+sp.Rational(1, 2**17) < sp.Rational(9, 16)
    assert sp.Rational(9, 16)**2/(1-sp.Rational(1, 16**2)) == sp.Rational(27, 85)
    # The difference from the chosen loose bound is visibly nonnegative.
    gap = sp.Rational(2, 3)*(R**2+1)-sp.Rational(27, 85)*(R+1)**2
    assert sp.expand(gap-(sp.Rational(27, 85)*(R-1)**2
                          +sp.Rational(8, 255)*(R**2+1))) == 0
    assert sp.expand(3*(R**2-(128+192*R+sp.Rational(2, 3)*(R**2+1)))
                     -(R**2-576*R-386)) == 0
    assert 577**2-576*577-386 > 0
    assert 2*577-576 > 0
    assert 1/(1+sp.Rational(65, 64)*577) > sp.Rational(1, 588)
    checked += 1
    x, y, z, c, nu, e = sp.symbols("x y z c nu e", positive=True)
    odd_residual = rho*c*x**2-(4*nu*c+e*z)*y
    flux = rho*c*x**2*y+c*x**2*z
    eliminated = (4*nu*c+e*z)*y*(y+z/rho)
    assert sp.simplify(flux-eliminated-odd_residual*(y+z/rho)) == 0
    checked += 1
    assert sp.expand(sp.Rational(81, 2)*R**3-(13*R**2+sp.Rational(3, 2)*R**3)
                     -13*R**2*(3*R-1)) == 0
    assert sp.simplify(1/(1+3*R)-1/(6*R)
                       -(3*R-1)/(6*R*(1+3*R))) == 0
    assert 6*4*sp.Rational(81, 2) == 972
    # Returning from tau=2t preserves the physical energy budget.
    h1_sq, b0 = sp.symbols("h1_sq b0", positive=True)
    assert sp.expand(2*(nu*h1_sq/2-b0/2)-(nu*h1_sq-b0)) == 0
    checked += 1
    return checked


def verify_weakfed_stationary_asymptotics():
    """Exact eliminated flux, rescaling and recurrence; no limit proof."""
    checked = 0
    r, s, c, x, z, I = sp.symbols("r s c x z I", positive=True)
    y = r*x**2/(s*z)
    flux = (c*x**2+s*y**2)*z
    assert sp.simplify(flux-(c*x**2*z+(r**2/s)*x**4/z)) == 0
    assert sp.expand((flux-I)*z
                     -(c*x**2*z**2-I*z+r**2*x**4/s)) == 0
    checked += 1
    m = sp.symbols("m", integer=True, nonnegative=True)
    A = (16**sp.Rational(1, 3)*I)**sp.Rational(1, 3)
    xm = A*16**(-m/sp.Integer(3))
    xn = A*16**(-(m+1)/sp.Integer(3))
    assert sp.simplify(16**m*xm**2*xn-I) == 0
    K = 16**(-m)/sp.Integer(64)
    ym = A*16**sp.Rational(1, 3)*16**(-4*m/sp.Integer(3))/64
    assert sp.simplify(K*xm**2/xn-ym) == 0
    assert sp.simplify(16**(-m/sp.Integer(3))
                       -(4**m)**(-sp.Rational(2, 3))) == 0
    checked += 1
    H = sp.Rational(261*86**2*257, 2)
    raw = sp.Rational(1, 2)*4**(-m)*86**2*16**(-2*m)*257*261
    assert sp.simplify(raw-H*1024**(-m)) == 0
    assert H*sp.Rational(1, 1024**3) < sp.Rational(1, 2)
    checked += 1
    # Finite telescoping of 2*w_j+w_{j+1}=epsilon_j.
    # Decay of the remaining term requires the separate O(j) proof.
    for depth in (1, 2, 5, 8):
        w = sp.symbols(f"w0:{depth+1}")
        rhs = sum(sp.Rational((-1)**j, 2**(j+1))*(2*w[j]+w[j+1])
                  for j in range(depth))
        rhs += sp.Rational((-1)**depth, 2**depth)*w[-1]
        assert sp.expand(rhs-w[0]) == 0
        checked += 1
    return checked


def verify_weakfed_tail_constants():
    """Check envelope arithmetic, not compactness or passage to a limit."""
    m = sp.symbols("m", integer=True, nonnegative=True)
    checked = 0
    for exponent in (sp.Rational(2, 3), sp.Rational(1, 2), sp.Integer(1)):
        ratio = sp.Integer(16)**(-exponent)
        tail = ratio**m/(1-ratio)
        assert sp.simplify(tail-tail.subs(m, m+1)-ratio**m) == 0
        # All three geometric sums have ratio <=1/4.
        assert bool(ratio <= sp.Rational(1, 4))
        assert 48*sp.Rational(4, 3) == 64
        checked += 1
    assert sp.Rational(257+1, 3) == 86
    odd_tail = sp.Integer(86)**2*16**(-2*m)/(1-sp.Rational(1, 16**2))
    assert sp.simplify(odd_tail-odd_tail.subs(m, m+1)
                       -86**2*16**(-2*m)) == 0
    C = sp.Integer(514)
    L = 1+sp.Rational(129, 128)*C
    assert L == sp.Rational(33217, 64)
    assert 1/(2*L) == sp.Rational(32, 33217)
    checked += 1
    return checked


def verify_pair_budget_and_nonclosure():
    """Relabel the original vector field; test only finite local algebra."""
    checked = 0
    nu, force = sp.symbols("nu force", nonnegative=True)
    for m in range(6):
        a = sp.symbols(f"a0:{2*m+4}", nonnegative=True)
        r = sp.symbols(f"r0:{m+1}", positive=True)
        p = lambda i: r[i//2] if i % 2 == 0 else sp.Rational(1, 2**i)
        q = lambda i: sp.Integer(4)**i if i % 2 == 0 else sp.Integer(0)
        rates = []
        for n in (2*m, 2*m+1):
            rhs = -a[n]*(p(n)*a[n+1]+q(n)*a[n+2]+nu*4**n)
            if n:
                rhs += p(n-1)*a[n-1]**2
            if n >= 2:
                rhs += q(n-2)*a[n-2]**2
            if n == 0:
                rhs += force
            rates.append(rhs)
        x, y, z = a[2*m:2*m+3]
        c, s = sp.Integer(16)**m, sp.Rational(1, 2**(2*m+1))
        incoming = (q(2*m-2)*a[2*m-2]**2
                    + p(2*m-1)*a[2*m-1]**2) if m else 0
        pair_x = incoming-x*(c*z+r[m]*y+nu*c)+(force if m == 0 else 0)
        pair_y = r[m]*x**2-s*y*z-4*nu*c*y
        assert sp.expand(rates[0]-pair_x) == 0
        assert sp.expand(rates[1]-pair_y) == 0
        outgoing_flux = (c*x**2+s*y**2)*z
        target = incoming*x-outgoing_flux-nu*c*(x**2+4*y**2)
        if m == 0:
            target += force*x
        assert sp.expand(x*rates[0]+y*rates[1]-target) == 0
        assert sp.expand(rates[1].subs(y, 0)-r[m]*x**2) == 0
        checked += 1
    # Printed MHD second-component vector field is tangent to b=0.
    # No uniqueness or invariance assertion for weak solutions is tested.
    w, kappa, mu, ell = sp.symbols("w kappa mu ell", real=True)
    u, un, b, bn = sp.symbols("u un b bn", real=True)
    mhd_b_rhs = -mu*ell*b-kappa*w*(u*bn-b*un)
    assert mhd_b_rhs.subs({b: 0, bn: 0}) == 0
    checked += 1
    # Two finite-support states with identical norms of every pair.
    state1 = [(0, 0), (sp.Rational(3, 5), sp.Rational(4, 5)),
              (sp.Rational(3, 5), sp.Rational(4, 5)), (0, 0)]
    state2 = [state1[0], (state1[1][1], state1[1][0]),
              state1[2], state1[3]]
    norms = lambda state: [x*x+y*y for x, y in state]
    assert norms(state1) == norms(state2)
    c, s = sp.Integer(16), sp.Rational(1, 8)
    flux = lambda state: (c*state[1][0]**2+s*state[1][1]**2)*state[2][0]
    gap = flux(state2)-flux(state1)
    assert gap == sp.Rational(21, 125)*(c-s)
    assert gap > 0
    checked += 1
    return checked


if __name__ == "__main__":
    print("Closed conservative truncations:", verify_closed_truncations(), "PASS")
    print("Internal cut flux identities:", verify_internal_cuts(), "PASS")
    print("Growth-lemma rescalings:", verify_rescaling(), "PASS")
    print("Coverage edge-length cases:", verify_coverage_constants(), "PASS")
    print("Sparse weak-cut identities:", verify_sparse_cut_support(), "PASS")
    print("Difference-energy identities:", verify_difference_energy(), "PASS")
    print("Stationary elimination cases:", verify_stationary_elimination(), "PASS")
    print("Relative Lyapunov identities:", verify_relative_lyapunov(), "PASS")
    print("Positive coordinate-change cases:", verify_positive_coordinate_changes(), "PASS")
    print("N1 face and support cases:", verify_n1_face_and_support(), "PASS")
    print("Fixed-window scalar checks:", verify_fixed_window_comparison(), "PASS")
    print("Closed-window bound checks:", verify_closed_window_bound(), "PASS")
    print("Strong-path obstruction checks:", verify_strong_path_obstruction(), "PASS")
    print("Weakly-fed outlet checks:", verify_weakly_fed_outlets(), "PASS")
    print("Pair budget and non-closure checks:", verify_pair_budget_and_nonclosure(), "PASS")
    print("Weakly-fed tail and window constants:", verify_weakfed_tail_constants(), "PASS")
    print("Weakly-fed stationary asymptotic identities:", verify_weakfed_stationary_asymptotics(), "PASS")
    print("Feeding perturbation bounds:", verify_feeding_perturbation_bounds(), "PASS")
    print("Scope: exact finite identities only; existence, infinite limits and novelty require separate review.")
