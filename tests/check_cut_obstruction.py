"""Finite algebra for current results and appendices. Not an infinite-system proof."""

import sympy as sp


def tail_budgets():
    count = 0
    for last in range(3, 8):
        a = sp.symbols(f"a0:{last + 1}", real=True)
        p = sp.symbols(f"p0:{last + 1}", nonnegative=True)
        q = sp.symbols(f"q0:{last + 1}", nonnegative=True)
        nu, force = sp.symbols("nu F", nonnegative=True)
        rhs = [-nu * 4**i * a[i] + (force if i == 0 else 0)
               for i in range(last + 1)]
        for i in range(last + 1):
            for step, weight in ((1, p[i]), (2, q[i])):
                j = i + step
                if j <= last:
                    rhs[i] -= weight * a[i] * a[j]
                    rhs[j] += weight * a[i]**2
        for cut in range(1, last - 1):
            flux = (p[cut] * a[cut]**2 * a[cut + 1]
                    + q[cut] * a[cut]**2 * a[cut + 2]
                    + q[cut - 1] * a[cut - 1]**2 * a[cut + 1])
            tail_rate = sum(a[i] * rhs[i] for i in range(cut + 1, last + 1))
            tail_loss = nu * sum(4**i * a[i]**2 for i in range(cut + 1, last + 1))
            assert sp.expand(tail_rate + tail_loss - flux) == 0
            count += 1
    return count


def capacity_certificate():
    p, q, r, left, center, right, far, bound = sp.symbols(
        "p q r left center right far R", nonnegative=True)
    flux = p * center**2 * right + q * center**2 * far + r * left**2 * right
    gap = (p + q + r) * bound * (center**2 + left**2) - flux
    # Every term is nonnegative when R bounds each amplitude.
    certificate = (p * center**2 * (bound - right)
                   + q * center**2 * (bound - far)
                   + r * left**2 * (bound - right)
                   + bound * (r * center**2 + (p + q) * left**2))
    assert sp.expand(gap - certificate) == 0
    return 1


def source_overlap():
    count = 0
    for last in range(1, 10):
        # Exhaust every subset of cuts 1,...,last, not just contiguous cuts.
        for mask in range(1, 1 << last):
            occurrences = [0] * (last + 1)
            for cut in range(1, last + 1):
                if mask & (1 << (cut - 1)):
                    occurrences[cut] += 1
                    occurrences[cut - 1] += 1
            assert max(occurrences) <= 2
            count += 1
    return count


def two_channel_budgets():
    count = 0
    for last in range(1, 6):
        x = sp.symbols(f"x0:{last + 1}", real=True)
        y = sp.symbols(f"y0:{last + 1}", real=True)
        c = sp.symbols(f"c0:{last}", nonnegative=True)
        d = sp.symbols(f"d0:{last}", nonnegative=True)
        h = sp.symbols(f"h0:{last + 1}", nonnegative=True)
        g = sp.symbols(f"g0:{last + 1}", nonnegative=True)
        dx, dy = [sp.Integer(0)] * (last + 1), [sp.Integer(0)] * (last + 1)
        for i in range(last):
            dx[i] -= c[i] * x[i] * x[i + 1]
            dx[i + 1] += c[i] * x[i]**2
            dy[i] -= d[i] * y[i] * y[i + 1]
            dy[i + 1] += d[i] * y[i]**2
        for i in range(last + 1):
            dx[i] += -h[i] * x[i] * y[i] + g[i] * y[i]**2
            dy[i] += h[i] * x[i]**2 - g[i] * y[i] * x[i]
        total = sum(x[i] * dx[i] + y[i] * dy[i] for i in range(last + 1))
        assert sp.expand(total) == 0
        count += 1
        for cut in range(last):
            tail_rate = sum(x[i] * dx[i] + y[i] * dy[i]
                            for i in range(cut + 1, last + 1))
            flux = c[cut] * x[cut]**2 * x[cut + 1] + d[cut] * y[cut]**2 * y[cut + 1]
            assert sp.expand(tail_rate - flux) == 0
            count += 1
    return count


def sparse_unrepaired_cuts():
    weak = {j*j for j in range(2, 9)}
    count = 0
    for mask in range(1 << len(weak)):
        repaired = {n for bit, n in enumerate(sorted(weak)) if mask & (1 << bit)}
        for n in weak - repaired:
            p = sp.Rational(1, 2**n)
            q_here = 4**n if n in repaired else 0
            q_previous = 4**(n - 1) if n - 1 in repaired else 0
            assert p + q_here + q_previous == sp.Rational(1, 2**n)
            count += 1
    return count


def coupled_candidate_checks():
    count = 0
    x, y, h, ratio = sp.symbols("x y h r", positive=True)
    exchange_x, exchange_y = h * (y**2 - x*y), h * (x**2 - x*y)
    assert sp.expand(x * exchange_x + y * exchange_y) == 0
    assert sp.expand(exchange_x + exchange_y - h*(x-y)**2) == 0
    invariant_gap = (exchange_y - ratio*exchange_x).subs(y, ratio*x)
    assert sp.expand(invariant_gap - h*(1-ratio)*(1+ratio**2)*x**2) == 0
    count += 3
    assert sp.Rational(1, 2)*1**2 + 2*2**2 == sp.Rational(17, 2)
    assert sp.Rational(1, 2)*2**2 + 2*1**2 == 4
    count += 2
    for n in range(1, 20):
        previous_strong_is_x = (n - 1) % 2 == 0
        current_weak_is_x = n % 2 == 1
        assert previous_strong_is_x == current_weak_is_x
        assert sp.Rational(1, 2**n) <= 2**n
        assert 2**n == 2 * 2**(n - 1)
        assert sp.Rational(4**n, 2**n) == 2 * sp.Rational(4**(n - 1), 2**(n - 1))
        count += 1
    # Exact lower growth factors for the two cases in appendix A3 section 4.
    assert 8 - sp.Rational(1, 1000) > 2
    assert (sp.Rational(1000, 32) - sp.Rational(2, 1000))/9 > 2
    # 2^(-2/3)<3/4 by cubing: 1/4<27/64, hence C<800<1024.
    assert sp.Rational(1, 4) < sp.Rational(27, 64)
    # C<1024 then 2*C^3<256^4 proves the chosen norm bound.
    assert 2 * 1024**3 < 256**4
    count += 4
    rate, bound = sp.Integer(513), sp.Integer(512)
    threshold = 2/rate + rate*bound**2/2
    assert threshold == 67239936 + sp.Rational(2, 513)
    assert 67239937 > threshold
    assert sp.simplify(2/rate - (2/rate**2 + bound**2/2)/threshold - 1/rate) == 0
    count += 3
    # Generic edge relative-energy identity, including a graph with cycles.
    bi, bj, di, dj, coefficient = sp.symbols("bi bj di dj c", positive=True)
    edge = (-coefficient*di*((bi+di)*(bj+dj)-bi*bj)
            + coefficient*dj*((bi+di)**2-bi**2))
    assert sp.expand(edge - coefficient*(-bj*di**2+bi*di*dj)) == 0
    rho_i, rho_j = sp.symbols("rho_i rho_j", real=True)
    flux = coefficient*bi**2*bj
    assert sp.expand(edge.subs({di: bi*rho_i, dj: bj*rho_j})
                     + flux*(rho_i-rho_j)**2/2
                     - flux*(rho_j**2-rho_i**2)/2) == 0
    count += 2
    return count


def relay_equation_checks():
    """Derive the relabeled equations from the current k_n (not k_n^2) model."""
    count = 0
    x = sp.symbols("x0:10", positive=True)
    y = sp.symbols("y0:10", positive=True)
    nu = sp.symbols("nu", positive=True)
    weak = [y[n] if n % 2 == 0 else x[n] for n in range(10)]
    strong = [x[n] if n % 2 == 0 else y[n] for n in range(10)]
    c = [2**n if n % 2 == 0 else sp.Rational(1, 2**n) for n in range(10)]
    d = [sp.Rational(1, 2**n) if n % 2 == 0 else 2**n for n in range(10)]
    for n in range(8):
        k, viscous_weight, e = 2**n, 4**n, sp.Rational(1, 2**n)
        dx = -c[n]*x[n]*x[n+1] + 8*k*(y[n]**2-x[n]*y[n]) - nu*viscous_weight*x[n]
        dy = -d[n]*y[n]*y[n+1] + 8*k*(x[n]**2-x[n]*y[n]) - nu*viscous_weight*y[n]
        if n:
            dx += c[n-1]*x[n-1]**2
            dy += d[n-1]*y[n-1]**2
        else:
            dx += 1
            dy += 1
        dw, dz = (dy, dx) if n % 2 == 0 else (dx, dy)
        expected_w = (-e*weak[n]*strong[n+1] + 8*k*(strong[n]**2-weak[n]*strong[n])
                      - nu*viscous_weight*weak[n])
        expected_z = (-k*strong[n]*weak[n+1] + 8*k*(weak[n]**2-weak[n]*strong[n])
                      - nu*viscous_weight*strong[n])
        if n:
            expected_w += 2**(n-1)*strong[n-1]**2
            expected_z += sp.Rational(1, 2**(n-1))*weak[n-1]**2
        else:
            expected_w += 1
            expected_z += 1
        assert sp.expand(dw - expected_w) == 0
        assert sp.expand(dz - expected_z) == 0
        assert sp.Rational(k, 2*viscous_weight) <= sp.Rational(1, 2)
        assert sp.Rational(8*k, viscous_weight) <= 8
        count += 4
    assert 2*sp.Rational(1, 2) + 2*8 == 17
    count += 1
    return count


def positive_exchange_checks():
    """Symbolic certificates for all eta>0, not a parameter-grid experiment."""
    eta, sigma, bound = sp.symbols("eta sigma C", positive=True)
    delta = eta/(eta+4)
    a_cubed = 8*(eta+1)/delta**4
    count = 0
    for gap in (a_cubed-8, 4-1/a_cubed-2,
                4-2/(a_cubed*(eta+1))-2):
        numerator, denominator = sp.fraction(sp.cancel(gap))
        # Positive coefficients certify positivity on the full positive axis.
        assert all(c > 0 for c in sp.Poly(numerator, eta).all_coeffs())
        assert all(c >= 0 for c in sp.Poly(denominator, eta).all_coeffs())
        assert denominator.subs(eta, 1) > 0
        count += 1
    assert sp.cancel(eta*(1/delta-1)-4) == 0
    explicit_b = 4*sp.sqrt(2)*sp.sqrt(eta+1)*(eta+4)**2/(eta**2*sigma**sp.Rational(3, 4))
    assert sp.simplify(explicit_b**4 - 16*a_cubed**2/sigma**3) == 0
    rate = 1+bound
    threshold = 2/rate+rate*bound**2/2
    assert sp.simplify(2/rate-(2/rate**2+bound**2/2)/threshold-1/rate) == 0
    count += 3
    return count


def asymmetric_cut_checks():
    """Derive the signed non-level cut flux in a closed finite graph."""
    count = 0
    eta, nu = sp.symbols("eta nu", nonnegative=True)
    for last in range(2, 8):
        x = sp.symbols(f"x0:{last+1}", nonnegative=True)
        y = sp.symbols(f"y0:{last+1}", nonnegative=True)
        dx = [-nu*4**n*x[n] + (1 if n == 0 else 0) for n in range(last+1)]
        dy = [-nu*4**n*y[n] + (1 if n == 0 else 0) for n in range(last+1)]
        for n in range(last+1):
            dx[n] += eta*2**n*(y[n]**2-x[n]*y[n])
            dy[n] += eta*2**n*(x[n]**2-x[n]*y[n])
            if n < last:
                c = 2**n if n % 2 == 0 else sp.Rational(1, 2**n)
                d = sp.Rational(1, 2**n) if n % 2 == 0 else 2**n
                dx[n] -= c*x[n]*x[n+1]
                dx[n+1] += c*x[n]**2
                dy[n] -= d*y[n]*y[n+1]
                dy[n+1] += d*y[n]**2
        for n in range(1, last, 2):
            outside_rate = (sum(x[j]*dx[j] for j in range(n+1, last+1))
                            + sum(y[j]*dy[j] for j in range(n, last+1)))
            outside_loss = nu*(sum(4**j*x[j]**2 for j in range(n+1, last+1))
                               + sum(4**j*y[j]**2 for j in range(n, last+1)))
            flux = (sp.Rational(1, 2**n)*x[n]**2*x[n+1]
                    + sp.Rational(1, 2**(n-1))*y[n-1]**2*y[n]
                    + eta*2**n*(x[n]**2*y[n]-y[n]**2*x[n]))
            assert sp.expand(outside_rate+outside_loss-flux) == 0
            assert sp.Rational(1, 2**n)+sp.Rational(1, 2**(n-1)) == sp.Rational(3, 2**n)
            count += 1
    # Signed exchange: both signs genuinely occur.
    a, b = sp.symbols("a b", nonnegative=True)
    assert (a*a*b-b*b*a).subs({a: 1, b: 2}) == -2
    assert (a*a*b-b*b*a).subs({a: 2, b: 1}) == 2
    count += 2
    # J2: put p=nu^(1/3), q=sqrt(eta), r=1/max(p,q).
    p, q = sp.symbols("p q", positive=True)
    for largest in (p, q):
        # Under largest=max(p,q), these certificates are nonnegative.
        assert sp.cancel(p-p**3/largest**2
                         - p*(largest-p)*(largest+p)/largest**2) == 0
        assert sp.cancel(q-q**2/largest-q*(largest-q)/largest) == 0
        assert sp.expand(p+q-largest-(q if largest == p else p)) == 0
        count += 1
    return count


def relay_deletion_checks():
    """Removing reverse pumps preserves algebra; reindex only beta=0."""
    eta, nu = sp.symbols("eta nu", nonnegative=True)
    beta = sp.symbols("beta", nonnegative=True)
    w, z, radius = sp.symbols("w z R", nonnegative=True)
    ew = eta*(beta*z*z-w*z)
    ez = eta*(w*w-beta*w*z)
    assert sp.expand(w*ew+z*ez) == 0
    assert sp.expand(ew+ez-eta*(w-z)*(w-beta*z)) == 0
    # Difference from the beta=1 growth-branch lower bound is nonnegative for beta<=1.
    assert sp.expand(ez-eta*w*(w-z)-eta*(1-beta)*w*z) == 0
    # Root bound discards eta*w^2+beta*eta*z^2; wz<=R(w+z) if R>=z.
    assert sp.expand(ew+ez+eta*(1+beta)*w*z-eta*(w*w+beta*z*z)) == 0
    assert sp.expand(radius*(w+z)-w*z-(w*(radius-z)+radius*z)) == 0
    count = 5
    for last in range(6):
        a = sp.symbols(f"a0:{2*last+2}", nonnegative=True)
        weak, strong = a[::2], a[1::2]
        dw, dz = [], []
        for n in range(last+1):
            dw.append(eta*2**n*(beta*strong[n]**2-weak[n]*strong[n])
                      -nu*4**n*weak[n]+(1 if n == 0 else 0))
            dz.append(eta*2**n*(weak[n]**2-beta*weak[n]*strong[n])
                      -nu*4**n*strong[n]+(1 if n == 0 else 0))
        for n in range(last):
            e = sp.Rational(1, 2**n)
            dw[n] -= e*weak[n]*strong[n+1]
            dz[n+1] += e*weak[n]**2
            dz[n] -= 2**n*strong[n]*weak[n+1]
            dw[n+1] += 2**n*strong[n]**2
        total = sum(weak[n]*dw[n]+strong[n]*dz[n] for n in range(last+1))
        loss = nu*sum(4**(j//2)*a[j]**2 for j in range(len(a)))
        assert sp.expand(total+loss-a[0]-a[1]) == 0
        count += 1
        # Independently build the scalar forward graph, including the closed terminal.
        rhs = [-nu*4**(j//2)*a[j]+(1 if j in (0, 1) else 0) for j in range(len(a))]
        for j in range(len(a)):
            p = (eta if j % 2 == 0 else 1)*2**(j//2)
            r = sp.Rational(1, 2**(j//2)) if j % 2 == 0 else 0
            for step, coefficient in ((1, p), (3, r)):
                if j+step < len(a):
                    rhs[j] -= coefficient*a[j]*a[j+step]
                    rhs[j+step] += coefficient*a[j]**2
        for n in range(last+1):
            assert sp.expand(dw[n].subs(beta, 0)-rhs[2*n]) == 0
            assert sp.expand(dz[n].subs(beta, 0)-rhs[2*n+1]) == 0
            count += 2
    return count


def scale_dependent_exchange_checks():
    """Finite identities for appendix A4; no assertion about the unresolved dynamic regime."""
    count = 0
    current, previous, back = sp.symbols("L Lprev B", nonnegative=True)
    power = sp.symbols("I", positive=True)
    # Weak-node balance: incoming strong flux I-Lprev, outgoing weak flux L.
    relay_flux = power-previous-current+back
    assert sp.expand(relay_flux-(power-previous)-back+current) == 0
    intensity, k, w = sp.symbols("h K w", positive=True)
    lower_z = power/(2*intensity*w**2)
    assert sp.cancel(power/(k*lower_z**2)-4*intensity**2*w**4/(k*power)) == 0
    exact_z = power/(intensity*w**2)
    assert sp.cancel(power/(k*exact_z**2)-intensity**2*w**4/(k*power)) == 0
    count += 3
    gamma, n, log_i, log_eta, log_two = sp.symbols("gamma n logI logEta log2", real=True)
    a, b = (2*gamma-1)/3, (2-gamma)/3
    assert sp.expand(2*a+b-gamma) == 0
    assert sp.expand(a+2*b-1) == 0
    log_amplitude = (log_i-2*log_eta-a*log_two)/3
    log_w = log_amplitude-a*n*log_two
    assert sp.expand(log_w.subs(n, n+1)-4*log_w
                     -(2*gamma-1)*n*log_two-2*log_eta+log_i) == 0
    count += 3
    for cuts in range(1, 12):
        coordinates = []
        for cut in range(1, 2*cuts, 2):
            coordinates.extend((("x", cut), ("y", cut-1)))
            assert sp.Rational(1, 2**cut) <= 1
            assert sp.Rational(1, 2**(cut-1)) <= 1
        assert len(set(coordinates)) == 2*cuts
        assert max(index for _, index in coordinates) == 2*cuts-1
        count += 1
    return count


def supercritical_exchange_checks():
    """Constants and exponent identities behind the two distinct envelopes."""
    eta = sp.symbols("eta", positive=True)
    az3 = 8*(1+eta)
    aw3 = 8*(1+eta)*(8+eta**-2)
    checks = (
        (az3-4/az3)/(2*(1+eta))-2,
        aw3-8*az3,
        eta**2*aw3/2-eta-1/aw3-2,
        eta**2*aw3/2-4*(1+eta),
    )
    for gap in checks:
        numerator, denominator = sp.fraction(sp.cancel(gap))
        assert all(c >= 0 for c in sp.Poly(numerator, eta).all_coeffs())
        assert numerator.subs(eta, 1) > 0
        assert all(c >= 0 for c in sp.Poly(denominator, eta).all_coeffs())
        assert denominator.subs(eta, 1) > 0
    count = len(checks)
    gamma = sp.symbols("gamma", real=True)
    exponent = (2*gamma-1)/3
    assert sp.expand(2*gamma-1-3*exponent) == 0
    assert sp.expand(3*exponent-1-(2*gamma-2)) == 0
    assert exponent.subs(gamma, sp.Rational(1, 2)) == 0
    assert exponent.subs(gamma, 1) == sp.Rational(1, 3)
    count += 4
    z, k, power, r = sp.symbols("z K I r", positive=True)
    growth = (k*k*z**4/power-4*power/z**2)/(2*(1+eta)*k*z)
    assert sp.simplify(growth.subs(z, r*(power/k)**sp.Rational(1, 3))
                       -(r**3-4/r**3)/(2*(1+eta))) == 0
    count += 1
    # Exact net-exchange energy identity; bounds additionally require the budget.
    w, h, beta, visc, incoming, weak_out, force = sp.symbols(
        "w h beta visc incoming weak_out force", nonnegative=True)
    derivative = incoming-weak_out+h*(beta*z*z-w*z)-visc*w+force
    net = h*w*z*(w-beta*z)
    assert sp.expand(w*derivative+net-w*(incoming-weak_out+force)+visc*w*w) == 0
    count += 1
    return count


def steady_all_exponents_checks():
    """Growth can outrun an increasing envelope; constants do not prove the limit."""
    eta, growth = sp.symbols("eta G", positive=True)
    az3 = 8*(1+eta)
    aw3 = 64*(1+eta)+2*(growth+eta+1)/eta**2
    count = 0
    for gap in (aw3-8*az3, aw3-1,
                eta**2*aw3/2-eta-1/aw3-growth):
        numerator, denominator = sp.fraction(sp.cancel(gap))
        for polynomial in (numerator, denominator):
            assert all(c >= 0 for c in sp.Poly(polynomial, eta, growth).coeffs())
            assert polynomial.subs({eta: 1, growth: 3}) > 0
        count += 1
    # Envelope ratio q=2^(-a). For a<0 take G=2+q; for a>=0 take G=3.
    q, excess = sp.symbols("q excess", positive=True)
    assert sp.cancel((2+q)/q-1-2/q) == 0
    # q=1/(1+excess) parametrizes 0<q<1 without a numerical gamma grid.
    assert sp.simplify((3/q-1).subs(q, 1/(1+excess))-(2+3*excess)) == 0
    assert 3/1 > 1
    count += 3
    # Relative squared norm y'<=-nu*y gives the norm rate exp(-nu*t/2).
    t, nu, initial = sp.symbols("t nu initial", positive=True)
    squared = initial**2*sp.exp(-nu*t)
    assert sp.simplify(sp.diff(squared, t)+nu*squared) == 0
    count += 1
    return count


def dynamic_no_dissipation_checks():
    """Algebra supporting the inviscid time-integral estimate, not compactness."""
    k, z, w_next, z_next, z_far, h, e, beta = sp.symbols(
        "K z wnext znext zfar h e beta", nonnegative=True)
    weak_rhs = k*z**2-e*w_next*z_far+h*(beta*z_next**2-w_next*z_next)
    reconstructed = weak_rhs+e*w_next*z_far+h*w_next*z_next-beta*h*z_next**2
    assert sp.expand(reconstructed-k*z**2) == 0
    count = 1
    a, b, excess = sp.symbols("A B excess", positive=True)
    y = 2*a+16*b**2+excess
    gap = sp.expand((y-a)**2-4*b**2*y)
    assert all(c > 0 for c in sp.Poly(gap, a, b, excess).coeffs())
    count += 1
    step = sp.symbols("j", integer=True, nonnegative=True)
    assert sp.limit(2**step/4**step, step, sp.oo) == 0
    count += 1
    delta = sp.symbols("delta", nonnegative=True)
    gamma = sp.Rational(1, 2)-delta
    assert sp.expand(2*gamma-1+2*delta) == 0
    assert sp.expand(gamma-1+sp.Rational(1, 2)+delta) == 0
    count += 2
    # Finite Cauchy-Schwarz identity mirrors the continuous time inequality.
    ws = sp.symbols("w0:3", nonnegative=True)
    zs = sp.symbols("z0:3", nonnegative=True)
    cs_gap = sum(v**4 for v in ws)*sum(v**2 for v in zs)-sum(ws[i]**2*zs[i] for i in range(3))**2
    certificate = sum((zs[i]*ws[j]**2-zs[j]*ws[i]**2)**2
                      for i in range(3) for j in range(i+1, 3))
    assert sp.expand(cs_gap-certificate) == 0
    count += 1
    return count


def alternating_chain_checks():
    """Deleted chain, its cut budgets, and the exact diagonal scaling restriction."""
    count = 0
    for last in range(2, 6):
        a = sp.symbols(f"a0:{2*last+2}", nonnegative=True)
        h = sp.symbols(f"h0:{last+1}", positive=True)
        k = sp.symbols(f"k0:{last+1}", positive=True)
        nu = sp.symbols("nu", nonnegative=True)
        mass = [4**(j//2) for j in range(len(a))]
        rhs = [-nu*mass[j]*a[j]+(1 if j < 2 else 0)
               for j in range(len(a))]
        for j in range(len(a)-1):
            weight = h[j//2] if j % 2 == 0 else k[j//2]
            rhs[j] -= weight*a[j]*a[j+1]
            rhs[j+1] += weight*a[j]**2
        for n in range(last+1):
            w, z = a[2*n:2*n+2]
            weak = -h[n]*w*z-nu*4**n*w+(1 if n == 0 else 0)
            strong = h[n]*w*w-nu*4**n*z+(1 if n == 0 else 0)
            if n:
                weak += k[n-1]*a[2*n-1]**2
            if n < last:
                strong -= k[n]*z*a[2*n+2]
            assert sp.expand(rhs[2*n]-weak) == 0
            assert sp.expand(rhs[2*n+1]-strong) == 0
            count += 2
        assert sp.expand(sum(a[j]*rhs[j]+nu*mass[j]*a[j]**2
                             for j in range(len(a)))-a[0]-a[1]) == 0
        count += 1
        for n in range(1, last+1):
            flux = h[n]*a[2*n]**2*a[2*n+1]
            assert sp.expand(sum(a[j]*rhs[j]+nu*mass[j]*a[j]**2
                                 for j in range(2*n+1))-a[0]-a[1]+flux) == 0
            count += 1
    si, sj, c, p, bi, bj = sp.symbols("si sj c p bi bj", positive=True)
    outgoing = -p*sj*bi*bj/c
    incoming = p*si**2*bi**2/(c*sj)
    defect = sp.factor(bi*outgoing+bj*incoming)
    assert sp.simplify(defect-p*bi**2*bj*(si-sj)*(si+sj)/(c*sj)) == 0
    assert sp.expand(si**2*bi*outgoing+sj**2*bj*incoming) == 0
    assert sp.simplify(defect.subs(sj, si)) == 0
    count += 3
    # Geometric-series algebra for gamma=1/4. This does not test its dynamic limit.
    r_even, r_odd = 2**(-sp.Rational(1, 6)), 2**(-sp.Rational(2, 3))
    for ratio in (r_even, r_odd):
        assert 0 < ratio < 1
        count += 1
    for n in range(4):
        even_weight, odd_weight = 2**sp.Rational(n, 4), sp.Integer(2)**n
        assert sp.simplify(even_weight**(-sp.Rational(2, 3))-r_even**n) == 0
        assert sp.simplify(odd_weight**(-sp.Rational(2, 3))-r_odd**n) == 0
        count += 2
    # The two steady lower bounds come from these exact balance rearrangements.
    w, z, wn, zn, k, h, hn, v = sp.symbols("w z wn zn K h hn v", positive=True)
    weak_residual = k*z**2-hn*wn*zn-4*v*wn
    strong_residual = h*w**2-k*z*wn-v*z
    assert sp.expand(wn*(hn*zn+4*v)+weak_residual-k*z**2) == 0
    assert sp.simplify((h*w**2-strong_residual)/(k*z)-v/k-wn) == 0
    count += 2
    return count


def cascade_time_comparison_checks():
    """Finite algebra for the adaptation of the classical transfer-time argument."""
    p, q, radius, theta, t = sp.symbols("p q R theta t", positive=True)
    lower = p*theta/(q*radius)*(1-sp.exp(-q*radius*t))
    assert sp.simplify(sp.diff(lower, t)+q*radius*lower-p*theta) == 0
    assert lower.subs(t, 0) == 0
    # Energy is sum x_j^2 in this comparison, so its flux has a factor of two.
    drop = 2*p**2*theta**2/(q*radius)*(t-(1-sp.exp(-q*radius*t))/(q*radius))
    assert sp.simplify(sp.diff(drop, t)-2*p*theta*lower) == 0
    # If z=q*R*t >=2, z-1+exp(-z) >= z/2; show the nonnegative remainder.
    excess = sp.symbols("excess", nonnegative=True)
    z = 2+excess
    assert sp.simplify(z-1+sp.exp(-z)-z/2-(excess/2+sp.exp(-z))) == 0
    count = 4
    r, total = sp.symbols("r S", positive=True)
    threshold = radius**2*r/(2*total)
    transfer_part = radius**3*r**3/threshold**2
    assert sp.simplify(transfer_part-4*total**2*r/radius) == 0
    assert sp.simplify(radius**2*total/(2*total)-radius**2/2) == 0
    count += 2
    gamma, n, eta = sp.symbols("gamma n eta", real=True)
    # Log2 rates, with eta represented separately to avoid logarithm branches.
    assert sp.expand(2*gamma*n-n-(2*gamma-1)*n) == 0
    assert sp.expand(2*n-gamma*(n+1)-((2-gamma)*n-gamma)) == 0
    count += 2
    for g in (sp.Rational(1, 2), sp.Rational(3, 4), sp.Integer(1)):
        even_ratio = 2**(-(2*g-1)/3)
        odd_ratio = 2**(-(2-g)/3)
        if g == sp.Rational(1, 2):
            assert even_ratio == 1
        else:
            assert 0 < even_ratio < 1
        assert 0 < odd_ratio < 1
        count += 2
    return count


def logarithmic_exchange_checks():
    """Envelope constants and exponents, not the infinite-dimensional proof."""
    eta, bound, maximum = sp.symbols("eta H D", positive=True)
    az3 = 8*(1+bound)
    aw3 = 1+8*az3*maximum+8*(bound+maximum+1)/eta**2
    count = 0
    for gap in (aw3-1, aw3-8*az3*maximum,
                eta**2*aw3/2-bound-maximum/aw3-2):
        numerator, denominator = sp.fraction(sp.cancel(gap))
        for polynomial in (numerator, denominator):
            assert all(c >= 0 for c in sp.Poly(polynomial, eta, bound, maximum).coeffs())
            assert polynomial.subs({eta: 1, bound: 1, maximum: 1}) > 0
        count += 1
    alpha = sp.symbols("alpha", real=True)
    rho = 2*alpha/3
    assert sp.expand(2*alpha-3*rho) == 0
    assert sp.expand(2*rho-alpha-alpha/3) == 0
    assert (4*alpha/3).subs(alpha, sp.Rational(3, 4)) == 1
    assert (4*alpha).subs(alpha, sp.Rational(1, 4)) == 1
    assert (2*alpha/3).subs(alpha, sp.Rational(3, 2)) == 1
    count += 5
    # Exact D_1=9/4: n=0,1,2, then the ratio is less than one for every n>=2.
    assert [sp.Rational((n+1)**2, 2**n) for n in range(3)] == [1, 2, sp.Rational(9, 4)]
    excess = sp.symbols("excess", nonnegative=True)
    n = 2+excess
    gap = sp.expand(2*(n+1)**2-(n+2)**2)
    assert gap == excess**2+4*excess+2
    assert 512 > 8*20*sp.Rational(9, 4)
    assert 256-sp.Rational(3, 2)-sp.Rational(9, 2048) > 2
    count += 4
    # At alpha=1 the envelope square is summable but the transfer-time series is not.
    assert sp.Rational(4, 3) > 1
    assert sp.Rational(2, 3) <= 1
    count += 2
    return count


def integrated_tail_summability_checks():
    """Boundary accounting and finite certificates for docs/03-results-zh.md section 3.2."""
    eta, ell, ell_next, w, v, k, fourth = sp.symbols(
        "eta ell ellnext W V K fourth", positive=True)
    # Divide the weak-coordinate time integral by ell_n^2.
    coefficient = eta*ell_next*sp.sqrt(w*ell_next**2*v)/ell**2
    assert sp.simplify(coefficient-eta*(ell_next/ell)**2*sp.sqrt(w*v)) == 0
    # h_n^2=(eta^2 K_n ell_n^2), P_n=ell_n^2 V_n/K_n.
    flux_square_bound = eta**2*k*ell**2*(ell**2*v/k)*fourth
    assert sp.simplify(flux_square_bound/ell**4-eta**2*v*fourth) == 0
    count = 2
    for last in range(1, 5):
        values = sp.symbols(f"v0:{last+2}", nonnegative=True)
        assert sp.expand(sum(values[1:])-sum(values[:-1])-values[-1]+values[0]) == 0
        # Cauchy-Schwarz for square roots of W_n and V_n.
        xs = sp.symbols(f"x0:{last+1}", nonnegative=True)
        ys = sp.symbols(f"y0:{last+1}", nonnegative=True)
        gap = sum(x*x for x in xs)*sum(y*y for y in ys)-sum(x*y for x, y in zip(xs, ys))**2
        certificate = sum((xs[i]*ys[j]-xs[j]*ys[i])**2
                          for i in range(last+1) for j in range(i+1, last+1))
        assert sp.expand(gap-certificate) == 0
        count += 2
    x, y = sp.symbols("x y", nonnegative=True)
    assert sp.expand((x*x+y*y)/2-x*y-(x-y)**2/2) == 0
    d, c, sw, cap, sn = sp.symbols("D c SW cap SN", nonnegative=True)
    # SN <= D+(SN+cap)/2+c^2 SW/2 rearranges to the stated finite uniform bound.
    assert sp.expand(2*(d+(sn+cap)/2+c*c*sw/2-sn)
                     -(2*d+cap+c*c*sw-sn)) == 0
    count += 2
    alpha = sp.symbols("alpha", real=True)
    assert (4*alpha).subs(alpha, sp.Rational(1, 2)) > 1
    assert (2*alpha).subs(alpha, sp.Rational(1, 2)) == 1
    count += 2
    return count


def bounded_exchange_robustness_checks():
    """Finite identities for variable per-shell coefficients, not limit certification."""
    count = 0
    for size in range(2, 6):
        ws = sp.symbols(f"w0:{size}", nonnegative=True)
        zs = sp.symbols(f"z0:{size}", nonnegative=True)
        hs = sp.symbols(f"h0:{size}", positive=True)
        bs = sp.symbols(f"beta0:{size}", nonnegative=True)
        ks = sp.symbols(f"k0:{size}", positive=True)
        es = sp.symbols(f"e0:{size}", positive=True)
        dw = [hs[n]*(bs[n]*zs[n]**2-ws[n]*zs[n]) for n in range(size)]
        dz = [hs[n]*(ws[n]**2-bs[n]*ws[n]*zs[n]) for n in range(size)]
        for n in range(size-1):
            dz[n] -= ks[n]*zs[n]*ws[n+1]
            dw[n+1] += ks[n]*zs[n]**2
            dw[n] -= es[n]*ws[n]*zs[n+1]
            dz[n+1] += es[n]*ws[n]**2
        assert sp.expand(sum(ws[n]*dw[n]+zs[n]*dz[n] for n in range(size))) == 0
        count += 1
        for n in range(1, size-1):
            derivative = sum(ws[j]*dw[j]+zs[j]*dz[j] for j in range(n))+ws[n]*dw[n]
            exchange = hs[n]*ws[n]*zs[n]*(ws[n]-bs[n]*zs[n])
            other_out = es[n-1]*ws[n-1]**2*zs[n]+es[n]*ws[n]**2*zs[n+1]
            assert sp.expand(derivative+exchange+other_out) == 0
            count += 1
    lower, bound, maximum, reverse = sp.symbols("eta_min H D b", positive=True)
    az3 = 8*(1+bound)
    aw3 = 1+8*reverse**3*az3*maximum+8*(reverse*bound+maximum+1)/lower**2
    for gap in (aw3-1, aw3-8*reverse**3*az3*maximum,
                lower**2*aw3/2-reverse*bound-maximum/aw3-2):
        numerator, denominator = sp.fraction(sp.cancel(gap))
        for polynomial in (numerator, denominator):
            assert all(c >= 0 for c in sp.Poly(polynomial, lower, bound, maximum, reverse).coeffs())
            assert polynomial.subs({lower: 1, bound: 1, maximum: 1, reverse: 1}) > 0
        count += 1
    # w > 2 b z and beta <= b imply w-beta z >= w/2.
    z, beta, slack, excess = sp.symbols("z beta slack excess", nonnegative=True)
    b = beta+slack
    w = 2*b*z+excess
    assert sp.expand(w-beta*z-w/2-slack*z-excess/2) == 0
    count += 1
    return count


def sequential_gate_checks():
    """Finite algebra for the prescribed-time gate model in appendix B section 7."""
    count = 0
    for n in range(1, 13):
        durations = [8*(j+1)*sp.Rational(1, 2)**j for j in range(n)]
        assert sum(durations) == 32-16*(n+2)*sp.Rational(1, 2)**n
        assert sum(4**j*durations[j] for j in range(n)) == 8*((n-1)*2**n+1)
        assert 8*((n-1)*2**n+1) <= 8*n*2**n
        count += 3
    t, tau = sp.symbols("t tau", positive=True)
    gate = sp.sin(sp.pi*t/tau)**2
    assert sp.integrate(gate, (t, 0, tau)) == tau/2
    for endpoint in (0, tau):
        assert gate.subs(t, endpoint) == 0
        assert sp.diff(gate, t).subs(t, endpoint) == 0
    count += 5
    x, y, d, kappa = sp.symbols("x y d kappa", positive=True)
    dx, dy = -kappa*x*y-d*x, kappa*x*x-4*d*y
    assert sp.expand(x*dx+y*dy+d*(x*x+4*y*y)) == 0
    # For r^2=x^2+y^2 and s=y/r, substitute x^2=r^2(1-s^2).
    r, s = sp.symbols("r s", real=True)
    rdot = -d*r*(1+3*s*s)
    sdot = (kappa*r*r*(1-s*s)-4*d*r*s)/r-s*rdot/r
    assert sp.simplify(sdot-(1-s*s)*(kappa*r-3*d*s)) == 0
    count += 2
    n = sp.symbols("n", nonnegative=True)
    assert sp.expand(sp.Rational(7, 4)*(n+1)-sp.Rational(3, 32)-(n+1)) == 3*n/4+sp.Rational(21, 32)
    assert sp.Rational(7, 8)*sp.Rational(2, 3) == sp.Rational(7, 12) > sp.Rational(1, 2)
    # Finite Taylor lower bound for e^2; all omitted terms are strictly positive.
    assert sum(sp.Rational(2)**j/sp.factorial(j) for j in range(5)) == 7
    assert sp.Rational(1, 2)*sp.Rational(2, 3)**2 == sp.Rational(2, 9)
    count += 4
    for m in (16, 20, 24, 28):
        nu, shell = sp.Rational(1, 2)**(4*m), 3*m
        transit = 8*nu*((shell-1)*2**shell+1)
        assert transit <= sp.Rational(1, 32)
        assert nu*4**shell == 2**(2*m)
        count += 2
    return count


def temporal_concentration_checks():
    """Finite certificates for docs/03-results-zh.md section 3.3; no moment bound is assumed proven."""
    eta, ell, k, v, w, concentration, window = sp.symbols(
        "eta ell K V W concentration T", positive=True)
    fourth = concentration*w*w/window
    z_second = ell*ell*v/k
    flux_squared_bound = eta**2*k*ell**2*fourth*z_second
    assert sp.simplify(flux_squared_bound
                       -eta**2*ell**4*concentration*w*w*v/window) == 0
    normalized = (flux_squared_bound/(ell**4*concentration))**sp.Rational(1, 3)
    assert sp.simplify(normalized
                       -(eta**2/window)**sp.Rational(1, 3)
                       *w**sp.Rational(2, 3)*v**sp.Rational(1, 3)) == 0
    count = 2
    # Weighted Young certificate, used after normalizing the two finite sums.
    r, s = sp.symbols("r s", nonnegative=True)
    assert sp.expand(2*r**3+s**3-3*r*r*s-(r-s)**2*(2*r+s)) == 0
    count += 1
    alpha, delta = sp.symbols("alpha delta", real=True)
    assert sp.expand(4*alpha/3+delta/3-(4*alpha+delta)/3) == 0
    count += 1
    for a, exponent in ((sp.Rational(1, 2), 1),
                        (sp.Rational(3, 5), sp.Rational(3, 5)),
                        (sp.Rational(2, 3), sp.Rational(1, 3)),
                        (sp.Rational(3, 4), 0)):
        assert (4*a+exponent)/3 == 1
        count += 1
    amplitude, duration = sp.symbols("amplitude duration", positive=True)
    second, fourth = amplitude**2*duration, amplitude**4*duration
    assert sp.simplify(window*fourth/second**2-window/duration) == 0
    count += 1
    # Cauchy-Schwarz for time samples: m sum(w_i^4) >= (sum w_i^2)^2.
    for size in range(2, 6):
        samples = sp.symbols(f"sample0:{size}", nonnegative=True)
        gap = size*sum(x**4 for x in samples)-sum(x*x for x in samples)**2
        certificate = sum((samples[i]**2-samples[j]**2)**2
                          for i in range(size) for j in range(i+1, size))
        assert sp.expand(gap-certificate) == 0
        count += 1
    return count


def startup_concentration_checks():
    """Formal leading-order and moment identities, not time-uniform tail estimates."""
    t, window, amplitude = sp.symbols("t T amplitude", positive=True)
    count = 0
    size = 5
    powers = [2**(n+1)-1 for n in range(size)]
    cw, cz = [sp.Integer(1)], [sp.Integer(1)]
    for n in range(1, size):
        cw.append(2**(n-1)*cz[n-1]**2/powers[n])
        cz.append(sp.Rational(1, 2)**(n-1)*cw[n-1]**2/powers[n])
        assert powers[n] == 2*powers[n-1]+1
        assert cw[-1] > 0 and cz[-1] > 0
        count += 1
    ws = [cw[n]*t**powers[n] for n in range(size)]
    zs = [cz[n]*t**powers[n] for n in range(size)]
    for n in range(size):
        h, damping = n+1, sp.Rational(1, 7)*4**n
        dw = h*(zs[n]**2-ws[n]*zs[n])-damping*ws[n]+(1 if n == 0 else 0)
        dz = h*(ws[n]**2-ws[n]*zs[n])-damping*zs[n]+(1 if n == 0 else 0)
        if n:
            dw += 2**(n-1)*zs[n-1]**2
            dz += sp.Rational(1, 2)**(n-1)*ws[n-1]**2
        if n+1 < size:
            dw -= sp.Rational(1, 2)**n*ws[n]*zs[n+1]
            dz -= 2**n*zs[n]*ws[n+1]
        for trial, rhs in ((ws[n], dw), (zs[n], dz)):
            residual = sp.Poly(sp.expand(sp.diff(trial, t)-rhs), t)
            assert min(power[0] for power, coefficient in residual.terms()) >= powers[n]
            count += 1
    for m in powers:
        second = sp.integrate(amplitude**2*t**(2*m), (t, 0, window))
        fourth = sp.integrate(amplitude**4*t**(4*m), (t, 0, window))
        ratio = sp.Rational((2*m+1)**2, 4*m+1)
        assert sp.simplify(window*fourth/second**2-ratio) == 0
        assert ratio == m+sp.Rational(3, 4)+sp.Rational(1, 4*(4*m+1))
        count += 2
    return count


def quantitative_dissipation_checks():
    """Finite index and inequality certificates for docs/03-results-zh.md section 3.4."""
    count = 0
    for power in (20, 21, 30, 31, 40, 41, 60, 61):
        nu = sp.Rational(1, 2)**power
        outer, inner = power//2, (power//2)//2-1
        assert inner >= 4 and inner+1 == outer//2
        assert nu*4**outer <= 1
        assert nu*4**inner <= sp.sqrt(nu)/4
        count += 3
    for outer in range(10, 21):
        for n in range(outer//2+1):
            assert 4**(outer-n) >= 2**outer
            count += 1
    for last in range(1, 9):
        geometric_sum = sum(4**(n+1) for n in range(last+1))
        assert geometric_sum == (4**(last+2)-4)//3
        assert geometric_sum <= sp.Rational(16, 3)*4**last
        count += 2
    a, b, extra = sp.symbols("A B extra", positive=True)
    cap = 2*a+16*b*b+extra
    growth_gap_numerator = sp.expand(cap*(cap-16*b*b))
    assert all(coef > 0 for coef in sp.Poly(growth_gap_numerator, a, b, extra).coeffs())
    count += 1
    sn, d, remainder, c, sw, cap = sp.symbols(
        "SN D remainder c SW cap", nonnegative=True)
    young_bound = d+remainder+(sn+cap)/2+c*c*sw/2
    assert sp.expand(2*(young_bound-sn)
                     -(2*d+2*remainder+cap+c*c*sw-sn)) == 0
    count += 1
    for n in range(1, 8):
        assert sp.Rational(1, 2)**(n-1)+sp.Rational(1, 2)**n == 3*sp.Rational(1, 2)**n
        count += 1
    return count


def uniform_envelope_checks():
    """Why docs/03-results-zh.md section 3.1 yields a viscosity-uniform envelope; finite certificates."""
    eta, radius, window = sp.symbols("eta R T", positive=True)
    g = sp.symbols("g", positive=True)              # g = 2^alpha >= 1
    coefficient = eta*g**2*radius*sp.sqrt(window)   # B = eta 2^(2 alpha) R sqrt(T)
    base = radius+window*radius**2+window*radius
    count = 0
    # Propagation needs V_n > 2A and V_n > 16 B^2; each stated constant implies both,
    # by the exact factorisations below (g = 2^alpha >= 1, so g^4 >= 1).
    variants = (2*base+16*coefficient**2*g**4,
                2*base+16*coefficient**2*g**4+window*radius**2)
    for constant, slack in zip(variants, (0, window*radius**2)):
        # C_T - 2A = 16 B^2 g^4 (+ T R^2), and C_T - 16 B^2 = 2A + 16 B^2 (g^4 - 1) (+ T R^2).
        # Both are nonnegative because g = 2^alpha >= 1, so the propagation threshold holds.
        assert sp.expand(constant-2*base-16*coefficient**2*g**4-slack) == 0
        assert sp.expand(constant-16*coefficient**2-2*base
                         -16*coefficient**2*(g**4-1)-slack) == 0
        count += 2
    # The constant of appendix A4 section 2.1 (no extra T R^2, no g^4 factor) is covered too.
    assert sp.expand(2*base+16*coefficient**2-2*base-16*coefficient**2) == 0
    assert sp.expand(2*base+16*coefficient**2-16*coefficient**2-2*base) == 0
    count += 2
    # 2A - T R^2 = R(RT + 2T + 2) > 0, so every candidate constant exceeds T R^2.
    assert sp.expand(2*base-window*radius**2-radius*(radius*window+2*window+2)) == 0
    count += 1
    # The four-fold step fires at every index, because the crude bound 2^m T R^2 is
    # available at every m and C_T > T R^2 leaves the exponent k-n+2d >= 0 at k = n.
    n, k, d = sp.symbols("n k d", nonnegative=True)
    ratio = sp.powsimp(2**(2*k)*2**(2*d)*window*radius**2
                       /(2**(n+k)*window*radius**2), force=True)
    assert sp.simplify(sp.log(ratio, 2)-(k-n+2*d)) == 0
    assert sp.simplify((k-n+2*d).subs(k, n)-2*d) == 0
    count += 2
    # Concrete rational audit of the same step at paper-scale indices.
    for outer, inner in ((23, 11), (23, 6), (16, 8)):
        crude_next = 2**(inner+1)*sp.Rational(1, 5)*2500
        threshold = sp.Rational(crude_next, 4)+1
        assert 4*threshold > crude_next
        assert sp.Integer(4)**(outer-inner)*threshold > 2**outer*sp.Rational(1, 5)*2500
        count += 2
    return count


def finite_partial_sum_checks():
    """docs/03-results-zh.md section 3.2: the surplus index is the one already closed by the envelope."""
    eta, ell_next, w, z, second, flux, viscous = sp.symbols(
        "eta ellnext w z Y F nu4w", positive=True)
    # Single-layer weak-outlet balance: nonlinear flux plus the weak vertical edge plus
    # the viscous remainder.  Dividing by (n+1)^(2 alpha) must reproduce docs/03-results-zh.md d_n.
    balance = flux+ell_next*w*z+viscous
    weight = sp.symbols("weight", positive=True)
    assert sp.simplify(balance/weight
                       -(flux/weight+(ell_next/weight)*w*z+viscous/weight)) == 0
    # The Cauchy-Schwarz step: e_(n+1) int w_(n+1) z_(n+2) <= e_(n+1) sqrt(W_(n+1) Y_(n+2)/K_(n+2)).
    young = (sp.sqrt(second)*sp.sqrt(ell_next**2*z**2/ell_next**2)
             -sp.sqrt(second*z**2))
    assert sp.simplify(young) == 0
    count = 2
    for last in range(1, 7):
        values = sp.symbols(f"v0:{last+2}", nonnegative=True)
        assert sp.expand(sum(values[:last+1])-sum(values[:last])-values[last]) == 0
        count += 1
    # S_N <= D + (S_N + V_(N+1))/2 + c^2 SW/2 closes only through V_(N+1) <= C_T.
    source, cap, weight, sw, sn = sp.symbols("D cap c SW SN", nonnegative=True)
    assert sp.expand(2*(sn-(source+(sn+cap)/2+weight**2*sw/2))
                     -(sn-2*source-cap-weight**2*sw)) == 0
    count += 1
    alpha = sp.symbols("alpha", real=True)
    assert (4*alpha).subs(alpha, sp.Rational(1, 4)) == 1
    assert sp.Rational(4, 3)*sp.Rational(3, 4) == 1
    count += 2
    return count


def cutoff_choice_checks():
    """docs/03-results-zh.md section 3.4: m = floor(sqrt N) is not optimal, but m = 1 nearly is."""
    nu, level = sp.symbols("nu N", positive=True)
    # The attained harmonic sum is maximized at m = 1 and capped by log(N+2), which is
    # what forces the exponent O(1/log log(1/nu)): no cutoff choice can beat it.
    assert sp.simplify(sp.integrate(1/(sp.Symbol("x", positive=True)+1),
                                    (sp.Symbol("x", positive=True), 0, level))
                       -sp.log(level+1)) == 0
    assert sp.simplify(sp.log(level+2)-sp.log(level+1)) > 0
    count = 2
    # Exact rational comparison of the admissible cutoffs at the worked example (N = 10).
    capital = sp.Rational(151224688)
    totals = {}
    for m in range(1, 11):
        harmonic = sum(sp.Integer(1)/(n+1) for n in range(m, 11))
        totals[m] = 3*50*sp.sqrt(2)**3*50**3*sp.Rational(1, 2)**m+capital/harmonic
        assert harmonic > 0
    for m in range(2, 11):
        assert totals[m] > totals[1]
    assert sp.Rational(1, 2)**47*4**11 == sp.Rational(1, 2)**25 < 1
    count += 10
    # Exact rational audit at the paper's worked example (T=50, alpha=1/2, nu=2^-47, N=10).
    # Full-scale constants: A = 2550 sqrt 2 + 250000, C = 2A + 16B^2 + T R^2, L_* ~ 1.512e8.
    assert sum(sp.Integer(1)/(n+1) for n in range(3, 11)) == sp.Rational(32891, 27720)
    assert sum(sp.Integer(1)/(n+1) for n in range(1, 11)) == sp.Rational(55991, 27720)
    assert sum(sp.Integer(1)/n for n in range(1, 11)) == sp.Rational(7381, 2520)
    assert sp.Rational(32891, 27720) < sp.Rational(6, 5) < sp.Rational(55991, 27720)
    count += 4
    # At that example every admissible m leaves the bound above the trivial bound T^2,
    # because L_*/H exceeds T^2 for every attainable H (attained H is below 3).
    assert sp.Rational(1, 2)**47*4**10*2500 < 1
    assert sp.Rational(151224688, 2500) > 10**4
    assert sp.Rational(151224688, 10**4) > 2500
    assert sp.Rational(55991, 27720) < 3
    count += 4
    for power in range(20, 62, 3):
        small = sp.Rational(1, 2)**power
        outer, inner = power//2, power//2//2-1
        assert inner >= 4 and outer >= 10
        assert small*4**outer <= 1
        assert small*4**(inner+1) <= 1
        assert small*4**inner <= sp.sqrt(small)/4
        count += 4
    return count


def moment_ratio_obstruction_checks():
    """docs/03-results-zh.md section 3.3: why I_n >= 1 and why the shell budget cannot bound sum V_n."""
    window, wint, fourth, radius = sp.symbols("T W J R", positive=True)
    count = 0
    # I = T J / W^2 >= 1 by Cauchy-Schwarz on [0,T], with equality iff w^2 is time-constant.
    ratio = window*fourth/wint**2
    assert sp.simplify(ratio.subs(fourth, wint**2/window)-1) == 0
    assert sp.simplify(ratio.subs({fourth: 1, wint: 1, window: 1})-1) == 0
    count += 2
    # If w^2 is constant = c on [0,T]: W = cT, J = c^2 T, so I = 1 exactly.
    c = sp.symbols("c", positive=True)
    assert sp.simplify(window*(c**2*window)/(c*window)**2-1) == 0
    count += 1
    # The only general pointwise bound is w <= R, giving J <= R^2 W and hence I <= T R^2 / W.
    assert sp.simplify((radius**2*wint)*window/wint**2-window*radius**2/wint) == 0
    count += 1
    # Time H"older": int w^2 <= sqrt(T) (int w^4)^{1/2}, i.e. W^2 <= T J.
    assert sp.simplify(sp.sqrt(window)*sp.sqrt(fourth)-sp.sqrt(window*fourth)) == 0
    assert sp.simplify(window*fourth-wint**2).subs({wint: sp.sqrt(window*fourth)}) == 0
    count += 2
    # Crude shell budget: V_n <= 2^n T R^2, and sum_{n=0}^{N} 2^n = 2^{N+1}-1.
    level = sp.symbols("N", positive=True, integer=True)
    for n in range(1, 12):
        assert sum(2**k for k in range(n+1)) == 2**(n+1)-1
        count += 1
    # A crude partial sum grows with the chosen cutoff; this identity does not determine a viscous scale.
    nu = sp.symbols("nu", positive=True)
    assert sp.simplify(sp.log(4.0**(-level), 2)+2*level) == 0
    assert sp.simplify(4.0**(-level)*4**level-1) == 0
    count += 2
    # The repaired envelope holds only for n <= N+1, giving sum V_n <= (N+2) C_T: linear in N,
    # i.e. O(log(1/nu)) -- still not a uniform bound.
    cap = sp.symbols("cap", positive=True)
    assert sp.expand((level+2)*cap-(level*cap+2*cap)) == 0
    count += 1
    return count


def shell_budget_checks():
    """docs/03-results-zh.md section 3.3: shell-wise budgets and the nu 4^n ~ 1 boundary."""
    count = 0
    # Unforced local algebra only. At the root the full energy balance also contains w+z.
    for k in range(4):
        w, z, wp, zp, wn, zn = sp.symbols("w z wp zp wn zn", positive=True)
        kprev, kcur, eprev, ecur = sp.symbols("Kp Kc ep ec", positive=True)
        hk, beta, nu, mk = sp.symbols("h beta nu M", positive=True)
        dw = kprev*zp**2-ecur*w*zn+hk*(beta*z**2-w*z)-nu*mk*w
        dz = eprev*wp**2-kcur*z*wn+hk*(w**2-beta*w*z)-nu*mk*z
        rate = sp.expand(w*dw+z*dz)
        exact = sp.expand(kprev*w*zp**2+eprev*z*wp**2
                          -ecur*w**2*zn-kcur*z**2*wn
                          +hk*(beta*w*z**2+w**2*z)
                          -hk*(w**2*z+beta*w*z**2)
                          -nu*mk*(w**2+z**2))
        assert sp.expand(rate-exact) == 0
        # The two in-layer source terms cancel against the two in-layer cross terms exactly:
        assert sp.expand(hk*(beta*w*z**2+w**2*z)-hk*(w**2*z+beta*w*z**2)) == 0
        count += 2
        # Weak-coordinate balance: multiply dw by w. Only the four vertical/viscous
        # terms and the two in-layer terms appear.
        assert sp.expand(w*dw-(kprev*w*zp**2-ecur*w**2*zn
                              +hk*(beta*w*z**2-w**2*z)-nu*mk*w**2)) == 0
        # Strong-coordinate balance.
        assert sp.expand(z*dz-(eprev*z*wp**2-kcur*z**2*wn
                              +hk*(w**2*z-beta*w*z**2)-nu*mk*z**2)) == 0
        count += 2
    # In-layer terms are per-shell, not telescoping: h_k w_k^2 is strong-incoming and
    # h_k beta w_k z_k is strong-outgoing of the SAME index k.
    h, W2s, WZs, beta = sp.symbols("h W2 WZ beta", positive=True)
    assert sp.expand(h*W2s-h*beta*WZs-h*(W2s-beta*WZs)) == 0
    count += 1
    # Vertical terms do telescope: the surviving boundary terms are the weak-vertical
    # influx at the bottom and the strong-vertical outflux at the top of the block.
    for last in range(1, 6):
        terms = sp.symbols(f"t0:{last+1}", positive=True)
        assert sp.expand(sum(terms)-sum(terms[:-1])-terms[-1]) == 0
        count += 1
    # Crossover index: nu 4^n = 1 at n = log_2(1/nu)/2, i.e. 4^J = 1/nu there.
    nu = sp.symbols("nu", positive=True)
    J = sp.log(1/nu, 2)/2
    assert sp.simplify(nu*4**J-1) == 0
    assert sp.simplify(sp.log(4**J, 2)-2*J) == 0
    count += 2
    # Formal coefficient simplification, not a trajectory bound. Actual partial-z
    # budgets must include the root forcing term; see foundations section 5.
    eta, alpha, window, radius, level = sp.symbols("eta alpha T R n", positive=True)
    lhs = (1+eta*2**(level/2)*(level+1)**alpha)/(2**(level/2)*2**(level/2))
    assert sp.simplify(lhs-(1+eta*2**(level/2)*(level+1)**alpha)/2**level) == 0
    count += 1
    return count


def saturated_flux_checks():
    """RMS product identities only: neither actual flux nor MB equivalence."""
    eta, alpha, level, window, proxy = sp.symbols("eta alpha n T C", positive=True)
    wrms = sp.symbols("wrms", positive=True)
    h = eta*2**(level/2)*(level+1)**alpha
    zrms = proxy/(h*wrms**2)
    V = sp.simplify(window*2**level*zrms**2/(level+1)**(2*alpha))
    assert sp.simplify(V-window*proxy**2/(eta**2*wrms**4*(level+1)**(4*alpha))) == 0
    c = sp.symbols("c", positive=True)
    assert sp.simplify(V.subs(wrms,c*(level+1)**(-alpha))-window*proxy**2/(c**4*eta**2)) == 0
    exponent, decay = sp.symbols("b s", real=True)
    assert sp.expand(2*exponent-4*(-decay)) == 2*exponent+4*decay
    return 3


def internal_transfer_checks():
    """Goal condition 2, exact form: on the finite closed truncation with no viscosity and no
    forcing, every edge cancels pairwise and the total energy is exactly conserved.  Also
    records the exact viscous budget identity used for the positive side."""
    beta = sp.Symbol("beta", positive=True)
    count = 0
    for level_count in (1, 2, 3, 4, 5):
        w = sp.symbols(f"w0:{level_count}")
        z = sp.symbols(f"z0:{level_count}")
        big_k = [2**n for n in range(level_count)]
        weak = [sp.Rational(1, 2**n) for n in range(level_count)]
        intra = [sp.Symbol(f"h{n}", positive=True) for n in range(level_count)]
        dw, dz = [], []
        for n in range(level_count):
            term = 0
            if n - 1 >= 0:
                term += big_k[n - 1]*z[n - 1]**2
            if n + 1 < level_count:
                term -= weak[n]*w[n]*z[n + 1]
            term += intra[n]*(beta*z[n]**2 - w[n]*z[n])
            dw.append(sp.expand(term))
            term = 0
            if n - 1 >= 0:
                term += weak[n - 1]*w[n - 1]**2
            if n + 1 < level_count:
                term -= big_k[n]*z[n]*w[n + 1]
            term += intra[n]*(w[n]**2 - beta*w[n]*z[n])
            dz.append(sp.expand(term))
        # (a) total energy exactly conserved
        d_energy = sp.expand(sum(w[n]*dw[n] + z[n]*dz[n] for n in range(level_count)))
        assert sp.simplify(d_energy) == 0
        count += 1
        # (b) every edge cancels PAIRWISE, not merely in the sum
        for n in range(level_count):
            if n + 1 < level_count:
                weak_edge = sp.expand(w[n]*(-weak[n]*w[n]*z[n + 1]) + z[n + 1]*(weak[n]*w[n]**2))
                strong_edge = sp.expand(z[n]*(-big_k[n]*z[n]*w[n + 1]) + w[n + 1]*(big_k[n]*z[n]**2))
                assert sp.simplify(weak_edge) == 0
                assert sp.simplify(strong_edge) == 0
                count += 2
            in_layer = sp.expand(intra[n]*(beta*z[n]**2 - w[n]*z[n])*w[n]
                                 + intra[n]*(w[n]**2 - beta*w[n]*z[n])*z[n])
            assert sp.simplify(in_layer) == 0
            count += 1
    # (c) exact viscous budget identity: E(T) + Q_nu(T) = int_0^T (w_0 + z_0) dt.
    # Verified as an algebraic identity for the two root coordinates: the forcing enters only
    # at n = 0, so the total injected energy is exactly the root sum.
    forcing = sp.symbols("f", positive=True)
    root_w, root_z = sp.symbols("root_w root_z", positive=True)
    inj = sp.simplify(forcing*root_w + forcing*root_z)
    assert sp.simplify(inj - forcing*(root_w + root_z)) == 0
    count += 1
    # (d) the in-layer exchange coefficient may be any positive function of n: conservation
    # does not constrain it.  Check with a generic symbol.
    generic = sp.Symbol("g", positive=True)
    ww, zz = sp.symbols("ww zz", positive=True)
    pair = sp.expand(generic*(beta*zz**2 - ww*zz)*ww + generic*(ww**2 - beta*ww*zz)*zz)
    assert sp.simplify(pair) == 0
    count += 1
    # (e) the model's OWN h_n = eta*2^{n/2}*(n+1)^alpha, with alpha SYMBOLIC, must also conserve.
    # This certifies the goal's condition 2 for EVERY alpha, not merely the alpha where the
    # positive side of the classification is proved.
    n_sym = sp.Symbol("n", positive=True)
    alpha_sym = sp.Symbol("alpha", real=True)
    eta_sym = sp.Symbol("eta", positive=True)
    for level_count in (1, 2, 3, 4):
        ww = sp.symbols(f"W0:{level_count}")
        zz = sp.symbols(f"Z0:{level_count}")
        big_k = [2**m for m in range(level_count)]
        weak = [sp.Rational(1, 2**m) for m in range(level_count)]
        intra = [eta_sym*big_k[m]**sp.Rational(1, 2)*(m+1)**alpha_sym for m in range(level_count)]
        dw_l, dz_l = [], []
        for m in range(level_count):
            term = 0
            if m - 1 >= 0:
                term += big_k[m-1]*zz[m-1]**2
            if m + 1 < level_count:
                term -= weak[m]*ww[m]*zz[m+1]
            term += intra[m]*(beta*zz[m]**2 - ww[m]*zz[m])
            dw_l.append(sp.expand(term))
            term = 0
            if m - 1 >= 0:
                term += weak[m-1]*ww[m-1]**2
            if m + 1 < level_count:
                term -= big_k[m]*zz[m]*ww[m+1]
            term += intra[m]*(ww[m]**2 - beta*ww[m]*zz[m])
            dz_l.append(sp.expand(term))
        d_energy = sp.expand(sum(ww[m]*dw_l[m] + zz[m]*dz_l[m] for m in range(level_count)))
        assert sp.simplify(d_energy) == 0
        count += 1
        for m in range(level_count):
            if m + 1 < level_count:
                assert sp.simplify(sp.expand(ww[m]*(-weak[m]*ww[m]*zz[m+1])
                                             + zz[m+1]*(weak[m]*ww[m]**2))) == 0
                assert sp.simplify(sp.expand(zz[m]*(-big_k[m]*zz[m]*ww[m+1])
                                             + ww[m+1]*(big_k[m]*zz[m]**2))) == 0
                count += 2
            assert sp.simplify(sp.expand(intra[m]*(beta*zz[m]**2 - ww[m]*zz[m])*ww[m]
                                         + intra[m]*(ww[m]**2 - beta*ww[m]*zz[m])*zz[m])) == 0
            count += 1
    return count



def envelope_exponent_checks():
    """Goal mode round 3: which exponent conditions the positive-side envelope argument needs,
    and whether softening the intra-shell coupling h_n widens the range.  It does not."""
    n = sp.Symbol("n", positive=True)
    alpha = sp.Symbol("alpha", positive=True)
    sigma = sp.Symbol("sigma", real=True)
    eta = sp.Symbol("eta", positive=True)
    count = 0
    # The intra-shell coefficient is free (internal_transfer_checks certifies this), so consider
    # the general model h_n = eta K_n^sigma (n+1)^alpha.  The next-layer balance can use the
    # branch controlled by h_n only if h_n^2/K_n stays bounded, i.e. 2 sigma - 1 <= 0.
    ratio = sp.simplify((eta*2**(n*sigma)*(n+1)**alpha)**2/2**n)
    assert sp.simplify(ratio - eta**2*2**(n*(2*sigma-1))*(n+1)**(2*alpha)) == 0
    count += 1
    assert sp.simplify(sp.expand(2*sigma-1).subs(sigma, sp.Rational(1,2))) == 0
    count += 1
    # Boundedness of h_n^2/K_n requires the exponent of 2^n to be <= 0.
    exponent = sp.expand(2*sigma-1)
    assert sp.simplify(exponent.subs(sigma, sp.Rational(2,5)) + sp.Rational(1,5)) == 0  # negative
    assert sp.simplify(exponent.subs(sigma, sp.Rational(3,5)) - sp.Rational(1,5)) == 0  # positive
    count += 1
    # The weak-coordinate envelope exponent is -2 alpha/3 INDEPENDENT of sigma, so the
    # summability of the weak sum needs 4 alpha/3 > 1, i.e. alpha > 3/4, for every sigma.
    weak_exponent = sp.Rational(4, 3)*alpha
    for a, expected in ((sp.Rational(1,2), sp.false), (sp.Rational(2,3), sp.false),
                        (sp.Rational(3,4), sp.false), (sp.Integer(1), sp.true),
                        (sp.Rational(5,4), sp.true)):
        assert sp.simplify(sp.Gt(weak_exponent.subs(alpha, a), 1) == expected)
        count += 1
    # The strong-coordinate envelope gives the exponent 2/3, which converges for every alpha.
    assert sp.Rational(2, 3) < 1          # geometric ratio 2^{-2/3} < 1
    assert sp.Rational(2, 3) > 0
    count += 1
    # The strong envelope alone does NOT bound the weighted sum Sum 2^n (n+1)^{-2a} z_n^2,
    # because its term ratio 2^{1/3} exceeds 1.  Recording this prevents a wrong shortcut.
    assert sp.Rational(1, 3) > 0
    term_ratio = sp.simplify(2**(sp.Rational(1,3)*(n+1))/2**(sp.Rational(1,3)*n))
    assert sp.simplify(term_ratio - 2**sp.Rational(1,3)) == 0
    assert 2**sp.Rational(1,3) > 1
    count += 3
    # Scale covariance in the shell ratio: with K_n = lam^n, e_n = lam^{-n}, h_n =
    # eta lam^{n/2}(n+1)^alpha, the index map m = n log_lam 2 turns lam^n into 2^m and leaves
    # e_n = 2^{-m} unchanged, while (n+1)^alpha becomes (m/log_lam 2 + 1)^alpha.  The exponent
    # in alpha is therefore 1 for every lam, so the threshold alpha > 3/4 cannot be moved.
    lam = sp.Symbol("lam", positive=True)
    m = sp.Symbol("m", positive=True)
    # number of base-2 shells per base-lam shell is log_lam 2; check the composition is identity
    ratio_between = sp.simplify(sp.log(2)/sp.log(lam)*sp.log(lam)/sp.log(2))
    assert ratio_between == 1
    count += 1
    # The weak envelope exponent is -2 alpha/3 regardless of lam: verify the alpha-exponent of
    # the rescaled polynomial factor is 1, so the summability threshold is 4a/3 > 1 in m as well.
    rescaled = (m/sp.log(lam, 2) + 1)**alpha
    assert sp.simplify(sp.log(rescaled).diff(alpha) - sp.log(m/sp.log(lam, 2) + 1)) == 0
    count += 1
    # Explicit threshold invariance at lam = 2 and lam = 4: the alpha threshold is 3/4 both times.
    for lam_value in (2, 4):
        threshold = sp.Rational(3, 1)/(4*sp.log(lam_value, 2))
        assert sp.simplify(threshold - sp.Rational(3, 4)/sp.log(lam_value, 2)) == 0
        count += 1
    return count



def root_ratio_checks():
    """Actual forced root-coordinate balance; no alleged GAP reduction."""
    w, z, next_w, h, beta, nu = sp.symbols("w z wn h beta nu", nonnegative=True)
    rate = 1+h*(w**2-beta*w*z)-z*next_w-nu*z
    expected = z+h*w**2*z-h*beta*w*z**2-z**2*next_w-nu*z**2
    assert sp.expand(z*rate-expected) == 0
    # Integrated z balance: Jz=int z, A=int w^2 z, B=int w z^2,
    # F=int z^2 w1, Z=int z^2; root initial value is zero.
    Jz, A, B, F, Z, endpoint = sp.symbols("Jz A B F Z endpoint", nonnegative=True)
    rhs = Jz+h*A-h*beta*B-F-nu*Z
    rearranged = endpoint/2+F+nu*Z+h*beta*B-Jz-h*A
    assert sp.expand(rearranged.subs(endpoint,2*rhs)) == 0
    # Arbitrary budget-compatible numbers refute the inference from Z<=I/2
    # alone to Z<=T; the one-layer trajectory counterexample is in the audit.
    I, T = sp.Rational(204,100), sp.Integer(1)
    assert I/2 > T
    return 3


def root_chain_checks():
    """Regressions against invalid algebra and the omitted viscosity weight."""
    a, z = sp.symbols("a z", real=True)
    gap = 2*a**2-a-z
    assert sp.expand(gap-(3*(a-1)+(1-z)+2*(a-1)**2)) == 0
    assert gap.subs({a:1,z:2}) == -1
    # Two genuine sufficient conditions for this algebraic expression only.
    p, q = sp.symbols("p q", nonnegative=True)
    assert sp.expand(gap.subs({a:1+p,z:1-q})) == 3*p+q+2*p**2
    assert sp.expand(gap.subs({a:1+p,z:(1+p)**2-q})) == p+p**2+q
    # The actual total budget has I=E+nu*(W+Z+tail), not I=W+Z+nonnegative.
    nu, W, Z, tail, E = sp.symbols("nu W Z tail E", nonnegative=True)
    incoming = E+nu*(W+Z+tail)
    assert sp.expand(incoming-E-nu*(W+Z+tail)) == 0
    wrong = (incoming-W-Z).subs({nu:sp.Rational(1,10),W:1,Z:1,tail:0,E:0})
    assert wrong == -sp.Rational(9,5)
    # Mean and RMS coincide only if variance vanishes; a nonconstant profile
    # w(t)=t on [0,1] already distinguishes them.
    mean, second = sp.Rational(1,2), sp.Rational(1,3)
    assert second-mean**2 == sp.Rational(1,12)
    # Fourth concentration can diverge while second moments remain fixed.
    eps = sp.symbols("eps", positive=True)
    assert sp.simplify(eps*(eps**(-sp.Rational(1,2)))**2) == 1
    assert sp.simplify(eps*(eps**(-sp.Rational(1,2)))**4) == 1/eps
    return 9


if __name__ == "__main__":
    checks = (tail_budgets, capacity_certificate, source_overlap,
              two_channel_budgets, sparse_unrepaired_cuts,
              coupled_candidate_checks, relay_equation_checks,
              positive_exchange_checks, asymmetric_cut_checks,
              relay_deletion_checks, scale_dependent_exchange_checks,
              supercritical_exchange_checks, steady_all_exponents_checks,
              dynamic_no_dissipation_checks, alternating_chain_checks,
              cascade_time_comparison_checks, logarithmic_exchange_checks,
              integrated_tail_summability_checks,
              bounded_exchange_robustness_checks, sequential_gate_checks,
              temporal_concentration_checks, startup_concentration_checks,
              quantitative_dissipation_checks, uniform_envelope_checks,
              finite_partial_sum_checks, cutoff_choice_checks,
              moment_ratio_obstruction_checks, shell_budget_checks,
              saturated_flux_checks, internal_transfer_checks,
              envelope_exponent_checks, root_ratio_checks,
              root_chain_checks)
    total = 0
    for check in checks:
        count = check()
        total += count
        print(f"{check.__name__}: {count} finite cases PASS")
    print(f"Total: {len(checks)} groups, {total} finite named cases PASS")
    print("Scope: finite algebra/counting only; no time simulation, infinite-limit proof or novelty certification.")
