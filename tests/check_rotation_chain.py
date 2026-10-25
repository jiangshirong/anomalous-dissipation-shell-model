#!/usr/bin/env python3
"""机制附录 B.10 的 R 模型：有限代数核验与严格三层轨道证书。

默认运行两部分；--algebra-only 仅运行 Fraction 有限恒等式检查，
--certificate-only 仅运行严格整数区间 Taylor 轨道认证。

代数检查覆盖原始内边的双线性能量配对、三层封闭方程与末层通量，
并以 chi=3/5,4/5,7/9、lambda=chi**(-3) 的有理代理核验方向消项
及内层稳态平衡。代理使用未归一 e=(1,q)，完整保留 |e|²=1+q²。
这些有限检查不是实际二进壳 chi=2**(-1/3) 的全部代数验证；真实
参数和无限层定理的证明见正文，不由有限样本测试代替。

仅使用 Python 标准库。区间端点是 2**(-256) 的整数倍，所有认证判断
均由整数或 Fraction 完成；Decimal 仅用于向外舍入后的可读输出。

证明前提（由主数学证明独立建立，不能由本脚本代替）：
  (1) 给定五维多项式 ODE 的解存在且范数守恒；
  (2) 每个中心范数 <= R=101/100 时，36 阶 Taylor 解析余项不超过
      delta = R*(C*R*h)**37/(1-C*R*h)，C=5，h=1/64；
  (3) 真轨道误差的单步传播因子不超过 A=64/59。
本脚本严格验证每步中心范数、区间运算/舍入误差，并按上述前提生成
T=2 与 T=14 处五个坐标、末层能量和回流符号的可复现证书。
本证书只涉及有限封闭无黏三层轨道，不证明无限维异常耗散或原创性。
"""

from argparse import ArgumentParser
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from time import perf_counter


BITS = 256
SCALE = 1 << BITS
ORDER = 36
H_DEN = 64
STEPS = 896
R = Fraction(101, 100)
C = 5
A = Fraction(64, 59)
ROUND_RADIUS = Fraction(1, 1 << 200)


def dot2(x, y):
    return x[0] * y[0] + x[1] * y[1]


def add2(x, y):
    return x[0] + y[0], x[1] + y[1]


def scale2(c, x):
    return c * x[0], c * x[1]


def edge_bilinear(first_donor, second_donor, second_receiver, coefficient, a):
    """B(A,V) on one edge, with <V,B(A,V)>=0.

    Only A at the donor enters this edge. J(x,y)=(-y,x).
    Setting A=V gives the original R-model quadratic interaction.
    """
    donor = scale2(-coefficient * second_receiver[0], first_donor)
    receiver = (
        coefficient * dot2(first_donor, second_donor)
        - a * coefficient * sum(first_donor) * second_receiver[1],
        a * coefficient * sum(first_donor) * second_receiver[0],
    )
    return donor, receiver


def check_algebra():
    """Finite exact-rational checks; not a symbolic or infinite-dimensional proof."""
    a = Fraction(1, 2)
    pairings = 0
    coordinate_checks = 0
    energy_checks = 0
    flux_checks = 0
    for j in range(1, 13):
        values = tuple(Fraction((j + 2) * (k + 1) % 17 - 8, j + k + 2) for k in range(9))
        first_donor = values[0:2]
        second_donor = values[2:4]
        second_receiver = values[4:6]
        coefficient = Fraction(j + 1, j + 3)
        donor, receiver = edge_bilinear(first_donor, second_donor, second_receiver, coefficient, a)
        assert dot2(second_donor, donor) + dot2(second_receiver, receiver) == 0
        pairings += 1

        r, p, q, s, t = values[4:9]
        u0, u1, u2 = (r, Fraction(0)), (p, q), (s, t)
        d0, incoming1 = edge_bilinear(u0, u0, u1, Fraction(1), a)
        outgoing1, d2 = edge_bilinear(u1, u1, u2, Fraction(2), a)
        d1 = add2(incoming1, outgoing1)
        assembled = d0 + d1 + d2
        expected = (
            -r * p, Fraction(0),
            r * r - a * r * q - 2 * p * s,
            a * r * p - 2 * q * s,
            2 * (p * p + q * q) - 2 * a * (p + q) * t,
            2 * a * (p + q) * s,
        )
        for actual, target in zip(assembled, expected):
            assert actual == target
            coordinate_checks += 1
        assert dot2(u0, d0) + dot2(u1, d1) + dot2(u2, d2) == 0
        energy_checks += 1
        assert 2 * dot2(u2, d2) == 4 * s * (p * p + q * q)
        flux_checks += 1

    cancellation_checks = 0
    profile_coordinate_checks = 0
    for chi in (Fraction(3, 5), Fraction(4, 5), Fraction(7, 9)):
        lam = chi ** -3
        q = a * chi / (1 - a * chi)
        e = (Fraction(1), q)
        je = (-q, Fraction(1))
        assert q - a * chi * (1 + q) == 0
        cancellation_checks += 1
        # e is intentionally unnormalised: retaining its norm squared is essential.
        cancelled = add2(
            add2((dot2(e, e), Fraction(0)), scale2(a * chi * sum(e), je)),
            scale2(-1, e),
        )
        for component in cancelled:
            assert component == 0
            cancellation_checks += 1
        for n in range(1, 5):
            previous = scale2(chi ** (n - 1), e)
            current = scale2(chi ** n, e)
            following = scale2(chi ** (n + 1), e)
            _, incoming = edge_bilinear(previous, previous, current, lam ** (n - 1), a)
            outgoing, _ = edge_bilinear(current, current, following, lam ** n, a)
            for component in add2(incoming, outgoing):
                assert component == 0
                profile_coordinate_checks += 1

    counts = (pairings, coordinate_checks, energy_checks, flux_checks,
              cancellation_checks, profile_coordinate_checks)
    assert counts == (12, 72, 12, 12, 9, 24)
    print("PASS: Fraction finite algebra checks =", sum(counts))
    print("  edge bilinear pairings =", pairings)
    print("  three-shell coordinate identities =", coordinate_checks)
    print("  three-shell energy identities =", energy_checks)
    print("  last-shell flux identities =", flux_checks)
    print("  rational-proxy direction cancellations =", cancellation_checks)
    print("  rational-proxy interior-profile coordinates =", profile_coordinate_checks)
    print("  proxy chi = 3/5,4/5,7/9; these are not the actual binary-shell chi=2^(-1/3).")


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[1], -x[0]


def sub(x, y):
    return x[0] - y[1], x[1] - y[0]


def mul(x, y):
    corners = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return min(corners) // SCALE, -((-max(corners)) // SCALE)


def div_pos(x, d):
    assert isinstance(d, int) and d > 0
    return x[0] // d, -((-x[1]) // d)


def times_pos(x, d):
    assert isinstance(d, int) and d > 0
    return d * x[0], d * x[1]


def convolution(x, y, j):
    lower = 0
    upper = 0
    for k in range(j + 1):
        v = mul(x[k], y[j - k])
        lower += v[0]
        upper += v[1]
    return lower, upper


def convolution_square(x, j):
    lower = 0
    upper = 0
    for k in range((j + 1) // 2):
        v = mul(x[k], x[j - k])
        lower += 2 * v[0]
        upper += 2 * v[1]
    if j % 2 == 0:
        v = mul(x[j // 2], x[j // 2])
        lower += v[0]
        upper += v[1]
    return lower, upper


def one_step(center):
    # Coefficients u_j in u(t)=sum_j u_j*t**j, not derivatives u^(j).
    r, p, q, s, t = [[(v, v)] for v in center]
    pq = [add(p[0], q[0])]
    for j in range(ORDER):
        rp = convolution(r, p, j)
        rq = convolution(r, q, j)
        rr = convolution_square(r, j)
        ps = convolution(p, s, j)
        qs = convolution(q, s, j)
        pp = convolution_square(p, j)
        qq = convolution_square(q, j)
        pqt = convolution(pq, t, j)
        pqs = convolution(pq, s, j)
        # a=1/2: F=(-rp,r²-rq/2-2ps,rp/2-2qs,
        #           2(p²+q²)-(p+q)t,(p+q)s).
        f = (
            neg(rp),
            sub(sub(rr, div_pos(rq, 2)), times_pos(ps, 2)),
            sub(div_pos(rp, 2), times_pos(qs, 2)),
            sub(times_pos(add(pp, qq), 2), pqt),
            pqs,
        )
        for series, fj in zip((r, p, q, s, t), f):
            series.append(div_pos(fj, j + 1))
        pq.append(add(p[-1], q[-1]))
    endpoint = []
    for series in (r, p, q, s, t):
        value = series[-1]
        for coefficient in reversed(series[:-1]):
            value = add(div_pos(value, H_DEN), coefficient)
        endpoint.append(value)
    next_center = tuple((lo + hi) // 2 for lo, hi in endpoint)
    max_deviation = max(max(c - lo, hi - c) for c, (lo, hi) in zip(next_center, endpoint))
    # sqrt(5)*max coordinate error <= 3*max coordinate error.
    assert 3 * max_deviation * (1 << 200) <= SCALE, "rounding radius exceeds 2^-200"
    norm_squared_numerator = sum(v * v for v in next_center)
    assert 10000 * norm_squared_numerator <= 10201 * SCALE * SCALE, "center exits radius R"
    return next_center, max_deviation, norm_squared_numerator


def square_interval(x):
    lo, hi = x
    assert lo <= hi
    if lo >= 0:
        return lo * lo, hi * hi
    if hi <= 0:
        return hi * hi, lo * lo
    return Fraction(0), max(lo * lo, hi * hi)


def format_bound(x, upper):
    # Both decimal division and quantization use the requested outward direction.
    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        y = Decimal(x.numerator) / Decimal(x.denominator)
        return format(y.quantize(Decimal("0.000000000001")), "f")


def format_interval(x):
    return "[" + format_bound(x[0], False) + ", " + format_bound(x[1], True) + "]"


def check_certificate():
    started = perf_counter()
    h = Fraction(1, H_DEN)
    local_tail = R * (C * R * h) ** (ORDER + 1) / (1 - C * R * h)
    center = (SCALE, 0, 0, 0, 0)
    max_round_units = 0
    max_center_norm_sq_num = SCALE * SCALE
    records = {}
    for step in range(1, STEPS + 1):
        center, round_units, norm_sq_num = one_step(center)
        max_round_units = max(max_round_units, round_units)
        max_center_norm_sq_num = max(max_center_norm_sq_num, norm_sq_num)
        if step in (128, STEPS):
            radius = (local_tail + ROUND_RADIUS) * (A ** step - 1) / (A - 1)
            intervals = tuple((Fraction(c, SCALE) - radius, Fraction(c, SCALE) + radius) for c in center)
            last_layer = tuple(sum(x[i] for x in map(square_interval, intervals[3:5])) for i in (0, 1))
            middle_layer = tuple(sum(x[i] for x in map(square_interval, intervals[1:3])) for i in (0, 1))
            records[step] = intervals, radius, last_layer, middle_layer
            print("T =", Fraction(step, H_DEN), "steps =", step, flush=True)
            print("global error radius <=", format_bound(radius, True), flush=True)
            for name, value in zip(("r", "p", "q", "s", "t"), intervals):
                print(name, format_interval(value), flush=True)
            print("s^2+t^2", format_interval(last_layer), flush=True)
            print("p^2+q^2", format_interval(middle_layer), flush=True)
            # Emit exact centers and exact radius ingredients for independent reuse.
            print("exact center numerators over 2^256:", center, flush=True)
            print("exact error formula: (delta+2^-200)*((64/59)^%d-1)/(5/59)" % step, flush=True)
    at2 = records[128]
    at14 = records[STEPS]
    assert at2[1] < Fraction(1, 10 ** 34)
    assert at14[1] < Fraction(91, 10 ** 10)
    assert at2[2][0] > Fraction(1, 2), "no certified half-energy transfer at T=2"
    assert at14[0][3][1] < 0, "no certified negative receiver s at T=14"
    assert at14[3][0] > 0, "no certified nonzero middle-shell energy at T=14"
    # For this closed system d(s²+t²)/dt = 4*s*(p²+q²).
    s_bounds = at14[0][3]
    middle_bounds = at14[3]
    flux = (4 * s_bounds[0] * middle_bounds[1], 4 * s_bounds[1] * middle_bounds[0])
    assert flux[1] < 0
    print("T=14 derivative of last-shell norm squared:", format_interval(flux))
    print("all 896 center norms squared <=", format_bound(Fraction(max_center_norm_sq_num, SCALE * SCALE), True))
    print("largest endpoint rounding bound 3*maxdev/2^256: numerator =", 3 * max_round_units)
    print("PASS: exact rational assertions; T=2 last-layer norm squared > 1/2; T=14 strict backflow.")
    print("elapsed_seconds = %.3f" % (perf_counter() - started))


def main():
    parser = ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--algebra-only", action="store_true", help="仅运行有限 Fraction 恒等式检查")
    selection.add_argument("--certificate-only", action="store_true", help="仅运行三层严格轨道证书")
    args = parser.parse_args()
    if not args.certificate_only:
        check_algebra()
    if not args.algebra_only:
        check_certificate()


if __name__ == "__main__":
    main()
