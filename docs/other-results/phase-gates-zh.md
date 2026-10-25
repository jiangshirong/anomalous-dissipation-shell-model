# 自治相位门 A：同控制器能量下的快慢耗散对照

快侧 Ω_n=1024(n+1)²，慢侧 Ω_n=K_n^(−2)。两侧载荷、控制器初值及标准黏性相同；转子初始半径 εK_n^(−2)，转子能量为 ε²K_n^(−4)/2。快侧固定窗实际载荷损失有正下界，慢侧每个固定窗总损失趋零。原 Ω_n=1 模型仍未解决。

<!-- source: docs/appendices/mechanisms-zh.md#b12-来流驱动的旋转门同控制器能量下的快慢耗散对照 -->
## B.12. 来流驱动的旋转门：同控制器能量下的快慢耗散对照

以下给出 A 的明确快慢变体的完整内部证明。角速度由当前振幅驱动，初始门全部关闭；控制器计入能量，其零黏损失趋零，快侧正下界来自输运振幅的实际黏性损失。慢侧每个固定窗总损失趋零，而两侧每层控制器能量在所有时间完全相同。[原单位角速度及记忆门](memory-gates-zh.md)的全链穿越问题仍未解决。自治旋转、快旋转平均及修正能量已有 Tao 先例，不能将这些工具称为首次提出；具体覆盖边界见本文件末节。

<!-- source: docs/appendices/mechanisms-zh.md#b121-模型解类与定理 -->
### B.12.1 模型、解类与定理

令 K_n=2^n、Ω_n=Ω(n+1)²、Ω=1024，n≥0，固定 ε>0。每层有非负输运振幅 u_n 和实转子 (p_n,q_n)，定义

```math
\begin{aligned}
g_n&=\frac{q_n^2}{p_n^2+q_n^2},\\
u_n'&=K_{n-1}g_{n-1}u_{n-1}^2
-K_ng_nu_nu_{n+1}-\nu K_n^2u_n,\\
p_n'&=-\Omega_nK_nu_nq_n-\nu K_n^2p_n,\qquad
q_n'=\Omega_nK_nu_np_n-\nu K_n^2q_n.
\end{aligned}
```

删除负下标项，无外力；共同初值为 u_0=1、u_n=0（n≥1）、p_n=εK_n^(−2)、q_n=0。全部参数与初值独立于 ν。封闭截断保留 0,…,J 层，删除跨 J 的完整输运边，保留本层转子，不引入末端吸收。逐边输运与转子旋转分别保能量；非黏部分二次齐次，但含归一化比值，并非多项式二次壳模型。

转子半径为 εK_n^(−2)exp(−νK_n²t)，任意有限时刻严格正。写成相位形式，

```math
g_n=\sin^2\theta_n,\qquad
\theta_n'=\Omega_nK_nu_n,\qquad \theta_n(0)=0.
```

相位不是预排时钟；u_n=0 时本层相位停止。定义 H=ℓ²、||u||_V²=ΣK_n²u_n²，以及载荷的实际损失

```math
Q_{\nu,u}(T)=\nu\int_0^T\sum_{n\ge0}K_n^2u_n^2\,dt.
```

**定理。** 每个 ν>0 存在唯一指定初值的非负能量类解；u 及 (p,q) 属于 C([0,∞);ℓ²)∩L²_loc V，并有

```math
\frac12\|u(t)\|_H^2+Q_{\nu,u}(t)=\frac12.
```

置 χ=2^(−1/3)，定义明确常数

```math
\begin{gathered}
h=\chi(1-\chi^2)^2(1-\chi^4),\qquad c_A=h/32,\qquad
T_A=4/c_A,\\
A_0=\chi^{-4},\qquad B_0^2=(1-\chi^4)^{-1},\\
\varepsilon_A=
\min\left\{\frac1{4A_0},\frac1{64B_0^2T_A}\right\}>0 .
\end{gathered}
```

则所有 0<ν≤1 都满足 Q_(ν,u)(T_A)>ε_A。转子损失另满足

```math
Q_{\nu,\mathrm{rot}}(T)
=\frac{\varepsilon^2}{2}\sum_{n\ge0}K_n^{-4}
(1-e^{-2\nu K_n^2T})
\le\frac43\varepsilon^2\nu T\longrightarrow0 .
```

总能量预算包括两种损失，初始总能量为 1/2+8ε²/15。因此正异常耗散不是由初始控制器能量消失制造。未宣称上述常数最优。

**慢侧对照。** 只把 Ω_n 改为 K_n^(−2)，其余图、作用公式、初值、ε 和标准黏性完全相同，则每个固定 T 的总实际损失趋零。两侧所有 ν>0、t≥0、n≥0 的转子能量均为 ε²K_n^(−4)e^(−2νK_n²t)/2，转子累计损失也完全相同。这个确定性自治对照改变的是状态驱动的相位速度；不是瞬时门系数同分布的 D1 对照，也未证明临界角速度。

<!-- source: docs/appendices/mechanisms-zh.md#b122-固定正黏性构造唯一性与控制器能量 -->
### B.12.2 固定正黏性构造、唯一性与控制器能量

有限系统的正性及载荷预算给 0≤u_n≤1。高层初值为零，舍弃向外损失后用 Duhamel，依次得到

```math
u_n(t)\le\frac1{2\nu K_n}\quad(n\ge1),\qquad
u_n(t)\le\frac1{2\nu^3K_n^3}\quad(n\ge2).
```

所以 u_n≤C_νK_n^(−3)，C_ν=max{8,(2ν³)^(−1)}，对时间及截断一致。这个估计允许依赖 ν，仅用于固定黏性的构造。各坐标及相位在固定有限窗口紧；上述尾界给 C_tH、C_tV 的强收敛，还控制 ΣK_n⁴u_n² 的尾。因此可传递方程、初值及精确载荷预算，无需预设黏性一致正则性。跨切面通量也由该尾界趋零。

唯一性在相位的加权空间中闭合。对两解 (u,θ)、(v,φ)，令 ξ_n=(θ_n−φ_n)/(Ω_nK_n)，则 ξ'=u−v。固定相位时，振幅非线性的 H 差估计由 sup K_nu_n、sup K_nv_n 有界控制。sin² 是 1-Lipschitz，相位差产生的系数由

```math
\sup_n\Omega_nK_n^2v_n^2,\qquad
\sup_n\Omega_nK_n^2v_nv_{n+1}
\le C_\nu^2\sup_n\Omega_nK_n^{-4}<\infty
```

控制。有限项移位不改变 ℓ² 范数，故有常数 L_ν,M_ν<∞ 使

```math
\frac12\frac d{dt}\|u-v\|_H^2
+\nu\|u-v\|_V^2
\le L_\nu\|u-v\|_H^2
+M_\nu\|u-v\|_H\|\xi\|_H.
```

合并 ξ'=u−v 后用 Gronwall 得唯一性。任意指定初值的非负能量类解同样满足 Duhamel 尾界，故这不是只对某个抽取子序列的唯一性。

由相位重构 p_n=εK_n^(−2)e^(−νK_n²t)cosθ_n、q_n=εK_n^(−2)e^(−νK_n²t)sinθ_n。旋转项不做净功，半径公式逐坐标成立；求和给上述转子预算和 Q_(ν,rot) 上界。初始半径序列属于 V，但不属于黏性算子定义域 D(diag K_n²)，因此不声称转子在 t=0 已为 H 强可微；使用 C_tH∩L²_tV、导数 L²_tV′ 即可。比值在零半径处的模型奇异性并未被这一指定初值构造消除。

<!-- source: docs/appendices/mechanisms-zh.md#b123-只校正内边的精确平均恒等式 -->
### B.12.3 只校正内边的精确平均恒等式

以下在已构造的无限解上取有限矩。令

```math
\begin{gathered}
a_n=\chi^n,\quad w_n=K_n^{2/3}=K_na_n,\quad
z_n=\sqrt\chi K_n^{1/3}u_n,\qquad
M_N=\sum_{n=0}^Na_nu_n,\\
\mathcal H_N=z_0^2+\sum_{n<N}(z_{n+1}-z_n)^2,\qquad
\mathcal E_N=w_N(u_N^2+u_{N+1}^2),\qquad
V_N=\sum_{n=0}^Na_nK_n^2u_n .
\end{gathered}
```

内边的平均二次式恰为

```math
D_N=\sum_{n<N}K_na_nu_n(\chi u_n-u_{n+1})
=\frac12(\mathcal H_N-z_N^2).
```

利用相位的实际方程，

```math
(\sin2\theta_n)'=2\Omega_nK_nu_n\cos2\theta_n,
```

所以无须除以 u_n，甚至在 u_n=0 时仍有

```math
K_na_nu_n(\chi u_n-u_{n+1})(g_n-\tfrac12)
=-\frac{a_n}{4\Omega_n}
(\chi u_n-u_{n+1})(\sin2\theta_n)'.
```

只对 n<N 使用该恒等式。令

```math
C_N=\sum_{n<N}\frac{a_n}{4\Omega_n}
(\chi u_n-u_{n+1})\sin2\theta_n,\qquad
S_N=M_N+C_N .
```

则精确地

```math
S_N'=\frac14(\mathcal H_N-z_N^2)
-w_Ng_Nu_Nu_{N+1}-\nu V_N+R_N,\qquad
R_N=\sum_{n<N}\frac{a_n}{4\Omega_n}
(\chi u_n'-u_{n+1}')\sin2\theta_n .
```

跨 N 的真实边没有平均化，其贡献完整保留。校正项只含内边，因而 R_N 不引入 u_(N+2)。此外

```math
|C_N|\le\delta M_N,\qquad
\delta=\frac{\chi+\chi^{-1}}{4\Omega}<10^{-3},
\qquad S_N(0)=1.
```

因此 S_N≥0，且 S_N≤(1+δ)(1−χ²)^(−1/2)，对 N、ν、时间一致。

<!-- source: docs/appendices/mechanisms-zh.md#b124-用离散-hardy-估计吸收自驱动相位误差 -->
### B.12.4 用离散 Hardy 估计吸收自驱动相位误差

置 b_n=a_n/(n+1)²、U_N=Σ_(n≤N)w_nu_n²/(n+1)²、L=χ+4/χ<7。相邻移位权比较给

```math
|R_N|\le\frac{L}{4\Omega}\sum_{n=0}^Nb_n|u_n'|.
```

入流贡献至多 χU_N。对内边出流，令 A_n=√w_n u_n/(n+1)，则

```math
\sum_{n<N}b_nK_nu_nu_{n+1}
=\sum_{n<N}\chi\frac{n+2}{n+1}A_nA_{n+1}
\le2\chi U_N.
```

末边至多上述切面量的一半，黏性贡献至多 νV_N。因此

```math
|R_N|\le\frac{L}{4\Omega}
\left(3\chi U_N+\frac12\mathcal E_N+\nu V_N\right).
```

有限 Hardy 不等式给

```math
\sum_{n=0}^N\frac{z_n^2}{(n+1)^2}
\le4\left[z_0^2+\sum_{n<N}(z_{n+1}-z_n)^2\right].
```

这里可对差分 d_0=z_0、d_j=z_j−z_(j−1) 应用累计和形式的 Hardy 不等式；其常数 4 与 N 无关。故 U_N 至多上述差分平方和的 4/χ 倍，得到

```math
|R_N|\le\frac{3L}{\Omega}\mathcal H_N
+\frac{L}{8\Omega}\mathcal E_N
+\frac{L}{4\Omega}\nu V_N.
```

用 3/4<χ<4/5、Ω=1024，有

```math
S_N'\ge
\left(\frac14-\frac{3L}{\Omega}\right)\mathcal H_N
-\left(\frac\chi4+\frac12+\frac{L}{8\Omega}\right)\mathcal E_N
-\left(1+\frac{L}{4\Omega}\right)\nu V_N
\ge\frac18\mathcal H_N-\mathcal E_N-2\nu V_N.
```

最后，由 M_N=χ^(−1/2)Σχ^(2n)z_n，将 z_n 展成从根起的差分和，再用 Cauchy–Schwarz：

```math
M_N^2\le
\frac{\mathcal H_N}{\chi(1-\chi^2)^2(1-\chi^4)}.
```

结合 S_N≤(1+δ)M_N，以 c_A=h/32 得到

```math
S_N'\ge c_AS_N^2-\mathcal E_N-2\nu V_N.
```

(n+1)² 的角速度因子用于产生 Hardy 可吸收权；不是已经证明单位角速度版本同样成立。

<!-- source: docs/appendices/mechanisms-zh.md#b125-实际黏性预算与固定窗口正下界 -->
### B.12.5 实际黏性预算与固定窗口正下界

只为估计选 N=N_ν=⌊(3/4)log₂(1/ν)⌋，模型及初值不随 ν 改变。与真实权 K_n² 相比，

```math
\int_0^t\mathcal E_N\,ds\le A_0Q_{\nu,u}(t),\qquad
\nu\int_0^tV_N\,ds\le B_0\sqrt{tQ_{\nu,u}(t)}.
```

第一式使用 χ⁴<νK_N^(4/3)≤1 和 K_(N+1)²=4K_N²；第二式是有限和 Cauchy–Schwarz 及 νΣ_(n≤N)K_n^(4/3)≤(1−χ⁴)^(−1)。没有对发散的无穷正幂权求和。

于是所有 t≤T 有

```math
S_N(t)\ge1+c_A\int_0^tS_N(s)^2\,ds
-A_0Q_{\nu,u}(T)-2B_0\sqrt{TQ_{\nu,u}(T)}.
```

若 Q_(ν,u)(T_A)≤ε_A，两误差各至多 1/4，故 S_N≥1/2+c_A∫S_N²。对应 Riccati 下解在 2/c_A 前爆破，与 S_N 的一致上界矛盾，而 T_A=4/c_A。因此 Q_(ν,u)(T_A)>ε_A，对全部 0<ν≤1 成立。N=0 时空内边和、空校正及上述末端估计均有效。

这是先构造固定正黏性的无限解，再对其有限矩应用真实预算的证明；没有先令有限截断黏性消失，也没有通过控制器耗能或人为末端吸收取得下界。该特定自治门模型的内部定理闭合，物理适用性、全局优先权和实质原创性评价仍须与具体既有定理及等价推广区分。

<!-- source: docs/appendices/mechanisms-zh.md#b126-同控制器能量的慢侧严格无异常耗散 -->
### B.12.6 同控制器能量的慢侧严格无异常耗散

取 Ω_n=K_n^(−2)>0。固定正黏性构造仍有 B.12.2 的 Duhamel 尾界。相位空间现在为 η_n=K_n²，成员性由

```math
\|\theta(t)\|_\eta
=\left\|\int_0^t u(s)\,ds\right\|_H\le t
```

给出，而不是把逐坐标上界错误地当成可和界。相位差变量 ξ_n=K_n(θ_n−φ_n) 仍满足 ξ'=u−v，系数差强迫由 v_n²、v_nv_(n+1) 控制，唯一性同样成立。

门开启速度的实际方程给

```math
0\le\theta_n(t)=K_n^{-1}\int_0^tu_n(s)\,ds\le t/K_n,
\qquad
0\le g_n(t)\le t^2/K_n^2.
```

跨 N 的真实载荷能量流因而不超过 t²/K_N。尾部 n≥N+1 的载荷初值为零，尾预算和头部能量界给

```math
Q_{\nu,u,>N}(T)\le\frac{T^3}{3K_N},\qquad
Q_{\nu,u,\le N}(T)\le\nu K_N^2T.
```

先构造固定 ν 的无限链，再固定 N 令 ν↓0、最后 N→∞，得到 Q_(ν,u)(T)→0。载荷与控制器分别估计，避免把非零的初始控制器尾当作零；完整系统有

```math
Q_{\nu,\mathrm{total}}(T)
\le\nu K_N^2T+\frac{T^3}{3K_N}
+\frac43\varepsilon^2\nu T\longrightarrow0.
```

转子半径演化不依赖 Ω_n，所以快慢两侧控制器的逐层能量及真实损失完全相同。两侧角速度系数都严格正，整个系统自治，未冻结某个门或预排开关。此对照不证明每层瞬时通量相同，也不证明控制器振幅相位的全状态分布相同；差异恰在相位动力学和由此控制的实际输运。


<!-- source: docs/04-evidence-zh.md#831-记忆门探针与自治快慢旋转门的新定理 -->
## 文献与证据范围

**单一先例的具体边界。** 已核对 [Tao，1402.0290v3，§5.1、§5.4](https://arxiv.org/html/1402.0290v3#S5.SS4)、[§5.5 Theorem 5.3 与式(5.22)、(5.34)–(5.35)](https://arxiv.org/html/1402.0290v3#S5.SS5) 及 [Theorem 4.2](https://arxiv.org/html/1402.0290v3#S4)。原文已有自治泵／转子、延迟后快速传递，用旋转恒等式及导数纠正构造修正能量；全部控制变量计入能量，而且控制变量可以很小。不能称“微能量控制器”“自治门”或“原函数平均纠正”为本项目新概念，也不能说原文漏算控制器能量。

Tao 的指定结论是特定二次电路的有限时爆破，未直接给 A 的共同数据、固定窗、消失黏性真实损失下界或同控制器逐层能量的快慢对照。当前 A 的归一化角门调节独立载荷，且另证控制器损失趋零；这是本次实际完成的证明义务，不单凭外观差别认证原创性。非多项式比值及零半径奇异性是实质模型代价。

<!-- source: docs/04-evidence-zh.md#832-a-的-bmr-覆盖边界与原创性核查收束 -->
### BMR 包络的具体适用性

已核对 [BMR 1403.2852v2 的方程(1)](https://arxiv.org/html/1403.2852v2#S1.SS1)、[Theorems 5、13 与包络定义](https://arxiv.org/html/1403.2852v2#S3)、[§3.3.1](https://arxiv.org/html/1403.2852v2#S3.SS3.SSS1)。结论明确：A 的快慢两侧都进入固定正黏性的能量包络框架，指定定理没有直接给消失黏性耗散对照。

以 j=n+1、X_j=u_(j−1)，将 BMR 的耗散序列改记 G_j，取

```math
\beta=2,\qquad k_j=4^j,\qquad G_j=4/\nu,\qquad
\phi_j(t)=2^{-j-1}\sin^2\theta_{j-1}(t).
```

则 φ_jk_j=K_(j−1)sin²θ_(j−1)，k_j/G_j=νK_(j−1)²，精确匹配载荷方程。G_j 的正性、单调性、倒数和发散及定理中的加权条件成立。初门为零不违背假设：需要正的是 G，φ 只需有界。沿已有门轨道的包络适用；§3.3.1 为联合存在唯一性另要求系数的统一 Lipschitz 条件，不能据此跳过 B.12 自己的增广系统构造。

这里的断点是黏性一致性。用共同上界 Φ=1/4，其 Definition 6 的包络递推因子为

```math
C_j(v)=\left(1+\frac{2}{G_j\Phi v}\right)^{-1}
      =\left(1+\frac{2\nu}{v}\right)^{-1}.
```

对共同根脉冲，y_1=y_2=2；固定 j 时递推给 y_j(ν)→2，故这套包络本身不提供 ν 一致的尾紧性。快慢 Ω 信息在共同系数上界中也消失。这只说明指定定理及该包络不能直接推出对照，不断言无法发展更强方法。将转子能量一并纳入逐层半径仍可得到有界径向系数，因此“计入全部控制器能量”也不是逃离 BMR 框架的理由。固定 ν 的能量等式不意味着真实黏性积分在 ν↓0 时趋零。
