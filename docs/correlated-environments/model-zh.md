# D1：相关随机系数环境的异常耗散分类

相同单边完整时间路径分布，可以对应相反的零黏实际耗散。这里用平稳独立有限块的共同重置／私有翻转环境，给出任意有限块长、正有限速率和重置比例的必要充分分类。

[模型与统一定理](model-zh.md) · [完整证明](proofs-zh.md) · [推论与边界](consequences-zh.md) · [文献及证据](literature-zh.md) · [补充方法](supplementary-methods-zh.md)

## 方程与实际黏性损失

<!-- source: docs/appendices/mechanisms-zh.md B.11.1 -->

令 $K_n=2^n$，$n\ge0$。研究无迫的标量正壳链，初值为根层单位脉冲：

```math
\dot u_n=c_{n-1}(t)u_{n-1}^2-c_n(t)u_nu_{n+1}-\nu K_n^2u_n,
\qquad u_n(0)=\delta_{n0},\qquad
c_n(t)=K_n+(K_n-1)\zeta_n(t).
```

删除负下标项；有限 Galerkin 系统删除完整末边。$c_n\in\{1,2K_n-1\}$，$c_0=1$。对每个固定 $\nu>0$，先构造无限层解，再取零黏极限。使用

```math
H=\ell^2,\quad \|u\|_V^2=\sum_{n\ge0}K_n^2u_n^2,\qquad
Q_\nu(T)=\nu\int_0^T\sum_{n\ge0}K_n^2u_n(t)^2\,dt,
\qquad \frac12\sum_{n\ge0}u_n(t)^2+Q_\nu(t)=\frac12.
```

振幅非负，属于 $C([0,T];H)\cap L^2(0,T;V)$。环境只改变系数，振幅保持连续，重置不注入能量。构造、适应性、唯一性和预算见完整证明第一节。

## 平稳块环境与保持的统计量

<!-- source: docs/appendices/mechanisms-zh.md B.11.22 -->

将非根边分为相邻、非空、全部有限的区块，选择确定参数：

```math
m_1=1,\quad m_{j+1}=m_j+\ell_j,\quad
I_j=\{m_j,\ldots,m_{j+1}-1\},\quad
\ell_j\in\mathbb N_+,\quad 0<r_j<\infty,\quad \theta_j\in[0,1].
```

参数及块划分均独立于黏性。在以下单块定义中以 $\theta=\theta_j$ 记其重置比例。

在第 j 块中：

- 共同重置以速率 θr_j 发生，每次独立公平抽取一个符号，把本块全部 ζ_n 同时设为该符号；抽到原符号的重置也保留为时钟事件。
- 每条边另有相互独立的翻转时钟，速率为 (1−θ)r_j/2，每次只把自己的 ζ_n 变号。
- 各块的全部时钟、标记及平稳初始环境独立。

这些操作只改变系数 c_n=K_n+(K_n−1)ζ_n，不重置振幅 u，也不向能量账注入能量。每个有限块的总事件率有限；整个无限环境不要求具有有限总跳率。

若 z^(n) 表示仅翻转第 n 符号，1、−1 表示整块同号向量，块生成元为

```math
\mathcal G_j f(z)=
\theta r_j\left[\frac{f(\mathbf1)+f(-\mathbf1)}2-f(z)\right]
+\frac{(1-\theta)r_j}{2}
\sum_{n\in I_j}\bigl[f(z^{(n)})-f(z)\bigr].
```

**平稳初始化。** θ>0 时，可使用双向时间的平稳标记 Poisson 重置和独立翻转时钟：在时刻 t，取最近一次共同重置的公平符号 S，随后按各私有时钟的翻转次数取奇偶。最近重置的年龄 B 服从 Exp(θr_j)，独立于 S 和重置后的各私有时钟。等价地，在 t=0 按这条年龄混合抽取块初值，再使用独立的未来时钟。θ=0 时直接从逐边独立公平符号启动。这给出了每个 θ 的明确平稳块过程，不把非平稳同步初值当成稳健性。

**单边完整路径分布不变。** 对任意只依 ζ_n 的函数 g，直接代入生成元得

```math
\mathcal G_j[g(z_n)]
=\frac{r_j}{2}\bigl[g(-z_n)-g(z_n)\bigr].
```

因此每条边单独看始终是翻转率 r_j/2 的公平平稳二态 Markov 链，其完整时间路径分布与 θ 无关；不仅是均值或自相关相同。特别地，对 n∈I_j，

```math
\mathbb E[\zeta_n(t)\zeta_n(0)]=e^{-r_j|t|}.
```

**确有持续失步。** 对同块不同边 n、k，年龄混合给

```math
\mathbb E[\zeta_n(t)\zeta_k(t)]
=\mathbb E e^{-2(1-\theta)r_jB}
=\frac{\theta}{2-\theta},\qquad
P(\zeta_n(t)\ne\zeta_k(t))=\frac{1-\theta}{2-\theta}.
```

θ=0 时由独立性得到相同公式的端点值。0<θ<1 时，每个有限块的所有符号配置均有严格正的平稳概率：条件于任意 B>0、S，各条边的奇偶概率都严格介于 0 与 1。θ=1 才退化为块内始终同步，且此时有效翻转率为 r_j/2，恢复 [同步块分类及无限尾块](supplementary-methods-zh.md#分块相关二态环境的精确分类) 的模型。

年龄变量换成 r_jB 后，其分布与 r_j 无关。因此固定 θ 和块划分时，改变 r_j 序列 不改变单时全场分布；固定 r_j 序列和块划分而改变 θ，则保持每条边的完整时间路径分布，但通常改变单时全场联合分布。

## 统一费用判据与完整分类

<!-- source: docs/appendices/mechanisms-zh.md B.11.39 -->

区块 I_j=[m_j,m_(j+1)) 相邻、全部有限且相互独立，m_1=1。现在允许任意确定参数

```math
\ell_j=m_{j+1}-m_j\in\mathbb N_+,\qquad
0<r_j<\infty,\qquad 0\le\theta_j\le1.
```

参数不依黏性或样本，不要求增长率、单调性或 θ_j 的正下界。共同公平重置率为 θ_jr_j，私有翻转率为 (1−θ_j)r_j/2；仍从各块的平稳律启动。每条边完整时间过程始终是翻转率 r_j/2 的公平二态链，与 θ_j 无关。

对 0<θ<1，定义

```math
\varepsilon=1-\theta,\qquad
\Phi_\theta(k)=k+\left\lceil\frac{4^k}{\theta(1-\theta)}\right\rceil,
\qquad
\mathcal N_j=\min\{h\ge1:\Phi_{\theta_j}^{\circ h}(m_j)\ge m_{j+1}\}.
```

定义块费用

```math
\mathcal C_j=
\begin{cases}
\ell_j/r_j,&\theta_j=0,\\[2mm]
\displaystyle\frac{\min\{\ell_j,\theta_j^{-1}\}}{r_j}
 +\frac{\mathcal N_j-1}{\theta_jr_j},&0<\theta_j<1,\\[3mm]
1/r_j,&\theta_j=1.
\end{cases}
\qquad \mathcal C=\sum_j\mathcal C_j.
```

当 0<θ_j<1 且 ℓ_j≤1/θ_j 时，必有 𝒩_j=1，费用就是 ℓ_j/r_j；当 ℓ_j>1/θ_j 时，费用恰为 𝒩_j/(θ_jr_j)。第一项负责新块入口，第二项负责巨块内部额外的传递阶段。

**定理。** 对上述全部参数：

1. 𝒞=∞ 时，几乎处处同时对全部有限 T 有 Q_ν(T)→0；固定 T 的期望也趋零。
2. 𝒞<∞ 时，对任意 a∈(0,1)，存在不依 ν 的环境时刻 τ_a，Eτ_a<∞，使共同概率一事件上全部 0<ν≤1、s≥0 满足

```math
Q_\nu(\tau_a+s)\ge\frac{1-a-e^{-2s/\nu}}2.
```

因此正侧包括确定有限窗的期望及共同正概率 AD，并有

```math
\lim_{T\to\infty}\liminf_{\nu\downarrow0}Q_\nu(T)
=\lim_{T\to\infty}\limsup_{\nu\downarrow0}Q_\nu(T)=\frac12
\quad\text{几乎处处}.
```

这是指定平稳块独立、共同精确重置／私有翻转环境族的完整分类。概率一结论中的极限顺序为先零黏、后长时间；固定预定窗口的正侧结论是期望和共同正概率下界。
