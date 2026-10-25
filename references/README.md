# 文献入口与使用定位

本目录提供原始来源和已记录的读取定位。第三方论文全文仅在本地保存；具体适用性与贡献比较集中在相应主题文档。

- [确定性双通道成果的文献比较](../docs/sparse-modification/literature-zh.md)
- [D1 的文献比较与定理适用性](../docs/correlated-environments/literature-zh.md)
- [其他模型的有效推导](../docs/other-results/earlier-models-zh.md)
- [关键纠错和证据边界](../docs/history/corrections-zh.md)

## 双通道网络及基础方法

| 原文 | 用途与已记录定位 |
|---|---|
| [Cheskidov–Friedlander，The vanishing viscosity limit for a dyadic model](https://arxiv.org/html/0810.3718) | §2、Theorems 2.2、3.4、4.2；已核对 Eqs. (2.1)、(2.3)、(2.4) 与 2.2 证明；几何单调性使用 c>3/2、相邻归一化阻尼比小于 2，不能直接换为本链交错序列 |
| [Keller–Pinchover–Pogorzelski，Criticality theory for Schrödinger operators on graphs](https://arxiv.org/abs/1708.09664) | 2017 预印本，[2020 期刊版](https://doi.org/10.4171/JST/286)；已读取[作者上传稿](https://www.researchgate.net/publication/319415923_Criticality_theory_for_Schrodinger_operators_on_graphs) §2.1–2.2、§4.2 Proposition 4.8 及证明；[双通道模型的图基态映射](../docs/sparse-modification/model-zh.md#公式来源已知图基态表示的直接应用)排除 (R) 代数的新方法主张。arXiv PDF 未成功读取；原文紧支撑静态表示不替代本项目无限维轨道分析 |
| [Friedlander–Glatt-Holtz–Vicol，Inviscid Limits for a Stochastically Forced Shell Model of Turbulent Flow](https://arxiv.org/html/1404.1098v1) | Eqs. (1.1)–(1.2)、Proposition 4.1、Theorems 4.2、5.1、6.1 及 §6 证明；随机统计稳态与确定性平衡点区分，长期耗散定理要求 1≤c<2 |
| [Barbato–Bianchi–Flandoli–Morandin，A dyadic model on a tree](https://arxiv.org/abs/1207.2846) | 一父树、Prop. 4.1、Theorem 6.4；作者稿 PDF 第 21–23 页曾目视核对 |
| [Bianchi–Morandin，Structure function and fractal dissipation for an intermittent inviscid dyadic model](https://arxiv.org/html/1607.01334v2) | Eq. (3)、Theorems 3.1–3.2；一般归一化系数有正下界，不能只看重复系数特例 |
| [Metta，Stationary solutions for dyadic mixed model of the Euler equation](https://arxiv.org/html/2008.03747v2) | Eq. (1)、§2、Theorem 2.2 及 §4；混合 KP–Obukhov 的无黏稳态分类 |
| [Jeong–Li，A Blow-Up Result for Dyadic Models of the Euler Equations](https://arxiv.org/html/1406.2423v1) | Theorem 2、Remark 4；高正则爆破与多邻居推广预期，不是 W/S 稳态对照 |
| [Tao，Finite time blowup for an averaged three-dimensional Navier–Stokes equation](https://arxiv.org/abs/1402.0290) | §5.1、Eq. (5.3)；泵单元的明确先例，不赋予本模型真实 NS 来源 |
| [Gürcan–Li–Morel，Turbulence as a Network of Fourier Modes](https://doi.org/10.3390/math8040530) | §§2、4 的二次网络表达；曾使用公开原论文镜像，不以未目视的复杂下标作证明 |
| [Gürcan，Dynamical Network Models of the Turbulent Cascade](https://arxiv.org/abs/2004.01134) | v2 §§II–III；网络重连／加边及数值级联，不是退化正性图的稳态分类 |
| [Plunian–Stepanov–Frick，Shell Models of Magnetohydrodynamic Turbulence](https://arxiv.org/pdf/1209.3844v2) | v2 PDF 第 40 页 Eq. (233)、第 68 页 Appendix C；这些页的符号曾目视核对 |
| [Plunian–Stepanov，A non local shell model of turbulent dynamo](https://arxiv.org/html/physics/0701141v1) | §§2.1、2.4–2.6；复振幅三元相互作用及随机强迫，非“非局部”概念首次 |
| [Dai–Friedlander，Dyadic models for ideal MHD](https://arxiv.org/html/2104.09440v2) | Eqs. (1.3)–(1.4)、§6；双分量已有先例，零分量向量场不同 |
| [Dai–Friedlander，Uniqueness and non-uniqueness results for dyadic MHD models](https://arxiv.org/html/2107.04073v2) | Eqs. (1.2a)–(1.2b)、Theorems 1.1–1.4；按作者稿题名和版本引用 |
| [Barbato–Morandin–Romito，Smooth solutions for the dyadic model](https://arxiv.org/abs/1007.3401) | Remark 3.4；有限跨度二次系统的局部构造背景 |
| [Montagnani，A stationary solution for the mixed turbulence shell model](https://arxiv.org/abs/2009.07952) | Theorem 5.2；不变分布意义的 stationary，不等于确定性平衡点 |
| [Ajzner–Alexakis，arXiv:2307.15505](https://arxiv.org/abs/2307.15505) | §III 的惯性区定常构造先忽略黏性、外力；与黏性解族极限区分 |
| [Tuteri–Chibbaro–Alexakis，arXiv:2601.04788v1](https://arxiv.org/html/2601.04788v1) | §2.2、Prop. 2.3.1；无黏多分支定常递推，非正黏性唯一性反例 |


## 进一步方法与模型来源

下面的定位承接原文献索引。只有主题文档中明确读取的方程、定理及证明段才作为覆盖判断依据；“来源属于某主题”不代表该文所有结论都已核查。

| 原始来源与读取定位 | 对应主题 |
|---|---|
| [Bianchi–Morandin §3.1](https://arxiv.org/html/1607.01334) | 尺度退化交换族的补充定位 |
| [树模型 §6.3](https://arxiv.org/html/1207.2846) | 尺度退化交换族的补充定位 |
| [Cheskidov–Dai–Friedlander，Dyadic models for fluid equations: a survey](https://arxiv.org/html/2209.10203) | 尺度退化交换族的补充定位 |
| [Cheskidov–Zaya，Regularizing Effect of the Forward Energy Cascade in the Inviscid Dyadic Model](https://arxiv.org/html/1310.7612v4) | 尺度退化交换族的补充定位 |
| [Barbato–Morandin，Positive and non-positive solutions for an inviscid dyadic model](https://arxiv.org/html/1201.2693v1) | 尺度退化交换族的补充定位 |
| [Barbato–Flandoli–Morandin，Energy dissipation and self-similar solutions for an unforced inviscid dyadic model](https://arxiv.org/html/0811.1689v1) | 尺度退化交换族的补充定位 |
| [Barbato–Morandin，Stochastic inviscid shell models: well-posedness and anomalous dissipation](https://ar5iv.labs.arxiv.org/html/1211.2910) | 尺度退化交换族的补充定位 |
| [Vincenzi–Gibbon，How turbulent convective dynamos work / shell-model dissipation bounds](https://arxiv.org/html/2009.01583v1) | 尺度退化交换族的补充定位 |
| [Cheskidov，Dissipation anomaly and anomalous dissipation](https://ar5iv.labs.arxiv.org/html/2311.04182) | 尺度退化交换族的补充定位 |
| [L'vov–Podivilov–Pomyalov–Procaccia–Vandembroucq，Improved shell model of turbulence](https://arxiv.org/abs/chao-dyn/9803025) | 尺度退化交换族的补充定位 |
| [Ditlevsen，Turbulence and Shell Models（Cambridge）](https://www.cambridge.org/us/universitypress/subjects/physics/nonlinear-science-and-fluid-dynamics/turbulence-and-shell-models) | 尺度退化交换族的补充定位 |
| [Gürcan，A nested polyhedra model of turbulence](https://arxiv.org/html/1603.01503) | 尺度退化交换族的补充定位 |
| [Chen–de Wit–Fruchart–Toschi–Vitelli，Odd viscosity suppresses intermittency](https://arxiv.org/pdf/2404.18894) | 尺度退化交换族的补充定位 |
| [Kiselev–Zlatoš，On Discrete Models of the Euler Equation](https://arxiv.org/html/math/0507129v1) | 尺度退化交换族的补充定位 |
| [Cheskidov–Foias，On global attractors of the 3D Navier–Stokes equations](https://arxiv.org/html/math/0608475v2) | 尺度退化交换族的补充定位 |
| [Barbato–Morandin–Romito，Global regularity for a logarithmically supercritical hyperdissipative dyadic equation](https://arxiv.org/html/1403.2852v2) | 尺度退化交换族的补充定位 |
| [Zhang，Sharp Regularity Thresholds for Viscous Katz–Pavlović Dyadic Models](https://arxiv.org/html/2609.21312v1) | 尺度退化交换族的补充定位 |
| [Heywood，A Curious Phenomenon in a Model Problem, Suggestive of the Hydrodynamic Inertial Range and Smallest Scale of Motion](https://link.springer.com/article/10.1007/s00021-003-0092-4) | 尺度退化交换族的补充定位 |
| [Mailybaev，RG approach to the inviscid limit for shell models of turbulence](https://doi.org/10.1088/1361-6544/adf407) | 尺度退化交换族的补充定位 |
| [作者稿 §§2–5](https://arxiv.org/html/2408.04659) | 尺度退化交换族的补充定位 |
| [Kankaria–Pal–Titi–Ray，Dynamical slowdown, bottlenecks, and multiscaling in Voigt-regularised turbulence](https://arxiv.org/html/2608.27355) | 尺度退化交换族的补充定位 |
| [L’vov–Pierotti–Pomyalov–Procaccia，Anomalous Scaling from Controlled Closure in a Shell Model of Turbulence](https://arxiv.org/html/chao-dyn/9810021) | 冻结排列与稳态路线的收敛审计 |
| [Barbato–Morandin，Stochastic inviscid shell models: well-posedness and anomalous dissipation](https://arxiv.org/html/1211.2910) | 冻结排列与稳态路线的收敛审计 |
| [De Rosa–Drivas–Inversi–Isett，Intermittency and Dissipation Regularity in Turbulence](https://arxiv.org/abs/2502.10032v2) | 中间区间矩条件 (MB) 的定向排重定位 |
| [De Rosa–Isett，Intermittency and lower dimensional dissipation in incompressible fluids](https://arxiv.org/abs/2212.08176v2) | 中间区间矩条件 (MB) 的定向排重定位 |
| [Cheskidov–Zaya，Regularizing Effect of the Forward Energy Cascade in the Inviscid Dyadic Model](https://arxiv.org/html/1310.7612) | 中间区间矩条件 (MB) 的定向排重定位 |
| [Tao，arXiv:1402.0290v3](https://arxiv.org/html/1402.0290v3) | 时间闸门的来源补查 |
| [Gürcan，arXiv:2004.01134v2](https://arxiv.org/html/2004.01134v2) | 时间闸门的来源补查 |
| [Zimin–Hussain，Wavelet based model for small-scale turbulence，Physics of Fluids 7 (1995), 2925–2927](https://doi.org/10.1063/1.868669) | N1 原始来源缺口 |
| [出版方 PDF 入口](https://pubs.aip.org/aip/pof/article-pdf/7/12/2925/19197588/2925_1_online.pdf) | N1 原始来源缺口 |
| [3D Vector Wavelet-Based Subgrid Scale Model for LES of Nonequilibrium Turbulence，ADA299639](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA299639.xhtml) | N1 原始来源缺口 |
| [Internet Archive](https://archive.org/details/DTIC_ADA299639) | N1 原始来源缺口 |
| [Cheskidov–Friedlander–Pavlović，math/0610815，§3、§4 Lemma 4.1、Theorems 4.2/4.4](https://arxiv.org/html/math/0610815) | 纯置换结果及自由衰减选题的来源边界 |
| [Barbato–Bianchi–Flandoli–Morandin 1207.2846v1](https://arxiv.org/html/1207.2846v1) | 完整图的多分量归约及纯置换覆盖边界 |
| [NUS IMS 报告摘要档案](https://imsarchives.nus.edu.sg/oldwww/Programs/wbfst/abstracts.html) | Heywood 单一来源定位的补充 |
| [Tao 的门电路 §5 及 Theorem 4.2](https://arxiv.org/html/1402.0290v3#S5) | 旋转链的指定原文比较 |
| [BMR 的向量壳包络 §3.3.1](https://arxiv.org/html/1403.2852v2#S3.SS3.SSS1) | 旋转链的指定原文比较 |
| [Jeong–Li 的混合模型 Theorem 2 及矩证明](https://arxiv.org/html/1406.2423v1#S2) | 旋转链的指定原文比较 |
| [Cheskidov–Friedlander 的长期平均耗散 Theorem 4.2](https://arxiv.org/html/0810.3718#S4) | 旋转链的指定原文比较 |
| [Cheskidov–Dai–Palasek，2609.26790v1，§1.2 / Theorem 1.10](https://arxiv.org/html/2609.26790v1#S1.SS2) | 旋转链的指定原文比较 |
| [Barbato–Flandoli–Morandin 式(1)、Theorem 1](https://arxiv.org/html/1007.1004v1) | 同边际冻结刷新对照的窄范围来源 |
| [Andersen–Muratore-Ginanneschi §II、§V–VII](https://arxiv.org/html/chao-dyn/9902007v3) | 同边际冻结刷新对照的窄范围来源 |
| [BFM 1007.1004v1，§6 的 Propositions 17–18、Corollary 19 及 §6.1](https://arxiv.org/html/1007.1004v1#S6) | 概率一长期失能的单一来源补读 |
| [Luo–Wang 2207.09631v1，§1.4 与 §6](https://arxiv.org/html/2207.09631v1#S1.SS4) | D1 的有限覆盖核查 |
| [Tao §5.4](https://arxiv.org/html/1402.0290v3#S5.SS4) | 自治相位门与 Tao 旋转平均的具体先例 |
| [§5.5 Theorem 5.3、式(5.22)、(5.34)–(5.35)](https://arxiv.org/html/1402.0290v3#S5.SS5) | 自治相位门与 Tao 旋转平均的具体先例 |
| [Theorem 4.2](https://arxiv.org/html/1402.0290v3#S4) | 自治相位门与 Tao 旋转平均的具体先例 |
| [Genadot–Thieullen，Multiscale Piecewise Deterministic Markov Process in Infinite Dimension: Central Limit Theorem and Langevin Approximation，1211.1894v1，§2.2 Theorem 4](https://arxiv.org/html/1211.1894v1#S2.SS2) | D1 的一般平均化与 dyadic 拼接核查 |
| [§3.2 Proposition 9、式(34)](https://arxiv.org/html/1211.1894v1#S3.SS2) | D1 的一般平均化与 dyadic 拼接核查 |
| [Liu–Röckner–Sun–Xie，Strong averaging principle for slow-fast stochastic partial differential equations with locally monotone coefficients，1907.03260v2，式(1.1)、§2 式(2.1)、Theorem 2.4、Remark 2.5](https://arxiv.org/html/1907.03260v2#S2) | D1 的一般平均化与 dyadic 拼接核查 |
| [Barbato–Flandoli–Morandin，0811.1689v1，式(1)、Theorem 6、Lemma 7](https://arxiv.org/html/0811.1689v1#S3) | D1 的一般平均化与 dyadic 拼接核查 |
| [Theorem 8](https://arxiv.org/html/0811.1689v1#S4) | D1 的一般平均化与 dyadic 拼接核查 |
| [BFM，Anomalous dissipation in a stochastic inviscid dyadic model，1007.1004v1，§4–5、Propositions 7–9](https://arxiv.org/html/1007.1004v1#S5) | 刷新率完整分类的方法先例 |
| [Barbato–Morandin，Stochastic inviscid shell models: well-posedness and anomalous dissipation，1211.2910v2，§1 式(1)](https://arxiv.org/html/1211.2910v2#S1) | 分块相关分类的定向来源 |
| [§6 Proposition 14、Theorem 6](https://arxiv.org/html/1211.2910v2#S6) | 分块相关分类的定向来源 |
| [Biroli–Larralde–Majumdar–Schehr，Extreme statistics and spacing distribution in a Brownian gas correlated by resetting，2211.00563v2，式(3)–(6)](https://arxiv.org/html/2211.00563v2) | 共同重置与独立失步的方法先例 |
| [de Mauro–Biroli–Majumdar–Schehr，Dynamically emergent correlations in Brownian particles subject to simultaneous non-Poissonian resetting protocols，2509.06658v2，PDF 第 1–5 页、§I–II.A、式(1)–(15)](https://arxiv.org/pdf/2509.06658) | 渐弱共同重置分类的方法边界 |
| [Das–Giuggioli，Discrete space-time resetting model: Application to first-passage and transmission statistics，2209.08330v1，§2 式(1)–(2)、§2.3 式(5)–(8)](https://arxiv.org/html/2209.08330v1#S2) | 有限巨块阻断的定向方法比较 |
| [Paul C. Bressloff，Diffusive search for a stochastically-gated target with resetting，2009.02739v1，§2 式(2.1)–(2.9)、§3 式(3.1)–(3.11)](https://arxiv.org/html/2009.02739v1) | 指数迭代几何分类的门控方法比较 |

本次仅整理已有来源记录，没有启动新的文献检索或重新验收已收束的原创性判断。检索未命中不作为排除先例的证据。


## 早期三元相互作用公式的读取定位

以下仅保存旧报告中实际读取的公式定位；旧报告的全局等价排除不沿用。

| 来源 | 旧读取位置与范围 |
|---|---|
| [Benzi 等，chao-dyn/9510010](https://arxiv.org/abs/chao-dyn/9510010) | 作者 TeX 的 eq:shells，原存档源文件第503–513行；指定壳模型三元作用公式。 |
| [arXiv:1608.00742](https://arxiv.org/abs/1608.00742) | 作者 TeX 的 eq:sabra_helical_general_up / um，原存档第217–225行；指定螺旋壳相互作用。 |
| [arXiv:1607.05518](https://arxiv.org/abs/1607.05518) | 作者 TeX 的 eq:sabra_model3，原存档第754–764行；指定模型的非线性项。 |

源文件行号依当时下载版本，仅用于追溯；它们不构成所有坐标变换或一般模型的排除定理。
