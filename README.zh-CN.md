# 壳模型中的异常耗散：两项研究成果

[English](README.md)

本仓库集中呈现两项研究：确定性双通道壳网络中，零密度系数改动可以抑制异常耗散；随机系数标量壳链中，单边完整时间统计相同仍可因边间相关组织不同而产生相反耗散结局。两项研究均明确方程、解类、初边值条件和实际黏性损失，提供完整解析证明。

## 第一项成果：零密度改动的耗散对照

研究受迫、零初值的正性双通道网络，黏性权重为 K_n²、K_n=2^n。两组配置保持边图、外力、初值、黏性和弱层自然密度不变：周期交换配置在一个足够长的固定窗口有黏性一致的正损失下界；只在自然密度为零的连续块中削弱交换，即可使每个固定窗口的损失趋于零。改动层的相对系数变化可以很大。

- [模型、主结果与阅读路线](docs/sparse-modification/model-zh.md)
- [主定理的完整证明](docs/sparse-modification/proofs-zh.md)
- [进一步结果：弱块、冻结排列、纯置换、自由衰减与稳定性](docs/sparse-modification/extensions-zh.md)
- [对数交换族的有效结论与未解区间](docs/sparse-modification/logarithmic-zh.md)
- [文献定位与已有定理的适用范围](docs/sparse-modification/literature-zh.md)

## 第二项成果：相关随机环境的精确耗散分类

研究无迫、根层单位脉冲初值的正性标量链。每条耦合边在弱值与强值之间切换；独立有限块内的共同重置与私有翻转决定相关性。固定每条边完整的平稳时间路径分布，仍能得到实际零黏损失为零或为正的两种环境。

例如，每块含 j 条边、单边相关率为 j²：完全独立时每个固定窗口零损失；任意固定正比例的共同重置都会给出异常耗散。一般有限块长 ℓ_j、正有限相关率 r_j 及重置比例 θ_j∈[0,1] 的完整分类由下列费用级数给出：

```math
\mathcal C_j=\begin{cases}
\ell_j/r_j,&\theta_j=0,\\
\min\{\ell_j,\theta_j^{-1}\}/r_j+(\mathcal N_j-1)/(\theta_jr_j),&0<\theta_j<1,\\
1/r_j,&\theta_j=1,
\end{cases}
\qquad \mathcal C=\sum_j\mathcal C_j.
```

其中 𝒩_j 为模型定义的块内几何迭代次数。费用发散给几乎处处所有有限窗口损失趋零；费用收敛给确定窗口的正期望损失及共同正概率路径下界。概率一结论取样本依赖窗口；先取零黏下／上极限，再令时间趋于无穷，两者均达到全部初能。

- [模型、环境与统一分类定理](docs/correlated-environments/model-zh.md)
- [统一分类的完整证明](docs/correlated-environments/proofs-zh.md)
- [分类推论、几何反例与同步恢复阈值](docs/correlated-environments/consequences-zh.md)
- [补充方法与统一定理之外的有效结论](docs/correlated-environments/supplementary-methods-zh.md)
- [文献比较与具体数学增量](docs/correlated-environments/literature-zh.md)

D1 已按本项目固定的标准达到“具有实质原创性、可证明异常耗散的壳模型”这一目标。具体受控对照、统一控制与已有定理的比较保存在文献定位中。该判断的范围与依据保持明确。

## 其他有效成果与开放问题

- [自治旋转链：完整耗散证明与严格三层回流证书](docs/other-results/rotation-chain-zh.md)
- [自治相位门：相同控制器能量下的快慢耗散对照](docs/other-results/phase-gates-zh.md)
- [早期模型及有效技术推导](docs/other-results/earlier-models-zh.md)
- [局部结果、障碍与停止点](docs/other-results/partial-results-zh.md)
- [仍然开放的研究问题](docs/open-problems-zh.md)

## 证明、检查与来源

数学结论由研究正文中的解析论证支持。程序检查有限代数、有限证书和文档结构；其输出注明实际验证范围。文献适用性以已读取的具体方程和定理为依据。

```bash
python -m pip install -r requirements.txt
python tests/check_repository.py
python tests/check_identities.py
python tests/check_cut_obstruction.py
python tests/check_rotation_chain.py
python tests/check_repository_negative.py
```

旋转链完整证书使用精确有理数运算，耗时长于代数检查；仅检查代数可加 --algebra-only。默认仓库检查不依赖本地论文或下载文献。

- [文献入口](references/README.md)
- [关键纠错与结论使用边界](docs/history/corrections-zh.md)
- [研究简史](docs/history/research-history-zh.md)
- [原章节到新结构的迁移索引](docs/history/migration-map-zh.md)
- [目录设计与维护规则](docs/history/organization-zh.md)

## 本地材料与许可

论文、投稿信、编辑通信、第三方全文和构建工具集中在忽略的 local 目录。公开仓库包含论文所用研究结果及证明，不分发论文文件和投稿材料。

作者原创研究文档采用 [CC BY 4.0](LICENSES/CC-BY-4.0.txt)，原创程序采用 [MIT](LICENSES/MIT.txt)，具体适用范围见 [许可说明](LICENSE)。第三方材料保留其原许可。
