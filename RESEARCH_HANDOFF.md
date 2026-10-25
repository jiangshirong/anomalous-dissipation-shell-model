# 当前研究交接

更新日期：2026-10-25。本文件是当前结果索引；完整证明在主题正文，具体旧错误在纠错表。不要以旧章节或已撤回评价覆盖当前有效结论。

## 两项主成果

| 成果 | 当前有效结论 | 最短阅读路线 |
|---|---|---|
| 确定性双通道网络 | 共同零初值、外力、边图、黏性及弱层密度下，周期配置正耗散；零密度支撑的弱块修改可使所有固定窗损失趋零。 | [模型](docs/sparse-modification/model-zh.md) → [证明](docs/sparse-modification/proofs-zh.md) → [扩展](docs/sparse-modification/extensions-zh.md) |
| D1 相关随机系数 | 指定平稳独立有限块的任意正有限速率、有限几何与重置比例均由费用级数完整分类；单边完整时间路径律不变仍可给相反耗散。 | [模型](docs/correlated-environments/model-zh.md) → [统一证明](docs/correlated-environments/proofs-zh.md) → [推论](docs/correlated-environments/consequences-zh.md) |

第一成果的 GP 是稳态充要判据和实际动态正耗散充分条件，不能扩大成全部动态的充要分类。纯置换、自由衰减与固定黏性稳定性保留为进一步结果。[对数族](docs/sparse-modification/logarithmic-zh.md)的总损失在 1/2<α≤3/4 仍未分类；特别是 α=2/3 不能由已证通道定位代替。

D1 按已固定标准达到实质原创且可证明异常耗散的目标，[文献比较](docs/correlated-environments/literature-zh.md)保存具体依据。其正侧固定窗为正期望与同一正概率事件的路径下界；概率一结论允许样本依赖窗口，按先零黏 liminf/limsup、后长时的顺序恢复全部初能。[补充材料](docs/correlated-environments/supplementary-methods-zh.md)保留严格冻结、无限同步尾块、固定窗零损失事件、受控等价排除及强平均误差等独有结果。

## 其他有效研究

[旋转链 R](docs/other-results/rotation-chain-zh.md)有标准黏性下完整正耗散证明及有限三层回流证书；无限链真实低层回流未因有限证书而解决。[自治相位门 A](docs/other-results/phase-gates-zh.md)保留明确角速度与能量权重下的快慢对照，不能推广到未经证明的原参数。

[局部结果](docs/other-results/partial-results-zh.md)与[开放问题](docs/open-problems-zh.md)列出剩余缺口和停止原因。指定 D1 环境族的一般速率、任意有限几何及 θ 两端点已解决，不重新列为待办。后续研究需要新的概率结构或动力学控制。

## 文件与工作规则

- [研究简史](docs/history/research-history-zh.md)保存实质里程碑；[关键纠错](docs/history/corrections-zh.md)防止旧结论复用；[迁移表](docs/history/migration-map-zh.md)定位旧章节去向。
- 论文及投稿状态只在本地 local/papers 中维护。研究对话更新研究正文；投稿工作使用对应本地材料，不自动改已提交版本。
- 当前公开仓库保留真实历史，不包含论文、投稿材料或私有旧历史。默认不自动提交推送；以用户当轮明确授权为准。
- 本次整理起点为 a9c9d25a9e443428fcd738576b3cf28767be4e1a。整理前工作区含尚未同步的有效研究，必须一并保全；迁移不依赖旧远端快照代替当前工作区。
