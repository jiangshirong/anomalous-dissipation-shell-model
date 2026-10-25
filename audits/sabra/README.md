# 历史附录：一份 Sabra 实现的非线性能量审查

这项工作早于当前正性壳模型，二者不共用方程或能量归一化。这里保存的是对**指定表达式**的历史审查，不声称上游最新版本仍未修改。

来源是 `research-line/functional-stability-theory` 的两个脚本：`scripts/turbulence/compute_goy_shell_dfc.py` 和 `scripts/turbulence/compute_shell_dfc_waterline_ledger.py`。用户已提交 [Issue #1](https://github.com/research-line/functional-stability-theory/issues/1)；该历史检查独立于本仓库两项主成果。

当时被审查的表达式为

$$
N_n=i\left[k_nu_{n+1}^{\ast}u_{n+2}
-\tfrac12k_{n-1}u_{n+1}^{\ast}u_{n-1}
-\tfrac14k_{n-2}u_{n-1}u_{n-2}\right].
$$

目标能量为 $E=\sum_n|u_n|^2$。取三个相邻壳 $u_n=1,u_{n+1}=1,u_{n+2}=i$，其余为零，直接计算给出 $dE/dt=-7k_n/2\ne0$。这不依赖时间积分、黏性或外力。

该结论只针对上式及此能量定义，不替作者决定其意图，不直接修改 Sabra 约定，也不推断 FST 或三维 NS 理论的真伪。

## 可选复现

在仓库根目录运行：

```bash
python -m pip install -r requirements-audit.txt
python audits/sabra/check_nonlinearity.py
```

脚本是已保存表达式的独立转写，不下载或修改上游。它会在脚本旁的 `sabra_audit_results/` 生成 CSV、JSON 和报告；生成文件不入 Git。历史结果中 $k_0=1$、三个激活壳为 1、2、3 时总变化率为 −7。该程序只复现已保存表达式。

重复 Issue 草稿和重复报告不再维护；表达式、明确反例及独立复现程序在此保留。
