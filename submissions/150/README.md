# 素数采样悬赏 150 档证明

对满足 $150\leq\log p\leq651/4$ 的每个素数 $p$ 和题设允许的每个尺度 $X$，分别在 $3\mid p-1$ 与 $3\nmid p-1$ 两个分支验证 $\max_d(A-R_d-T_d)>0$。每格比较的三项使用同一截点 $d$。186 个相接闭区间覆盖整段，全段最小余量为 **0.000053428808**。右端 $651/4=162.75$ 与已有覆盖相接。[题面与验收条件](../../PROBLEM.md)

| $\log p$ 的闭区间 | 格数 | 两个分支中的最小余量 | 复验入口（压缩包内） |
| --- | ---: | ---: | --- |
| $[150,4801/32]$ | 1 | 0.013132132883 | `first-window/hyperbola-next/replay.py` |
| $[4801/32,152]$ | 64 | 0.000555406741 | `analytic-coverage/verify_cells.py` |
| $[152,651/4]$ | 121 | 0.000053428808 | 上游原版 `prime_sampling.cli.grade` 与 `native-continuation/submission.json` |

## 复验

从仓库根目录用 Python 3.10 或更高版本运行。默认复验只用标准库，在临时目录解包并校验全部文件的 SHA-256，再重放两条分支并检查闭区间衔接。Windows 可用 `--temp-parent D:\Temp` 指定已有的短目录。

```sh
python -B submissions/150/verify.py
```

预期输出末两行如下。完整的逐格记录见 [verification.txt](verification.txt)。

```text
Closed cells verified: 186 minimum margin in 1e-12: 53428808
Threshold 150: strict positivity on the entire closed interval [150,651/4].
```

[证书压缩包](certificate.zip) 的 SHA-256 为 `8387a28eb5e1e33348b8a5ebb1424d1dae8b846fa04d96b92ad7da04cfeca419`。内含复算所用的上游版本 `63a2d41e8124320784c244336ff36ac5f2285418`、证明、源码、谱表和标量预算。

终版中段账本为 `analytic-coverage/release/positive-0.json`，后段账本为 `native-continuation/submission.json`。部分子目录的 README 保留了生成这些数据时的阶段性结果。

## 证明与核对点

**主项 $A$。** 在 $(X/32,X]$ 上对 Chebyshev 函数 $\theta$ 分部积分，对数分母保留十四阶正项。每个参数窗都核对显式误差的适用范围，再将结果向下取整。实现位于 `main-density-new/verify.py`。[Broadbent 等的显式界](https://arxiv.org/abs/2002.11068)

**标量误差 $R_d$。** 按最小素因子 $q$ 写成 $m=qn$，在 $\log q=26$ 拼接前缀与尾项。前缀使用有限 Brun 树和 Selberg 权。对无平方因子的余因子 $n$，分给各素因子 $s$ 的权 $\log(s)/\log(n)$ 之和恰为一。这使尾项中的素因子计数可用非负系数多项式的卷积上界，含平方的整数单列计费。64 个中间格的审计逐一核对秩、起讫端点和截点，保留费用的向外取整。推导位于 `narrow-prefix/proof.md` 和 `analytic-coverage/log-population-proof.txt`。[Büthe 的显式估计](https://arxiv.org/abs/1511.02032)

**双线性误差 $T_d$。** 按较小素因子的指数剖分双曲条带，在同一模数分支和同一几何范围内给矩形定价。原生网格通过最短路拼接相邻微带。对矩形伸出真实条带的部分，用两段相位权收紧宏带费用。每个宏带的选择仅引用覆盖其完整实参数窗的证书。矩形分解和相位证明见 `hyperbola-next/geometry-proof.md`、`gram-cycle/corner-phase.md` 及 `gram-cycle/recursive-phase.md`。

例如，64 格中最小的余量出现在 $[4807/32,601/4]$ 的含 3 分支，截点为 $67/200$。整数账本给出

$$
1.949477571681 - 0.431942520000 - 1.516979644940 = 0.000555406741.
$$

## 从源码重建底层输入

要保留解包后的文件，可在仓库根目录运行 `python -m zipfile -e submissions/150/certificate.zip cert150`。以下命令均从解包根目录运行。生成矩形界和起筛前缀需要 C++17 编译器，这些重建比默认复验耗时更长。

```sh
python -B analytic-coverage/replay_tails.py
python -B narrow-prefix/release_audit.py --repo prime-sampling --release analytic-coverage/release --replay-prefixes
python -B analytic-coverage/verify_sources.py --step spectra
python -B analytic-coverage/verify_sources.py --step macros
```
