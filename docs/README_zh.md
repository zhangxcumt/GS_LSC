# GS-LSC 代码包 —— 中文快速上手

本目录是论文

> **Regional gravity field modelling by least-squares collocation with a block
> Gauss-Seidel iterative solver**（Zhang, Chang & Chang，投稿 *Computers &
> Geosciences*）

的配套代码包。英文说明见根目录 `README.md`，算法公式见 `docs/ALGORITHMS.md`。

## 1. 环境

* MATLAB R2018b 或更新版本，**不需要任何工具箱**。
* 复现论文全部实验需要约 32 GB 内存；只跑快速测试 8 GB 足够。
* 块 Gauss-Seidel 只需要 128 MB。

## 2. 一分钟快速测试

在 MATLAB 中：

```matlab
cd examples
run_quick_test
```

或在仓库根目录：

```matlab
run('examples/run_quick_test.m')
```

两个算例依次运行：

1. `example_known_solution.m`：由已知系数构造 `y = (C+D) alpha_true`，检验三种
   解法都能还原 `alpha_true`（秒级，不需要数据）；
2. `example_colorado_subset.m`：在科罗拉多真实数据的 5 895 个点子上跑三种解法，
   比较 222 个 GNSS/水准检核点上的高程异常（约 15 s）。

两项都通过时打印 `QUICK TEST PASSED`；任一项超差即报错退出。

## 3. 复现论文的表格与图

`scripts/` 下的脚本各自独立（开头自带 `addpath`），可从任意目录运行。结果写入
`results/`，图件写入 `figures/`。

| 论文条目 | 脚本 | 耗时（N = 50 902） |
|---|---|---|
| 表 3 标准 LSC | `run_standard_lsc` | 约 2 min，峰值 20 GB |
| 表 3 序贯 LSC | `run_sequential_lsc` | 约 25 min |
| 表 3 分块 Gauss-Seidel | `run_block_gs_lsc` | 约 35 min，128 MB |
| 表 3 汇总 | `make_table3` | 秒级 |
| 表 4 单次扫描耗时 | `run_table4_sweep_cost` | 约 1 min |
| 图 1 存储量 | `make_fig1_storage` | 秒级 |
| 图 4 GSVS17 误差 | `make_fig4_gsvs17` | 秒级 |
| 图 5 收敛过程 | `make_fig5_convergence` | 约 1.5 h |

一次跑完全部：`cd scripts; run_all`。

先用小样本验证整条链路（推荐）：

```powershell
$env:GSLSC_NMAX = 8000      # PowerShell
```
```
set GSLSC_NMAX=8000         # cmd
```
```matlab
cd scripts
run_all
```

此时每个脚本只用 8 000 个观测点，全部跑完约 5 min。

## 4. 论文中的数值设置

均在 `functions/gslsc_config.m` 中，可直接修改：

| 量 | 取值 |
|---|---|
| 协方差参数 (A, B, M, R_B) | (21, 24, 719, 6 369 500 m) |
| 观测噪声方差 D | 16 mGal² |
| 子区域 | 34.972–39.171 °N，108.203–104.003 °W |
| 子区域内观测数 | 50 902（论文中写作 5 × 10⁴） |
| GNSS/水准检核点 | 222 |
| 分块大小 n_b | 2000、4000 |
| 迭代停止阈值（预测量） | 0.01 cm |

## 5. 注意事项

* **N 的差别**：论文表 2 的经纬度范围（按该范围直接筛选）包含 **50 902** 个观测，
  脚本默认使用全部 50 902 个点；论文正文写作 5 × 10⁴。若要与论文数字严格对齐，
  在 `gslsc_config.m` 中设 `cfg.nMax = 50000`（固定随机种子抽样），或直接修改
  表 2 的观测数。
* **序贯 LSC 的内存**：`lsc_sequential` 会在递推状态 `P` 超过 `opt.maxMemGB`
  （默认 26 GB）时直接报错并提示改用分块 Gauss-Seidel —— 这正是论文要说明的
  结论，不是程序缺陷。
* **迭代法的残差**：`lsc_block_gauss_seidel` 只在迭代结束时给出精确残差
  `info.relResid`；扫描内部的逐块残差恒为零（每个块都被精确求解），不要用它
  判断收敛。停止判据应基于预测量的变化（`opt.thetaFun` + `opt.tolPred`）。
