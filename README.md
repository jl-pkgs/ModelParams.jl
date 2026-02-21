<h1>ModelParams.jl</h1>

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/ModelParams.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/ModelParams.jl/dev)
[![CI](https://github.com/jl-pkgs/ModelParams.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/ModelParams.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/ModelParams.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jl-pkgs/ModelParams.jl/tree/master)

> Dongdong Kong

## 1 ModelParams

> 参数自动堆叠

## 2 SCEUA Optimization algorithm

> 参数自动优化

- `n_pop` 是全班人数

- `n_complexes` ：分成几组

- `n_complex_size`： 每组人数

- `n_simplex_size` ：每次小组讨论参与人数

  > 每次讨论一定有“组内第一名”，其余名额从组内按“越优秀越容易被抽中”的概率补齐。

- `n_evolu_steps` ：每组每轮要讨论几次


### 2.1 高维参数（50-100维）取值建议

> 记 `d = 参数维度 (50 <= d <= 100)`

| 场景 | size_simplex | size_complex | n_evolu | n_complex | n_popu (= size_complex * n_complex) | maxn |
|---|---:|---:|---:|---:|---:|---:|
| 保守省预算 | `d + 1` | `d + 3 ~ 1.2d` | `size_simplex` | `3 ~ 4` | `4d ~ 5d` | `8n_popu ~ 12n_popu` |
| 平衡（推荐） | `d + 1` | `1.2d ~ 1.5d` | `0.5size_complex ~ size_complex` | `4 ~ 6` | `5d ~ 6d` | `12n_popu ~ 20n_popu` |
| 强全局搜索 | `d + 1` | `1.5d ~ 2d` | `size_complex ~ 1.5size_complex` | `6 ~ 8` | `6d ~ 8d` | `15n_popu ~ 25n_popu` |

- `d=50`（平衡）：`size_complex=60~75`, `n_complex=4~6`, `n_popu=250~300`, `maxn=3000~5000`
- `d=100`（平衡）：`size_complex=120~150`, `n_complex=4~6`, `n_popu=500~650`, `maxn=6000~10000`

> 注：若单次函数评估很昂贵（>1s），可将 `n_popu` 与 `maxn` 同时乘以 `0.7` 作为起点。

### 2.2 高维优化流程（50-100维）

> 50-100 维时，纯 SCE-UA 往往成本高、收敛慢。更高效的策略是：先降难度，再分阶段优化。

- 第 1 阶段（降难度）：做敏感性筛选、收紧边界、参数归一化（必要时对正值参数用 `log` 变换）。
- 第 2 阶段（全局粗搜）：用较保守预算的 SCE-UA 获取一批候选解，重点保证覆盖而非一次到最优。
- 第 3 阶段（局部精修）：以候选解为初值，在邻域内使用 `Optim.jl` 的 `LBFGS`（可导）、`CMA-ES`（不可导且参数相关性强）或`NelderMead`（不可导）做局部优化。

> 说明：本项目优化过程支持 `seed` 控制，结果可复现；因此不单列“稳健性复验”阶段。

#### 2.2.1 CMA-ES vs SCE-UA

- `CMA-ES` 会在线学习参数协方差与步长，能自动捕捉变量相关性和局部曲率。
- 在狭长谷地、病态缩放、强耦合参数场景中，`CMA-ES` 通常比 `SCE-UA` 的单纯形演化更快逼近局部极值。
- `SCE-UA` 更偏全局探索；进入局部盆地后，其更新方向与步长自适应能力通常弱于 `CMA-ES`。

#### 2.2.2 CMA-ES vs NelderMead（局部阶段）

- `NelderMead` 适合低维、低成本、快速试探。
- `CMA-ES` 在中高维、强相关、尺度差异大时通常更稳、更快收敛。

> `NelderMead` 不显式学习变量相关结构；`CMA-ES` 通过协方差自适应刻画局部几何。
> 对于 `d=50~100` 的局部精修，优先 `CMA-ES`；`NelderMead` 一般仅作基线或小规模对照。
