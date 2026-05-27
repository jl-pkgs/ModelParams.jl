<h1>ModelParams.jl</h1>

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/ModelParams.jl/dev)
[![CI](https://github.com/jl-pkgs/ModelParams.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/ModelParams.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/ModelParams.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jl-pkgs/ModelParams.jl/tree/master)

> Dongdong Kong

Parameter management and soil model infrastructure for physically-based models in Julia.

## Overview

ModelParams.jl provides three integrated layers:

1. **Parameter framework** — annotate model structs with calibration bounds and units; recursively collect, filter, and update parameters via a unified DataFrame API.
2. **Soil column model** — type-stable N-layer hydraulic and thermal profiles with SoA/AoS dual representation, multiple Ksat profile types, and a calibration-aware update pipeline.
3. **Optimization** — SCE-UA global optimizer and goodness-of-fit metrics (NSE, KGE, R², RMSE, …).

## Parameter Framework

### Annotating a struct

```julia
@bounds @units @with_kw mutable struct Campbell{T<:Real} <: AbstractRetention{T}
    θ_sat::T = 0.287 | (0.25, 0.50) | "m3 m-3"
    ψ_sat::T = -10.0 | (-100.0, -5.0) | "cm"
    Ksat::T  = 34.0  | nothing        | "cm h-1"   # nothing → excluded from calibration
    b::T     = 4.0   | (2.0, 15.0)   | "-"
end
```

The `| bound | unit` pipe syntax is handled at macro-expansion time — zero runtime cost.

### Collecting and updating parameters

```julia
model = ParamBEPS3{Float64,5}()

# Collect all calibratable parameters (bound ≠ nothing) as a DataFrame
df = parameters(model)      # columns: path, name, value, type, bound, unit

# Initial values and box constraints for an optimizer
x0, lb, ub, paths = get_opt_info(model)

# Write new values back into the model
update!(model, paths, x_new)
```

### filter_params — calibration subset

```julia
# Hydraulic parameters only, ψ_sat shared across all layers, b fixed
p = filter_params(model, :hydraulic;
    list_sameLayer = [:ψ_sat],   # one representative row; update_params! broadcasts
    list_fix       = [:b])        # excluded from calibration

# Then update + sync all caches in one call
update_params!(model, p.path, x_new; params=p, list_sameLayer=[:ψ_sat])
```

`filter_params` is generic — it works on any struct whose `parameters()` returns a standard DataFrame. `update_params!` is defined once on `AbstractSoilModel` and inherited by all subtypes.

## Soil Column Model

### Architecture

```
AbstractSoilModel{FT,N}
  SoilColumn{FT,N}
    hydraulic::HydraulicProfile
      profile::MultiLayer{FT,N,P}   ← SoA — parameter source of truth
      layers::Vector{P}              ← AoS cache — rebuilt after each update_params!
      kv::K                          ← KvExp | KvExpConst | KvExpPiecewise | KvLayers
      dz_cm::Vector{FT}
    thermal::ThermalProfile
      profile::MultiLayer{FT,N,P}
      layers::Vector{P}
```

`MultiLayer{FT,N,S}` is a generated Struct-of-Arrays container that automatically extracts all `AbstractFloat` fields from scalar type `S` into `Vector{FT}` arrays. Index with `profile[i]` to get back an `S` instance (AoS slice).

### Retention models

| Type | Parameters with bounds |
|---|---|
| `Campbell` | `θ_sat`, `ψ_sat`, `b` |
| `VanGenuchten` | `θ_sat`, `θ_res`, `α`, `n` (m derived) |

### Ksat profiles

| Type | Description |
|---|---|
| `KvLayers` | Per-layer constant |
| `KvExp` | Exponential decay: `kv·exp(−f·z)` |
| `KvExpConst` | Exponential above `z_exp`, constant below |
| `KvExpPiecewise` | Piecewise exponential segments |

`kv_layer_ksat` computes the thickness-weighted average Ksat for a layer, used by `_sync_ksat!` to pre-compute `profile.Ksat[i]` once before the solver hot path.

### update_params! pipeline

```
update!(model, paths, theta)          # write values to SoA profile
fill!(profile.field, value)           # broadcast list_sameLayer to all layers
update_hydraulic!(profile)            # VanGenuchten: m = 1 − 1/n
_sync_ksat!(kv, profile, dz_cm)       # Ksat ← layer-integrated kv profile
layers[i] = profile[i]               # rebuild AoS cache
```

### Extending with a new model

```julia
@bounds @with_kw mutable struct MyModel{FT,N,H,T} <: AbstractSoilModel{FT,N}
    hydraulic::H
    thermal::T
    extra_param::FT = 1.0 | (0.1, 10.0)
end
```

Inherits `filter_params` and `update_params!` automatically.

## Goodness-of-fit

```julia
GOF(obs, sim)
# → (NSE=…, R2=…, KGE=…, R=…, RMSE=…, MAE=…, bias=…, bias_perc=…, n_valid=…)
```

## SCE-UA Optimizer

```julia
result = sceua(cost_fn, x0, lb, ub; maxn=10000, kstop=5)
```

Shuffled Complex Evolution — global optimizer well-suited for hydrological model calibration.
