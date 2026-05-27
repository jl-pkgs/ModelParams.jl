# CLAUDE.md — ModelParams.jl

## Project structure

```
src/
  ModelParams.jl          # module entry, include order matters
  metadata.jl             # @metadata, @bounds, @units macros
  MultiLayer.jl           # MultiLayer{FT,N,S} SoA container + parameters/update!
  Params.jl               # get_params, parameters, update!, get_opt_info
  GOF.jl                  # goodness-of-fit metrics
  Optim/Optim.jl          # SCE-UA optimizer
  Kv_Profile.jl           # Kv struct types (KvExp, KvExpConst, KvExpPiecewise, KvLayers)
  Kv.jl                   # kv_at_depth, kv_layer_ksat dispatch
  Model_SoilDiffEqs.jl    # Campbell, VanGenuchten, HydraulicProfile, ThermalProfile
  SoilColumn.jl           # AbstractSoilModel, SoilColumn, filter_params, update_params!
  Retention_ParamTable.jl # get_soilpar(:Campbell/:VanGenuchten, soil_type) look-up table
  Retention_PTF.jl        # campbell_from_ptf + scalar PTF functions

test/
  test-Kv_profile.jl      # kv_at_depth / kv_layer_ksat for all Kv types
  test-SoilDiffEqs.jl     # retention types, HydraulicProfile, ThermalProfile
  test-SoilColumn.jl      # filter_params / update_params! on SoilColumn
  Model_SoilDiffEqs.jl    # integration test: SoilColumn round-trip
  Model_BEPS_base.jl      # ParamVeg definition (shared by BEPS tests)
  Model_BEPS_2.0.jl       # ParamBEPS2 definition + tests
  Model_BEPS_3.0.jl       # ParamBEPS3 <: AbstractSoilModel + tests
  test-type_stability.jl
  test-par_map.jl
  test-sceua.jl
```

## Key design patterns

### 1. `@bounds @units` annotation syntax

```julia
@bounds @units @with_kw mutable struct Foo{T} <: AbstractRetention{T}
    x::T = 1.0 | (0.0, 10.0) | "m"   # bound=(0,10), unit="m"
    y::T = 2.0 | nothing       | "-"  # bound=nothing → excluded from calibration
end
```

Macros generate zero-cost dispatch methods at compile time. `bounds(Foo, :x)` and `units(Foo, :x)` query them.

### 2. `MultiLayer{FT,N,S}` — SoA container

- Automatically extracts all `AbstractFloat` fields from scalar type `S` into `Vector{FT}` arrays.
- `profile[i]` → reconstructs a scalar `S` instance (AoS slice). Used to rebuild `.layers` cache.
- `get_params(x::MultiLayer)` returns one row per (field, layer) pair with `bound ≠ nothing`.

### 3. `parameters(model)` traversal

`get_params(x::S; path, with_unit)` recurses into struct fields:
- Sub-struct or `MultiLayer` (`has_definedbounds` / `isstructtype`) → `:predef` → recurse
- Scalar with `@bounds` annotation → `:macro` → emit one row
- `AbstractArray` / `Tuple` → `:macro` (no recursion)
- Explicit `| nothing` → `:skip`

Paths are `Vector{Symbol/Int}` like `[:hydraulic, :profile, :θ_sat, 3]`.

**Critical:** never define a 1-arg `get_params(::MyModel; kw...)`. It overrides the generic `get_params(x::S; ...)` and causes `parameters → get_params → parameters` infinite recursion. Use a different function name (e.g. `filter_params`) for model-specific filtering.

### 4. `filter_params` — calibration subset

Generic function (no type constraint) in `SoilColumn.jl`. Filters the `parameters()` DataFrame by:
- `mod`: `:hydraulic` | `:thermal` | `:veg` | `[:a,:b]` | `:all` → matches `path[1]`
- `inds`: keep only rows where `path[end] isa Integer && path[end] ∈ inds`
- `list_fix`: exclude rows by `name`
- `list_sameLayer`: dedup by `name`, keep first occurrence; `update_params!` broadcasts the value to all layers

KvExpPiecewise auto-fix: `hasfield(typeof(ps), :hydraulic) && ps.hydraulic.kv isa KvExpPiecewise` → `:z_exp` added to `fix_set`.

### 5. `update_params!` on `AbstractSoilModel`

Single shared implementation. Required fields on subtypes:
- `ps.hydraulic::HydraulicProfile` (has `.profile`, `.layers`, `.kv`, `.dz_cm`)
- `ps.thermal::ThermalProfile` (has `.profile`, `.layers`)

Pipeline order (must not reorder):
1. `update!(ps, paths, theta; params)` — write to SoA profile
2. `fill!(profile.field, value)` — broadcast `list_sameLayer`
3. `update_hydraulic!(profile)` — VG: `m = 1 − 1/n`; others: no-op
4. `_sync_ksat!(kv, profile, dz_cm)` — layer-integrated Ksat from kv profile
5. Rebuild AoS: `layers[i] = profile[i]` for hydraulic and thermal

### 6. Defining a new soil model

```julia
@bounds @with_kw mutable struct MyModel{FT,N,H,T} <: AbstractSoilModel{FT,N}
    hydraulic::H    # HydraulicProfile{FT,N,...}
    thermal::T      # ThermalProfile{FT,N,...}
    my_param::FT = 1.0 | (0.1, 10.0)
end
function MyModel{FT,N}(...) where {FT,N} ... end   # outer constructor
```

`filter_params` and `update_params!` are inherited for free.

## Testing

Run all tests:
```
julia --project -e "using Pkg; Pkg.test()"
```

Run a single file quickly:
```
julia --project -e "include(\"test/test-SoilColumn.jl\")"
```

Type-stability: `@inferred MyStruct{Float64,5}()` must return a fully concrete type (no free type parameters). Check with `@code_warntype`.

## Soil parameter initialisation

### `get_soilpar` — texture look-up table (`Retention_ParamTable.jl`)

Three calling forms (all equivalent):
```julia
get_soilpar(:Campbell, 7)              # Symbol (recommended)
get_soilpar(Val{:Campbell}, 7)         # Type (also works after recent change)
get_soilpar(7; retention="Campbell")   # string convenience wrapper
```
`soil_type` 1–12 follows USDA classification (1=Clay … 12=Sand). The `verbose=true` kwarg prints the texture name (requires `USDA` module in scope).

### `campbell_from_ptf` — pedotransfer functions (`Retention_PTF.jl`)

```julia
campbell_from_ptf(clay, silt, sand, bd, ph; ksat_method=:brakensiek)
# → Campbell{Float64}
```

Scalar PTF functions (all `@inline`, type-stable):

| Function | Output | Source |
|---|---|---|
| `θsat_toth(ph, bd, clay, silt)` | θ_sat [m³/m³] | Tóth et al. 2015 |
| `θres_rawls_brakensiek(sand, clay, θsat)` | θ_res [m³/m³] | Rawls & Brakensiek 1989 |
| `pore_size_index_brakensiek(sand, θsat, clay)` | λ [-] | Rawls & Brakensiek 1989 |
| `b_cosby(clay, sand)` | b [-] | Cosby et al. 1984 |
| `psi_sat_cosby(sand)` | ψ_sat [cm, negative] | Cosby et al. 1984 |
| `kv_brakensiek(θsat, clay, sand)` | Ksat [cm h⁻¹] | Brakensiek et al. 1984 |
| `kv_cosby(sand, clay)` | Ksat [cm h⁻¹] | Cosby et al. 1984 |

Unit conversion constant: `_MM_DAY_TO_CM_H = 1/240`.

## Parameter row counts (N=5 reference)

| Profile type | Bounded fields | Rows |
|---|---|---|
| CampbellLayers | θ_sat, ψ_sat, b | 15 |
| VanGenuchtenLayers | θ_sat, θ_res, α, n | 20 |
| KvLayers | kv | 5 |
| KvExpLayers | kv, f | 10 |
| ThermalMainLayers | κ, cv | 10 |
| ThermalBaseLayers | κ_dry, ρ_soil, V_SOM | 15 |
