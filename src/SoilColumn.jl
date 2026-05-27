export update_params!, filter_params, AbstractSoilModel


# Subtypes get filter_params and update_params! for free.
abstract type AbstractSoilModel{FT,N} end

##
@with_kw mutable struct SoilColumn{FT<:AbstractFloat,N,
    H<:HydraulicProfile{FT,N},T<:ThermalProfile{FT,N}} <: AbstractSoilModel{FT,N}
    hydraulic::H
    thermal::T
end

function SoilColumn{FT,N}(
    hydraulic=HydraulicProfile{FT,N}(),
    thermal=ThermalProfile{FT,N}()) where {FT<:AbstractFloat,N}
    SoilColumn{FT,N,typeof(hydraulic),typeof(thermal)}(hydraulic, thermal)
end

# Sync Ksat from kv profile into the SoA hydraulic profile.
function _sync_ksat!(kv::Union{AbstractKv,AbstractKvLayers},
    profile::AbstractRetentionLayers{FT,N},
    dz_cm::AbstractVector) where {FT,N}
    z_cm = FT(0)
    @inbounds for i in 1:N
        dz_i = isempty(dz_cm) ? FT(0) : FT(dz_cm[i])
        profile.Ksat[i] = kv_layer_ksat(kv, i, z_cm, z_cm + dz_i)
        z_cm += dz_i
    end
end


# mod:            :hydraulic | :thermal | [:hydraulic, :thermal] | :all
# list_fix:       field names excluded from calibration, matched via params.name
# list_sameLayer: field names shared across all hydraulic layers (e.g. Campbell ψ_sat);
#                 returns one deduped row per name; update_params! broadcasts it back.
#
# Generic: works on any struct whose parameters() returns a standard DataFrame.
# KvExpPiecewise.z_exp is auto-excluded when ps.hydraulic.kv isa KvExpPiecewise.
function filter_params(ps, mod::Union{Symbol,AbstractVector{Symbol}};
    inds=nothing,
    list_sameLayer::Vector{Symbol}=Symbol[],
    list_fix::Vector{Symbol}=Symbol[])

    params = parameters(ps)
    n = nrow(params)
    mask = trues(n)

    # module filter
    if mod !== :all
        mods = mod isa Symbol ? (mod,) : Tuple(mod)
        for (i, p) in enumerate(params.path)
            mask[i] && p[1] ∉ mods && (mask[i] = false)
        end
    end

    # layer index filter
    if !isnothing(inds)
        inds_set = Set(inds)
        for (i, p) in enumerate(params.path)
            mask[i] && isa(p[end], Integer) && p[end] ∉ inds_set && (mask[i] = false)
        end
    end

    # list_fix: z_exp is a design parameter for KvExpPiecewise, always excluded
    fix_set = Set(list_fix)
    hasfield(typeof(ps), :hydraulic) && ps.hydraulic.kv isa KvExpPiecewise &&
        push!(fix_set, :z_exp)
    for (i, name) in enumerate(params.name)
        mask[i] && name ∈ fix_set && (mask[i] = false)
    end

    # list_sameLayer: keep only first occurrence per name
    same_set = Set(list_sameLayer)
    seen = Set{Symbol}()
    for (i, name) in enumerate(params.name)
        mask[i] || continue
        name in same_set || continue
        name in seen ? (mask[i] = false) : push!(seen, name)
    end
    params[mask, :]
end


# Shared update logic for AbstractSoilModel subtypes.
# Requires ps.hydraulic::HydraulicProfile and ps.thermal::ThermalProfile.
#
# list_sameLayer: broadcast the source-layer value to all hydraulic profile layers
# list_fix:       excluded from params by filter_params, so update! never touches them
function update_params!(ps::AbstractSoilModel{FT,N}, paths, theta;
    params=nothing,
    list_sameLayer::Vector{Symbol}=Symbol[],
    list_fix::Vector{Symbol}=Symbol[]) where {FT<:Real,N}

    isnothing(params) && (params = filter_params(ps, :hydraulic; list_sameLayer, list_fix))
    length(paths) == length(theta) ||
        throw(ArgumentError("paths and theta must have the same length, got $(length(paths)) and $(length(theta))."))

    update!(ps, paths, theta; params)

    # broadcast list_sameLayer from source layer to all layers before update_hydraulic!
    for (i, name) in enumerate(params.name)
        name in list_sameLayer || continue
        I = params.path[i][end]
        I isa Integer || continue
        h = getproperty(ps.hydraulic.profile, name)
        fill!(h, h[I])
    end

    # update derived fields (e.g. VanGenuchten m = 1 − 1/n); must precede _sync_ksat!
    update_hydraulic!(ps.hydraulic.profile)
    _sync_ksat!(ps.hydraulic.kv, ps.hydraulic.profile, ps.hydraulic.dz_cm)

    # rebuild AoS caches from updated SoA profiles
    @inbounds for i in 1:N
        ps.hydraulic.layers[i] = ps.hydraulic.profile[i]
        ps.thermal.layers[i] = ps.thermal.profile[i]
    end
    return nothing
end
