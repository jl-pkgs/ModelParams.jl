export update_params!

# Sync Ksat from kv profile into the SoA hydraulic profile (writes to profile.Ksat[i]).
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


# list_fix:       field names excluded from calibration, matched via params.name
# list_sameLayer: field names shared across all hydraulic layers (e.g. Campbell ψ_sat);
#                 get_params keeps only the first row, update_params! broadcasts the value back
# mod:            :hydraulic, :thermal, [:hydraulic, :thermal], or :all
function get_params(ps::SoilModel, mod=:hydraulic; inds=nothing,
    list_sameLayer::Vector{Symbol}=Symbol[],
    list_fix::Vector{Symbol}=Symbol[],
    with_unit::Bool=true)

    # bypass SoilModel dispatch to call the generic struct-traversal get_params
    params = invoke(get_params, Tuple{Any}, ps; with_unit) |> DataFrame

    mods = mod === :all ? nothing :
           mod isa Symbol ? (mod,) : Tuple(mod)

    paths = isnothing(mods) ? params.path :
            filter(x -> x[1] ∈ mods, params.path)

    # z_exp is a manually-specified design parameter for KvExpPiecewise, not calibrated
    if ps.hydraulic.kv isa KvExpPiecewise
        paths = filter(x -> !(x[1] === :hydraulic && length(x) >= 3 &&
                              x[2] === :kv && x[3] === :z_exp), paths)
    end

    if !isnothing(inds)
        paths = filter(x -> !isa(x[end], Integer) || x[end] in inds, paths)
    end

    # filter params directly — avoids a recursive parameters(ps) call
    row_inds = indexin(paths, params.path)
    params2 = params[row_inds, :]

    if !isempty(list_fix)
        params2 = params2[map(n -> n ∉ list_fix, params2.name), :]
    end

    if !isempty(list_sameLayer)
        seen = Set{Symbol}()
        keep = trues(length(params2.name))
        for (i, name) in enumerate(params2.name)
            name in list_sameLayer || continue
            name in seen ? (keep[i] = false) : push!(seen, name)
        end
        params2 = params2[keep, :]
    end

    params2
end


function update_params!(ps::SoilModel{FT,N}, paths, theta;
    params=nothing,
    list_sameLayer::Vector{Symbol}=Symbol[],
    list_fix::Vector{Symbol}=Symbol[]) where {FT<:Real,N}

    isnothing(params) && (params = get_params(ps; list_sameLayer, list_fix))
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

    # sync Ksat from kv profile into SoA hydraulic profile
    _sync_ksat!(ps.hydraulic.kv, ps.hydraulic.profile, ps.hydraulic.dz_cm)

    # rebuild AoS caches from updated SoA profiles
    @inbounds for i in 1:N
        ps.hydraulic.layers[i] = ps.hydraulic.profile[i]
    end
    @inbounds for i in 1:N
        ps.thermal.layers[i] = ps.thermal.profile[i]
    end

    return nothing
end
