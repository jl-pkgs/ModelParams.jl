export kv_at_depth, effective_ksat, kv_layer_ksat

# ──────────────────────────────────────────────
# Unified dispatch: kv_at_depth(profile, layer, z_cm)
#
# Arguments:
#   layer  : 1-based layer index
#   z_cm   : center depth of layer [cm, positive downward]
# ──────────────────────────────────────────────
@inline kv_at_depth(kv::KvExp, ::Int, z_cm::Real) =
    kv.kv * exp(-kv.f * z_cm)

@inline kv_at_depth(kv::KvExpConst, ::Int, z_cm::Real) =
    kv.kv * exp(-kv.f * min(z_cm, kv.z_exp))

@inline function kv_at_depth(kv::KvExpPiecewise, ::Int, z_cm::Real)
    n = length(kv.kv)
    z_start = 0.0
    @inbounds for i in 1:n
        z_end = kv.z_exp[i]
        z_cm <= z_end && return kv.kv[i] * exp(-kv.f[i] * (z_cm - z_start))
        z_start = z_end
    end
    # below last segment boundary: constant tail
    @inbounds begin
        dz_last = kv.z_exp[n] - (n > 1 ? kv.z_exp[n-1] : 0.0)
        kv.kv[n] * exp(-kv.f[n] * dz_last)
    end
end

@inline kv_at_depth(kv::KvLayers, i::Int, ::Real) = kv.kv[i]

# @inline function kv_at_depth(kv::KvExpLayers, i::Int, z_cm::Real,
#                               nlayers_kv::Int, z_layered_cm::Real)
#     i <= nlayers_kv && return kv.kv[i]
#     # exponential continuation from the last layered-zone Ksat value
#     kv.kv[nlayers_kv] * exp(-kv.f[nlayers_kv] * (z_cm - z_layered_cm))
# end

# ──────────────────────────────────────────────
# Layer-integrated Ksat: kv_layer_ksat(profile, i, z1_cm, z2_cm)
#
# Returns the thickness-weighted average Ksat for layer i spanning [z1_cm, z2_cm] (cm, positive downward).
# Used by _sync_ksat! to pre-compute param[i].Ksat once, so the solver hot path
# needs no runtime ratio correction.
# ──────────────────────────────────────────────

# ∫[z1,z2] kv·exp(−f·z) dz / (z2−z1), falls back to centre-point when layer is thin.
@inline function _kv_integral(kv::T, f::T, z1::T, z2::T) where {T<:Real}
    dz = z2 - z1
    f * dz < T(1e-8) && return kv * exp(-f * (z1 + z2) / 2)
    kv / (f * dz) * (exp(-f * z1) - exp(-f * z2))
end

@inline function kv_layer_ksat(kv::KvExp, ::Int, z1_cm::T, z2_cm::T) where {T<:Real}
    _kv_integral(T(kv.kv), T(kv.f), z1_cm, z2_cm)
end

@inline function kv_layer_ksat(kv::KvExpConst, ::Int, z1_cm::T, z2_cm::T) where {T<:Real}
    kv0 = T(kv.kv)
    f = T(kv.f)
    z_exp = T(kv.z_exp)
    ksat_const = kv0 * exp(-f * z_exp)
    z1_cm >= z_exp && return ksat_const
    z2_cm <= z_exp && return _kv_integral(kv0, f, z1_cm, z2_cm)
    # layer straddles z_exp — weighted average of exponential and constant parts
    d_exp = z_exp - z1_cm
    d_const = z2_cm - z_exp
    dtot = z2_cm - z1_cm
    (_kv_integral(kv0, f, z1_cm, z_exp) * d_exp + ksat_const * d_const) / dtot
end

@inline function kv_layer_ksat(kv::KvExpPiecewise, ::Int, z1_cm::T, z2_cm::T) where {T<:Real}
    n = length(kv.kv)
    total = T(0)
    z_seg_start = T(0)
    done = false
    @inbounds for i in 1:n
        z_seg_end = T(kv.z_exp[i])
        if z1_cm < z_seg_end
            seg_z1 = max(z1_cm, z_seg_start)
            seg_z2 = min(z2_cm, z_seg_end)
            local_z1 = seg_z1 - z_seg_start
            local_z2 = seg_z2 - z_seg_start
            total += _kv_integral(T(kv.kv[i]), T(kv.f[i]), local_z1, local_z2) * (seg_z2 - seg_z1)
        end
        z_seg_start = z_seg_end
        if z2_cm <= z_seg_end
            done = true
            break
        end
    end
    # below all segment boundaries: constant tail
    if !done
        @inbounds begin
            dz_last = T(kv.z_exp[n]) - (n > 1 ? T(kv.z_exp[n-1]) : T(0))
            ksat_const = T(kv.kv[n]) * exp(-T(kv.f[n]) * dz_last)
        end
        total += ksat_const * (z2_cm - max(z1_cm, z_seg_start))
    end
    total / (z2_cm - z1_cm)
end

@inline kv_layer_ksat(kv::KvLayers, i::Int, ::T, ::T) where {T<:Real} = T(kv.kv[i])

# @inline function kv_layer_ksat(kv::KvExpLayers, i::Int, z1_cm::T, z2_cm::T,
#                                 nlayers_kv::Int, z_layered_cm::T) where {T<:Real}
#     i <= nlayers_kv && return T(kv.kv[i])
#     kv_base = T(kv.kv[nlayers_kv])
#     f_base  = T(kv.f[nlayers_kv])
#     _kv_integral(kv_base, f_base, z1_cm - z_layered_cm, z2_cm - z_layered_cm)
# end

"""
    effective_ksat(ps, i, z_cm) → Ksat [cm h⁻¹]

Return point-estimate Ksat at centre depth z_cm for layer i.
Primarily for inspection — the solver uses `ps.param_hydraulic[i].Ksat` (pre-integrated).
The default `kv_profile` is `KvLayers`, initialized from hydraulic Ksat.
"""
function effective_ksat(ps::SoilModel, i::Int, z_cm::Real)
    kv_at_depth(ps.hydraulic.kv, i, z_cm)
end
