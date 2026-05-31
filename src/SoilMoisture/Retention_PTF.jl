# https://github.com/Deltares/hydromt_wflow/blob/main/hydromt_wflow/workflows/ptf.py
export θsat_toth, θres_rawls_brakensiek, pore_size_index_brakensiek
export kv_brakensiek, kv_cosby, b_cosby, psi_sat_cosby
export campbell_from_ptf

const _MM_DAY_TO_CM_H = 1.0 / 240.0  # 1 mm day-1 = 1/240 cm h-1

"""
    θsat_toth(ph, bd, clay, silt) → θ_sat [m³/m³]

Saturated water content. Tóth et al. (2015), Eur. J. Soil Sci. 66, 226–238.
Inputs: pH [-], bulk density [g cm⁻³], clay [%], silt [%].
"""
@inline function θsat_toth(ph::T, bd::T, clay::T, silt::T) where {T<:Real}
    (0.5653
     - 0.07918  * bd^2
     + 0.001671 * ph^2
     + 5.438e-4 * clay
     + 0.001065 * silt
     + 0.06836
     - 1.382e-5 * clay * ph^2
     - 1.270e-5 * silt * clay
     - 4.784e-4 * bd^2 * ph^2
     - 2.836e-4 * silt * bd^2
     + 4.158e-4 * clay * bd^2
     - 0.01686  * bd^2
     - 3.541e-4 * silt
     - 3.152e-4 * ph^2)
end

"""
    θres_rawls_brakensiek(sand, clay, θsat) → θ_res [m³/m³]

Residual water content. Rawls & Brakensiek (1989).
Inputs: sand [%], clay [%], θ_sat [m³/m³].
"""
@inline function θres_rawls_brakensiek(sand::T, clay::T, θsat::T) where {T<:Real}
    (-0.0182482
     + 8.7269e-4 * sand
     + 5.13488e-3 * clay
     + 0.02939286 * θsat
     - 1.5395e-4 * clay^2
     - 1.0827e-3 * sand * θsat
     - 1.8233e-4 * clay^2 * θsat^2
     + 3.0703e-4 * clay^2 * θsat
     - 2.3584e-3 * θsat^2 * clay)
end

"""
    pore_size_index_brakensiek(sand, θsat, clay) → λ [-]

Brooks-Corey pore-size distribution index. Rawls & Brakensiek (1989).
Campbell b = 1/λ.
Inputs: sand [%], θ_sat [m³/m³], clay [%].
"""
@inline function pore_size_index_brakensiek(sand::T, θsat::T, clay::T) where {T<:Real}
    exp(-0.7842831
        + 0.0177544  * sand
        - 1.062498   * θsat
        - 5.304e-5   * sand^2
        - 2.73493e-3 * clay^2
        + 1.11134946 * θsat^2
        - 0.03088295 * sand * θsat
        + 2.6587e-4  * sand^2 * θsat^2
        - 6.10522e-3 * clay^2 * θsat^2
        - 2.35e-6    * sand^2 * clay
        + 7.98746e-3 * clay^2 * θsat
        - 6.74491e-3 * θsat^2 * clay)
end


"""
    kv_brakensiek(θsat, clay, sand) → K_sat [cm h⁻¹]

Saturated hydraulic conductivity. Brakensiek et al. (1984).
Inputs: θ_sat [m³/m³], clay [%], sand [%].
"""
@inline function kv_brakensiek(θsat::T, clay::T, sand::T) where {T<:Real}
    kv_mm_day = exp(19.52348   * θsat
                    - 8.96847
                    - 0.028212  * clay
                    + 1.8107e-4 * sand^2
                    - 9.4125e-3 * clay^2
                    - 8.395215  * θsat^2
                    + 0.077718  * sand * θsat
                    - 2.98e-3   * sand^2 * θsat^2
                    - 0.019492  * clay^2 * θsat^2
                    + 1.73e-5   * sand^2 * clay
                    + 0.02733   * clay^2 * θsat
                    + 1.434e-3  * sand^2 * θsat
                    - 3.5e-6    * clay^2 * sand) * 2.78e-6 * 1_000 * 3_600 * 24
    kv_mm_day * _MM_DAY_TO_CM_H
end


"""
    kv_cosby(sand, clay) → K_sat [cm h⁻¹]

Saturated hydraulic conductivity. Cosby et al. (1984), Water Resour. Res. 20(6), 682–690.
Inputs: sand [%], clay [%].
"""
@inline function kv_cosby(sand::T, clay::T) where {T<:Real}
    kv_mm_day = 60.96 * 10.0^(-0.6 + 0.0126 * sand - 0.0064 * clay) * 10.0
    kv_mm_day * _MM_DAY_TO_CM_H
end


"""
    b_cosby(clay, sand) → b [-]

Campbell shape parameter b. Cosby et al. (1984).
Inputs: clay [%], sand [%].
"""
@inline b_cosby(clay::T, sand::T) where {T<:Real} = 3.10 + 0.157 * clay - 0.003 * sand


"""
    psi_sat_cosby(sand) → ψ_sat [cm]

Air-entry matric potential at saturation. Cosby et al. (1984). Returns negative value.
Inputs: sand [%].
"""
@inline psi_sat_cosby(sand::T) where {T<:Real} = -10.1 * 10.0^(1.88 - 0.0131 * sand)


"""
    campbell_from_ptf(clay, silt, sand, bd, ph; ksat_method=:brakensiek) → Campbell{Float64}

Build Campbell hydraulic parameters from soil texture attributes using pedotransfer functions.
- θ_sat : Tóth et al. (2015)
- b     : 1 / pore_size_index, Rawls & Brakensiek (1989)
- ψ_sat : Cosby et al. (1984)
- K_sat  : Brakensiek et al. (1984) or Cosby et al. (1984)

Inputs: clay [%], silt [%], sand [%], bd [g cm⁻³], ph [-].
`ksat_method`: `:brakensiek` (default) or `:cosby`.
"""
function campbell_from_ptf(clay::Real, silt::Real, sand::Real, bd::Real, ph::Real;
                            ksat_method::Symbol=:brakensiek)
    clay, silt, sand, bd, ph = promote(Float64(clay), Float64(silt), Float64(sand),
                                       Float64(bd), Float64(ph))
    θ_sat = θsat_toth(ph, bd, clay, silt)
    λ     = pore_size_index_brakensiek(sand, θ_sat, clay)
    b     = 1.0 / λ
    ψ_sat = psi_sat_cosby(sand)
    K_sat  = if ksat_method === :brakensiek
        kv_brakensiek(θ_sat, clay, sand)
    elseif ksat_method === :cosby
        kv_cosby(sand, clay)
    else
        error("Unknown ksat_method: $ksat_method. Use :brakensiek or :cosby.")
    end
    Campbell{Float64}(; θ_sat, ψ_sat, K_sat, b)
end
