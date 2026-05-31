function get_soilpar(::Type{Val{:VanGenuchten}}, soil_type::Int=1; verbose=false)
  # Bonan 2019, Table 8.3
  soilparam = [
    # θ_sat, θ_res, α (cm⁻¹), n, K_sat (cm h⁻¹), ID_USDA, ID_Bonan, name
    0.38 0.068 0.008 1.09 0.2;   #  1, 11, Clay
    0.36 0.070 0.005 1.09 0.02;  #  2, 10, Silty clay
    0.38 0.100 0.027 1.23 0.12;  #  3,  9, Sandy clay
    0.41 0.095 0.019 1.31 0.26;  #  4,  8, Clay loam
    0.43 0.089 0.010 1.23 0.07;  #  5,  7, Silty clay loam
    0.39 0.100 0.059 1.48 1.31;  #  6,  6, Sandy clay loam
    0.43 0.078 0.036 1.56 1.04;  #  7,  5, Loam
    0.45 0.067 0.020 1.41 0.45;  #  8,  4, Silty loam
    0.41 0.065 0.075 1.89 4.42;  #  9,  3, Sandy loam
    0.45 0.067 0.020 1.41 0.45;  #  10, NA, Silty, nodata, used Silty loam (4)
    0.41 0.057 0.124 2.28 14.59; #  11, 2, Loamy sand
    0.43 0.045 0.145 2.68 29.7   #  12, 1, Sand
  ]
  verbose && println("[Soil Texture]: ", USDA.SoilTexture(soil_type))
  θ_sat, θ_res, α, n, K_sat = soilparam[soil_type, :]
  VanGenuchten(; θ_sat, θ_res, α, n, K_sat)
end

# function get_soilpar(::Type{VanGenuchten}, theta::AbstractVector)
#   θ_sat, θ_res, K_sat, α, n = theta[1:5]
#   VanGenuchten(; θ_sat, θ_res, α, n, K_sat)
# end

function get_soilpar(::Type{Val{:Campbell}}, soil_type::Int=1; verbose=false)
  # Bonan 2019, Table 8.3
  soilparam = [
    # θ_sat, ψ_sat (cm), b, K_sat (cm h⁻¹), ID_USDA, ID_Bonan, name
    0.482 -40.5 11.4 0.46;  # 1, 11, Clay
    0.492 -49.0 10.4 0.37;  # 2, 10, Silty clay
    0.426 -15.3 10.4 0.78;  # 3, 9, Sandy clay
    0.476 -63.0 8.52 0.88;  # 4, 8, Clay loam
    0.477 -35.6 7.75 0.61;  # 5, 7, Silty clay loam
    0.420 -29.9 7.12 2.27;  # 6, 6, Sandy clay loam
    0.451 -47.8 5.39 2.50;  # 7, 5, Loam
    0.485 -78.6 5.30 2.59;  # 8, 4, Silty loam
    0.435 -21.8 4.90 12.48; # 9, 3, Sandy loam
    0.485 -78.6 5.30 2.59;  # 10, NA, Silty, nodata, used Silty loam (4)
    0.410 -9.0 4.38 56.28;  # 11, 2, Loamy sand
    0.395 -12.1 4.05 63.36  # 12, 1, Sand
  ]
  verbose && println("[Soil Texture]: ", USDA.SoilTexture(soil_type))
  θ_sat, ψ_sat, b, K_sat = soilparam[soil_type, :]
  Campbell(; θ_sat, ψ_sat, b, K_sat)
end

get_soilpar(retention::Symbol, soil_type::Int=1; kwargs...) =
    get_soilpar(Val{retention}, soil_type; kwargs...)

# function get_soilpar(::Type{Campbell}, theta::AbstractVector)
#   θ_sat, ψ_sat, b, K_sat = theta[1:4]
#   Campbell(; θ_sat, ψ_sat, b, K_sat)
# end

# function get_soilpar(theta::AbstractVector; method_retention::String="van_Genuchten")
#   method_retention == "van_Genuchten" && (P = VanGenuchten)
#   method_retention == "Campbell" && (P = Campbell)
#   get_soilpar(P, theta)
# end

export get_soilpar
