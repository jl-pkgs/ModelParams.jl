# TODO: 这里使用枚举型更合适

module USDA

export soil_texture, SoilTexture
export CLAY, SILTY_CLAY, SANDY_CLAY, CLAY_LOAM, SILTY_CLAY_LOAM,
  SANDY_CLAY_LOAM, LOAM, SILTY_LOAM, SANDY_LOAM, SILT,
  LOAMY_SAND, SAND

# https://developers.google.com/earth-engine/datasets/catalog/OpenLandMap_SOL_SOL_TEXTURE-CLASS_USDA-TT_M_v02
@enum SoilTexture begin
  # USDA_UNDEFINED  # NA  
  CLAY = 1 # CL
  SILTY_CLAY = 2 # SICL
  SANDY_CLAY = 3 # SACL
  CLAY_LOAM = 4 # CLLO
  SILTY_CLAY_LOAM = 5 # SICLLO
  SANDY_CLAY_LOAM = 6 # SACLLO
  LOAM = 7 # LO
  SILTY_LOAM = 8 # SILO
  SANDY_LOAM = 9 # SALO
  SILT = 10 # SI
  LOAMY_SAND = 11 # LOSA
  SAND = 12 # SA
end

# const SoilTypes = [
#   "CLAY", "SILTY_CLAY", "SANDY_CLAY", "CLAY_LOAM", "SILTY_CLAY_LOAM",
#   "SANDY_CLAY_LOAM", "LOAM", "SILTY_LOAM", "SANDY_LOAM", "SILT",
#   "LOAMY_SAND", "SAND"]
# soil_texture(i::Int) = SoilTypes[I]

Base.to_index(i::SoilTexture) = Int(i)

# 土壤类型名称映射（英文 → 枚举）
const SOIL_TEXTURE_NAMES = Dict{String,SoilTexture}(
  "CLAY" => CLAY,
  "SILTY_CLAY" => SILTY_CLAY,
  "SANDY_CLAY" => SANDY_CLAY,
  "CLAY_LOAM" => CLAY_LOAM,
  "SILTY_CLAY_LOAM" => SILTY_CLAY_LOAM,
  "SANDY_CLAY_LOAM" => SANDY_CLAY_LOAM,
  "LOAM" => LOAM,
  "SILTY_LOAM" => SILTY_LOAM,
  "SANDY_LOAM" => SANDY_LOAM,
  "SILT" => SILT,
  "LOAMY_SAND" => LOAMY_SAND,
  "SAND" => SAND,
)

# 土壤类型名称映射（中文 → 枚举）
const SOIL_TEXTURE_NAMES_CN = Dict{String,SoilTexture}(
  "黏土" => CLAY,
  "粉砂质黏土" => SILTY_CLAY,
  "砂质黏土" => SANDY_CLAY,
  "黏壤土" => CLAY_LOAM,
  "粉砂质黏壤土" => SILTY_CLAY_LOAM,
  "砂质黏壤土" => SANDY_CLAY_LOAM,
  "壤土" => LOAM,
  "粉砂质壤土" => SILTY_LOAM,
  "砂质壤土" => SANDY_LOAM,
  "粉砂土" => SILT,
  "壤质砂土" => LOAMY_SAND,
  "砂土" => SAND,
)

# 数值 → 枚举的反向查找
const SOIL_TEXTURE_BY_ID = Dict{Int,SoilTexture}(
  Int(CLAY) => CLAY,
  Int(SILTY_CLAY) => SILTY_CLAY,
  Int(SANDY_CLAY) => SANDY_CLAY,
  Int(CLAY_LOAM) => CLAY_LOAM,
  Int(SILTY_CLAY_LOAM) => SILTY_CLAY_LOAM,
  Int(SANDY_CLAY_LOAM) => SANDY_CLAY_LOAM,
  Int(LOAM) => LOAM,
  Int(SILTY_LOAM) => SILTY_LOAM,
  Int(SANDY_LOAM) => SANDY_LOAM,
  Int(SILT) => SILT,
  Int(LOAMY_SAND) => LOAMY_SAND,
  Int(SAND) => SAND,
)

"""
    soil_type_from_name(name::Union{String,Symbol}) -> SoilTexture

根据名称获取土壤类型枚举。支持英文和中文名称。

# Arguments
- `name`: 土壤类型名称（如 "LOAM", "壤土", "SANDY_LOAM"）

# Returns
- `SoilTexture`: 对应的枚举值

# Examples
```julia
soil_type_from_name("LOAM")        # → LOAM (7)
soil_type_from_name("壤土")         # → LOAM (7)
soil_type_from_name(:SILTY_CLAY)   # → SILTY_CLAY (2)
```
"""
function soil_type_from_name(name::Union{String,Symbol})
  s = uppercase(string(name))

  # 尝试英文名称
  if haskey(SOIL_TEXTURE_NAMES, s)
    return SOIL_TEXTURE_NAMES[s]
  end

  # 尝试中文名称
  if haskey(SOIL_TEXTURE_NAMES_CN, name)
    return SOIL_TEXTURE_NAMES_CN[name]
  end

  error("Unknown soil texture name: '$name'. " *
        "Valid names: $(keys(SOIL_TEXTURE_NAMES)) or $(keys(SOIL_TEXTURE_NAMES_CN))")
end

"""
    soil_type_from_id(id::Int) -> SoilTexture

根据数值 ID 获取土壤类型枚举。

# Examples
```julia
soil_type_from_id(7)   # → LOAM
soil_type_from_id(1)   # → CLAY
```
"""
function soil_type_from_id(id::Int)
  haskey(SOIL_TEXTURE_BY_ID, id) || error("Unknown soil type ID: $id (valid: 1-12)")
  return SOIL_TEXTURE_BY_ID[id]
end

"""
    parse_soil_type(x::Union{Int,String,Symbol}) -> SoilTexture

通用土壤类型解析函数。支持数值、字符串、符号输入。

# Examples
```julia
parse_soil_type(7)            # → LOAM
parse_soil_type("LOAM")       # → LOAM
parse_soil_type("壤土")        # → LOAM
parse_soil_type(:SILTY_LOAM)  # → SILTY_LOAM
```
"""
parse_soil_type(x::Int) = soil_type_from_id(x)
parse_soil_type(x::Union{String,Symbol}) = soil_type_from_name(x)

# 导出新增的函数
export soil_type_from_name, soil_type_from_id, parse_soil_type
export SOIL_TEXTURE_NAMES, SOIL_TEXTURE_NAMES_CN

# https://github.com/NigelVanNieuwenhuizen/USDA-Soil-Texture-Calculator/blob/9cbf7ee58d384074de26e9a50fa4d1872e9a37f8/USDASoilTextureCalculator.py#L261C17-L286C46

"""
    soil_texture(sand::Float64, silt::Float64)

# Arguments
- `sand`: percentage of sand, %
- `silt`: percentage of silt, %

# Example usage:
```
sand = 60.0
silt = 20.0
println("Soil Texture: ", soil_texture(sand, silt))
```
"""
function soil_texture(sand::Float64, silt::Float64)
  clay = 100 - sand - silt

  if sand <= 45 && silt <= 40 && clay >= 40
    CLAY
  elseif sand <= 65 && sand >= 45 && silt <= 20 && clay >= 35 && clay <= 55
    SANDY_CLAY
  elseif sand <= 20 && silt >= 40 && silt <= 60 && clay >= 40 && clay <= 60
    SILTY_CLAY
  elseif sand >= 45 && sand <= 80 && silt <= 28 && clay >= 20 && clay <= 35
    SANDY_CLAY_LOAM
  elseif sand >= 20 && sand <= 45 && silt >= 15 && silt <= 53 && clay >= 27 && clay <= 40
    CLAY_LOAM
  elseif sand <= 20 && silt >= 40 && silt <= 73 && clay >= 27 && clay <= 40
    SILTY_CLAY_LOAM
  elseif sand >= 43 && sand <= 85 && silt <= 50 && clay <= 20
    SANDY_LOAM
  elseif sand >= 23 && sand <= 52 && silt >= 28 && silt <= 50 && clay >= 7 && clay <= 27
    LOAM
  elseif sand <= 50 && silt >= 50 && silt <= 88 && clay <= 27
    SILTY_LOAM
  elseif sand <= 20 && silt >= 80 && clay <= 12
    SILT
  elseif sand >= 70 && sand <= 90 && silt <= 30 && clay <= 15
    LOAMY_SAND
  elseif sand >= 85 && silt <= 15 && clay <= 10
    SAND
  else
    @error("Soil texture not found for sand: $sand, silt: $silt, clay: $clay")
  end
end


end

export USDA
