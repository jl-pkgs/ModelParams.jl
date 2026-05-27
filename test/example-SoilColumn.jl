using ModelParams, Parameters

FT = Float64
N = 5

# 构造默认 SoilColumn（Campbell 水力 + KvLayers + ThermalMain）
ps = SoilColumn{FT,N}()

# ── 1. 查看所有参数（30 行：15 Campbell + 5 KvLayers + 10 ThermalMain）
@time parameters(ps)

# ── 2. 仅提取水力参数用于率定（20 行）
params = filter_params(ps, :hydraulic)

# ── 3. ψ_sat 各层共享（Campbell b 排除）：16 行
params = filter_params(ps, :hydraulic;
    list_sameLayer = [:ψ_sat],   # ψ_sat 只保留第 1 层，update 时广播到所有层
    list_fix       = [:b])       # b 固定，不参与率定

# ── 4. 获取初值 / 边界（供优化器使用）
x0 = Float64.(params.value)
lb  = Float64[b[1] for b in params.bound]
ub  = Float64[b[2] for b in params.bound]

# ── 5. 用新参数更新模型（优化器迭代时调用）
x_new = x0 .* 1.05                        # 示例：把所有参数放大 5%
update_params!(ps, params.path, x_new;
    params,
    list_sameLayer = [:ψ_sat],
    list_fix       = [:b])

# 验证：ψ_sat 已广播到所有层
@show ps.hydraulic.profile.ψ_sat          # 各层应相同
@show ps.hydraulic.layers[1].ψ_sat        # AoS 缓存同步更新
