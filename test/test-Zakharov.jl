using ModelParams, Test, Parameters

# Zakharov: 最优值为1
function Zakharov(x::AbstractVector{Float64})
    y1 = 0.0
    y2 = 0.0
    @inbounds for i in eachindex(x)
        xi = x[i] - 0.5
        y1 += xi * xi
        y2 += (i + 1) * xi
    end
    y2_half = 0.5 * y2
    return y1 + y2_half^2 + y2_half^4 + 1
end


# SCEUA: 无法胜任当前任务

# @testset "Zakharov" 
begin
    n = 100
    lower = -1.0 * ones(n)
    upper = 1.0 * ones(n)
    x0 = zeros(n)
    @time u, feval, retcode = sceua(Zakharov, x0, lower, upper;
        f_reltol=1e-10, 
        maxn=Int(1e4), verbose=true)
end
