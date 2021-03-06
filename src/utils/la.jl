function inf_geometric_sum(r::Real, θ::Real)
    x = 1.0 - r * cos(θ)
    y = -r * sin(θ)
    return inv(apy2(x, y)), -atan(y, x)
end

function sqr2(x::Real, y::Real)
    isnan(x) && return x
    isnan(y) && return y
    # general case
    xabs = abs(x)
    yabs = abs(y)
    w = max(xabs, yabs)
    z = min(xabs, yabs)
    iszero(z) && return abs2(w)
    return abs2(w) * (1.0 + abs2(z / w))
end

function apy2(x::Real, y::Real)
    isnan(x) && return x
    isnan(y) && return y
    # general case
    xabs = abs(x)
    yabs = abs(y)
    w = max(xabs, yabs)
    z = min(xabs, yabs)
    iszero(z) && return w
    return w * sqrt(1.0 + abs2(z / w))
end
