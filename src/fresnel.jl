transmit(n::Real, κ::Real; forward::Bool=true) = forward ? transmit(1.00027, 0.0, n, κ) : transmit(n, κ, 1.00027, 0.0)

function transmit(ni::Real, κi::Real, nj::Real, κj::Real)
    nipnj = ni + nj
    κipκj = κi + κj

    m = 2.0 * apy2(ni, κi) / apy2(nipnj, κipκj)
    θ = atan(κi * nipnj - ni * κipκj, ni * nipnj + κi * κipκj)
    return m, θ
end

reflect(n::Real, κ::Real; forward::Bool=true) = forward ? reflect(1.00027, 0.0, n, κ) : reflect(n, κ, 1.00027, 0.0)

function reflect(ni::Real, κi::Real, nj::Real, κj::Real)
    nipnj = ni + nj
    κipκj = κi + κj
    nimnj = ni - nj
    κimκj = κi - κj

    x = nimnj * nipnj + κimκj * κipκj
    y = nipnj * κimκj - nimnj * κipκj

    m = flipsign(apy2(nimnj, κimκj) / apy2(nipnj, κipκj), x)
    iszero(y) && return m, 0.0
    return m, atan(abs(y), abs(x))
end

propagate(d::Real, k::Real) = 1.0, 1.00027 * d * k

function propagate(n::Real, κ::Real, d::Real, k::Real)
    dk = d * k
    return exp(-κ * dk), n * dk
end
