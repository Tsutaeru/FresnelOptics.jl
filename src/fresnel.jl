# f: linear frequency with standard unit (THz)
# c: speed of light, 2.99792458e-4 (m * THz)

export transmit, transmit!, reflect, reflect!, propagate, propagate!, ϵr2nκ!

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

propagate(d::Real, f::Real; forward::Bool=true) = forward ? (1.0, 1.00027 * d * f * 𝚷 / 2.99792458e-4) : (1.0, -1.00027 * d * f * 𝚷 / 2.99792458e-4)

function propagate(n::Real, κ::Real, f::Real, d::Real)
    dk = d * f * 𝚷 / 2.99792458e-4
    return exp(-κ * dk), n * dk
end

for fname in (:transmit, :reflect)
    fname! = Symbol(fname, '!')
    @eval begin
        function $(fname!)(rs::VecIO, θs::VecIO, n::VecI, κ::VecI; forward::Bool=true)
            @inbounds for i in eachindex(rs)
                rs[i], θs[i] = $(fname)(n[i], κ[i]; forward)
            end
        end
        function $(fname!)(rs::VecIO, θs::VecIO, n1::VecI, κ1::VecI, n2::VecI, κ2::VecI)
            @inbounds for i in eachindex(rs)
                rs[i], θs[i] = $(fname)(n1[i], κ1[i], n2[i], κ2[i])
            end
        end
    end
end

function propagate!(rs::VecIO, θs::VecIO, fs::VecI, d::Real; forward::Bool=true)
    ϕ = forward ? 1.00027 * d * 𝚷 : -1.00027 * d * 𝚷
    @simd for i in eachindex(θs)
        @inbounds rs[i], θs[i] = 1.0, ϕ * fs[i] / 2.99792458e-4
    end
end

function propagate!(rs::VecIO, θs::VecIO, n::VecI, κ::VecI, fs::VecI, d::Real)
    @inbounds for i in eachindex(rs)
        tmpdk = d * fs[i] * 𝚷 / 2.99792458e-4
        rs[i] = exp(-κ[i] * tmpdk)
        θs[i] = n[i] * tmpdk
    end
end

function ϵr2nκ!(ϵ1::VecIO, ϵ2::VecIO)
    @inbounds for i in eachindex(ϵ1)
        e1 = ϵ1[i]
        e2 = ϵ2[i]
        # modulus
        em = apy2(e1, e2)
        ϵ1[i] = sqrt(0.5 * (em + e1))
        ϵ2[i] = sqrt(0.5 * (em - e1))
    end
end
