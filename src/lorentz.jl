export lorentz_factors, LorentzModel, setup!, apply!

# Deprecated, will be moved to Osnaps.jl in the future.
lorentz_factors(::Type{T}; ϵc, ωp, fs, Ωs, γs) where T = lorentz_factors(T, ϵc, ωp, fs, Ωs, γs)

function lorentz_factors(::Type{Vector}, ϵc::Real, ωp::Real, fs::NTuple{N}, Ωs::NTuple{N}, γs::NTuple{N}) where N
    if @generated
        a = Vector{Any}(undef, 3*N+2)
        @inbounds a[1] = :ϵc
        @inbounds a[2] = :ωp
        @inbounds for i in 1:N
            a[i+2]     = :(fs[$i])
            a[i+2+  N] = :(Ωs[$i])
            a[i+2+2*N] = :(γs[$i])
        end
        return quote
            $(Expr(:meta, :inline))
            @inbounds return $(Expr(:vect, a...))
        end
    else
        a = Vector{Float64}(undef, 3*N+2)
        @inbounds a[1] = ϵc
        @inbounds a[2] = ωp
        @inbounds for i in 1:N
            a[i+2], a[i+2+N], a[i+2+2*N] = fs[i], Ωs[i], γs[i]
        end
        return a
    end
end

function lorentz_factors(::Type{Tuple}, ϵc::Real, ωp::Real, fs::NTuple{N}, Ωs::NTuple{N}, γs::NTuple{N}) where N
    if @generated
        a = Vector{Any}(undef, 3*N+2)
        @inbounds a[1] = :ϵc
        @inbounds a[2] = :ωp
        @inbounds for i in 1:N
            a[i+2]     = :(fs[$i])
            a[i+2+  N] = :(Ωs[$i])
            a[i+2+2*N] = :(γs[$i])
        end
        return quote
            $(Expr(:meta, :inline))
            @inbounds return $(Expr(:tuple, a...))
        end
    else
        return ntuple(i -> lorentz_factors(i, ϵc, ωp, fs, Ωs, γs), 3*N+2)
    end
end

function lorentz_factors(ix::Int, ϵc::Real, ωp::Real, fs::NTuple{N}, Ωs::NTuple{N}, γs::NTuple{N}) where N
    ix <= 0  && error("BoundsError: lorentz_factors(ix = $ix, ...)")
    ix == 1  && return ϵc
    ix == 2  && return ωp
    lx  = 2; rx = 2 + N
    ix <= rx && return fs[ix - lx]
    lx += N; rx += N
    ix <= rx && return Ωs[ix - lx]
    lx += N; rx += N
    ix <= rx && return γs[ix - lx]
    error("BoundsError: lorentz_factors(ix = $ix, ...)")
end

struct LorentzModel
    re_kernel::SubArray{Float64, 2, Array{Float64,3}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true}
    im_kernel::SubArray{Float64, 2, Array{Float64,3}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true}
    re_params::SubArray{Float64, 1, Array{Float64,2}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}
    im_params::SubArray{Float64, 1, Array{Float64,2}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}
    kernel::Array{Float64,3}
    params::Array{Float64,2}
    ratios::Array{Float64,1}
    nExp::Int
    nLor::Int

    function LorentzModel(nExp::Int, nLor::Int)
        kernel = Array{Float64,3}(undef, nExp, nLor, 2)
        params = Array{Float64,2}(undef, nExp, 2)
        ratios = Array{Float64,1}(undef, nLor)
        return new(view(kernel, :, :, 1), view(kernel, :, :, 2), view(params, :, 1), view(params, :, 2), kernel, params, ratios, nExp, nLor)
    end
end

function setup!(obj::LorentzModel, fs::VecI)
    re_params = obj.re_params
    im_params = obj.im_params
    @simd for i in eachindex(1:obj.nExp)
        @inbounds re_params[i], im_params[i] = abs2(fs[i]), fs[i]
    end
    return nothing
end

function setup!(obj::LorentzModel, fs::VecI, Ωs::VecI, γs::VecI)
    one2nExp  = eachindex(1:obj.nExp)
    ratios    = obj.ratios
    re_params = obj.re_params
    im_params = obj.im_params
    re_kernel = obj.re_kernel
    im_kernel = obj.im_kernel

    for j in 1:obj.nLor
        @inbounds Ωj2 = abs2(Ωs[j])
        @inbounds γj  = γs[j]
        @simd for i in one2nExp
            @inbounds re_kernel[i,j], im_kernel[i,j] = Ωj2 - re_params[i], γj * im_params[i]
        end
        @inbounds for i in one2nExp
            denom = sqr2(re_kernel[i,j], im_kernel[i,j])
            re_kernel[i,j] /= denom
            im_kernel[i,j] /= denom
        end
        @inbounds ratios[j] = fs[j]
    end
    return nothing
end

function apply!(ϵ1::VecIO, ϵ2::VecIO, obj::LorentzModel, ϵc::Real, ωp::Real)
    @simd for i in eachindex(1:obj.nExp)
        @inbounds ϵ1[i], ϵ2[i] = ϵc, 0.0
    end

    ωp2 = ωp * ωp
    gemv!('N', ωp2, obj.re_kernel, obj.ratios, true, ϵ1)
    gemv!('N', ωp2, obj.im_kernel, obj.ratios, true, ϵ2)
    return nothing
end
