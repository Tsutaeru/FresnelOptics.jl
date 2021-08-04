export LorentzModel, setup!, apply!

struct LorentzModel
    resonator::Matrix{Float64}
    frequency::Vector{Float64}
    strength::Vector{Float64}
    nExp::Int
    nLor::Int

    function LorentzModel(nExp::Int, nLor::Int)
        resonator = Matrix{Float64}(undef, nExp << 1, nLor)
        frequency = Vector{Float64}(undef, nExp << 1)
        strength  = Vector{Float64}(undef, nLor)
        return new(resonator, frequency, strength, nExp, nLor)
    end
end

function setup!(obj::LorentzModel, ωs::VecI)
    nExp = obj.nExp
    freq_real = @view obj.frequency[1:nExp]
    freq_imag = @view obj.frequency[nExp+1:end]
    @simd for i in eachindex(1:nExp)
        @inbounds freq_real[i], freq_imag[i] = abs2(ωs[i]), ωs[i]
    end
    return nothing
end

function setup!(obj::LorentzModel, ωs::VecI)
    nExp = obj.nExp
    freq_real = @view obj.frequency[1:nExp]
    freq_imag = @view obj.frequency[nExp+1:end]
    @simd for i in eachindex(1:nExp)
        @inbounds freq_real[i], freq_imag[i] = abs2(ωs[i]), ωs[i]
    end
    return nothing
end

function setup!(obj::LorentzModel, fs::VecI, Ωs::VecI, γs::VecI)
    nExp = obj.nExp
    nLor = obj.nLor

    strength  = obj.strength
    freq_real = @view obj.frequency[1:nExp]
    freq_imag = @view obj.frequency[nExp+1:end]
    resr_real = @view obj.resonator[1:nExp, :]
    resr_imag = @view obj.resonator[nExp+1:end, :]

    for j in 1:nLor
        @inbounds Ωj2 = abs2(Ωs[j])
        @inbounds γj  = γs[j]
        @simd for i in eachindex(1:nExp)
            @inbounds resr_real[i,j], resr_imag[i,j] = Ωj2 - freq_real[i], γj * freq_imag[i]
        end
        @inbounds for i in eachindex(1:nExp)
            denom = abs2(resr_real[i,j]) + abs2(resr_imag[i,j])
            resr_real[i,j] /= denom
            resr_imag[i,j] /= denom
        end
        @inbounds strength[j] = fs[j]
    end    
    return nothing
end

function apply!(ϵ12::VecIO, obj::LorentzModel, ϵc::Real, ωp::Real)
    nExp = obj.nExp
    ϵ1 = @view ϵ12[1:nExp]
    ϵ2 = @view ϵ12[nExp+1:end]

    @simd for i in eachindex(1:nExp)
        @inbounds ϵ1[i], ϵ2[i] = ϵc, 0.0
    end
    gemv!('N', abs2(ωp), obj.resonator, obj.strength, true, ϵ12)
    return nothing
end
