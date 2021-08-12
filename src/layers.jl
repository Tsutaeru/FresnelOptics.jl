export Layer

# For a linked-list type layred structure in the future.
mutable struct Layer
    power::SubArray{Float64, 2, Array{Float64,3}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true}
    phase::SubArray{Float64, 2, Array{Float64,3}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true}
    core::Array{Float64,3}
    prev::Layer
    next::Layer

    function Layer(sz::Int)
        # t_af, t_fa, r_af, r_fa, pr_f, pr_a
        core = Array{Float64,3}(undef, sz, 6, 2)
        node = new(view(core, :, :, 1), view(core, :, :, 2), core)
        node.prev = node
        node.next = node
        return node
    end
end

function setup!(o::Layer, n::VecI, κ::VecI, fs::VecI, d::Real)
    power = o.power
    phase = o.phase
    #### Forward/Backward Transmission
    transmit!(view(power, :, 1), view(phase, :, 1), n, κ; forward=true)
    transmit!(view(power, :, 2), view(phase, :, 2), n, κ; forward=false)
    #### Forward/Backward Reflection
    reflect!(view(power, :, 3), view(phase, :, 3), n, κ; forward=true)
    reflect!(view(power, :, 4), view(phase, :, 4), n, κ; forward=false)
    #### Propagation in layer/air
    propagate!(view(power, :, 5), view(phase, :, 5), n, κ, fs, d)
    propagate!(view(power, :, 6), view(phase, :, 6), fs, d)
end

function apply!(rs::VecIO, θs::VecIO, o::Layer)
    power = o.power
    phase = o.phase
    one2n = eachindex(rs)
    @simd for i in one2n
        @inbounds rs[i] = power[i,1] * power[i,2] * power[i,5] / power[i,6]
    end
    @simd for i in one2n
        @inbounds θs[i] = phase[i,1] + phase[i,2] + phase[i,5] - phase[i,6]
    end
    @inbounds for i in one2n
        tmpr = power[i,4] * power[i,5]
        tmpθ = phase[i,4] + phase[i,5]
        tmpr, tmpθ = inf_geometric_sum(tmpr * tmpr, 2.0 * tmpθ)
        rs[i] *= tmpr
        θs[i] += tmpθ
    end
end
