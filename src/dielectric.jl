export ϵr2nκ!

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
