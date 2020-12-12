function generate_hadamard(k::Integer)
        # pg 73
        # k must be a power of 2
        @assert ispow2(k)

        # base case
        if k == 2
                return [1 1; 1 -1] # return smallest Hadamard matrix
        else
                h = generate_hadamard(k ÷ 2)
                return hcat(vcat(h, h), vcat(h, -h))
        end
end

generate_hadamard(Integer(8))

function resolution_iv_design(k::Integer)
        # pg 74
        @assert ispow2(k)
        s = generate_hadamard(k)
        return vcat(s, -s)
end


# pg 75
function main_effect(contrast_vector, response_values)
        return dot(contrast_vector, response_values)/(length(contrast_vector))
end

function two_level_ff_design_variance(contrast_vector, response_values)
        return main_effect(contrast_vector,response_values)^2
end
