function _recursive_hadamard(k::Integer)
        """
        Generate a hadamard matrix via recursion.
        """
        # base case
        if k == 2
                return [1 1; 1 -1]
        else
                h = _recursive_hadamard(k ÷ 2)
                return hcat(vcat(h, h), vcat(h, -h))
        end
end


function _expanding_window_hadamard(k::Integer)
        """
        Generate hadamard matrix of size k using expanding window approach.
        """
        @assert ispow2(k)

        # intialize
        h = ones(Int64, k, k)
        h[2,2] = -1

        let
                bot_row = 2
                right_col = 2
                while bot_row < k

                        @inbounds cur = h[1:bot_row, 1:right_col]

                        new_bot = bot_row+1:bot_row*2
                        new_right = right_col+1:right_col*2
                        # update right
                        @inbounds h[1:bot_row, new_right] = cur

                        # update below
                        @inbounds h[new_bot, 1:right_col] = cur

                        # update diagonal
                        @inbounds h[new_bot, new_right] = - cur

                        # update window for each to 'copy' from.
                        bot_row *= 2
                        right_col *= 2
                end
        end
        return h
end


function _generate_hadamard(k)
        @assert ispow2(k)
        return _expanding_window_hadamard(k)
end


function resolution_iv_design(k::Integer)
        s = _generate_hadamard(k)
        return vcat(s, -s)
end


function ff_sample(num_parameters, levels_list=nothing)
        """
        Returns a 2-level fractional factorial design matrix, where
        number of rows is 2*number of parameters.

        See pg.74 in Global Sensitivity Analysis for Details.
        """

        # If k is not a power of 2, get the next larger number that is.
        k = Integer(round(2^ceil(log2(num_parameters))))

        design_matrix = resolution_iv_design(k)
        # Remove dummy columns
        design_matrix = design_matrix[:, 1:num_parameters]


        if levels_list == nothing
                design_matrix[design_matrix.==-1] .= 0
        else
                for (col_index, (low_value, high_value)) in enumerate(levels_list)
                        col = @view design_matrix[:,col_index] # get pointer to the column we want to update.
                        col[col .== -1] .= low_value
                        col[col .== 1] .= high_value
                end
        end
        return design_matrix
end


function ff_main_effect(contrast_vector, response_values)
        return dot(contrast_vector, response_values)/(length(contrast_vector))
end

function ff_variance(contrast_vector, response_values)
        return main_effect(contrast_vector,response_values)^2
end

using BenchmarkTools, Plots, Profile

@profile for i in 1:1000 _recursive_hadamard(1024) end
Juno.profiler()
Profile.clear()

@profile for i in 1:1000 _generate_hadamard(1024) end
Juno.profiler()
Profile.clear()


recursive_method_times = zeros(length(ks))
recursive_method_memory =  zeros(length(ks))
recursive_method_allocs =  zeros(length(ks))

expanding_window_method_times = zeros(length(ks))
expanding_window_method_memory = zeros(length(ks))
expanding_window_method_allocs = zeros(length(ks))

ks = 2 .^ (1:10)
for (index, k) in enumerate(ks)
        x = @benchmark _recursive_hadamard(Integer($k))
        recursive_method_times[index] = time(x)
        recursive_method_memory[index] = memory(x)
        recursive_method_allocs[index] = allocs(x)
        y = @benchmark _expanding_window_hadamard(Integer($k))
        expanding_window_method_times[index] = time(y)
        expanding_window_method_memory[index] = memory(y)
        expanding_window_method_allocs[index] = allocs(y)
end


plot(
        ks,
        recursive_method_times,
        label="Recursive Method",
        xaxis=:log,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Time)"
)
yaxis!("Time (ns)")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_times, label="Expanding Window Method")


plot(
        ks,
        recursive_method_memory,
        label="Recursive Method",
        xaxis=:log,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Memory)"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_memory, label="Expanding Window Method")


plot(
        ks,
        recursive_method_allocs,
        label="Recursive Method",
        xaxis=:log,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Allocations)"
)
yaxis!("Number of Allocations")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_allocs, label="Expanding Window Method")


ff_sample(5, [(-1, 1), (-1, 1), (-2, 5), (-1, 1), (-1, 1)])
