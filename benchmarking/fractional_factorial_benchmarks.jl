using BenchmarkTools, Plots, Profile
include("../src/fractional_factorial.jl")
include("../src/util.jl")


#############################
## Testing for correctness ##
#############################

# See if fractional factorial can estimate variances in sobol function.

sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(sobol_function_params1)

design_matrix = generate_ff_design_matrix(length(sobol_function_params1))
sample_matrix = generate_ff_sample_matrix(design_matrix, [(0.1, 0.8) for i in 1:8])

response_vec = run_model(sample_matrix, sobol_function1)

main_effects = ff_main_effects(design_matrix, response_vec)
println(round.(main_effects, digits=3))
variances = ff_variance(design_matrix, response_vec)
println(variances)

# Estimates and variances do follow the order of parameter importance, as expected.



## Flame graph
@profile for i in 1:1000 _recursive_hadamard(1024) end
Juno.profiler()
Profile.clear()

@profile for i in 1:1000 _generate_hadamard(1024) end
Juno.profiler()
Profile.clear()

##################
## Benchmarking ##
##################

ks = 2 .^ (1:10)
recursive_method_times = zeros(length(ks))
recursive_method_memory =  zeros(length(ks))
recursive_method_allocs =  zeros(length(ks))

expanding_window_method_times = zeros(length(ks))
expanding_window_method_memory = zeros(length(ks))
expanding_window_method_allocs = zeros(length(ks))

for (index, k) in enumerate(ks)
        println(index, " out of ", length(ks))
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
        markershape=:circle,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Time)"
)
yaxis!("Time (ns)")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_times, markershape=:circle, label="Expanding Window Method")
savefig("./18337FinalProject/plots/hadamard_time.png")


scatter(
        ks,
        recursive_method_memory,
        label="Recursive Method",
        xaxis=:log,
        markershape=:circle,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Memory)"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_memory, markershape=:circle, label="Expanding Window Method")
savefig("./18337FinalProject/plots/hadamard_memory.png")


plot(
        ks,
        recursive_method_allocs,
        label="Recursive Method",
        xaxis=:log,
        markershape=:circle,
        xticks = (ks, ks),
        legend=:topleft,
        title="Hadamard Generation \n Method Comparison (Allocations)"
)
yaxis!("Number of Allocations")
xaxis!("Number of Rows")
plot!(ks, expanding_window_method_allocs, markershape=:circle, label="Expanding Window Method")
savefig("./18337FinalProject/plots/hadamard_allocations.png")
