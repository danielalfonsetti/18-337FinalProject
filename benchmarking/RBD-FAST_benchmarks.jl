using Profile, BenchmarkTools, Plots, Random
Random.seed!(123)

include("../src/RBD-FAST.jl")
include("../src/util.jl")

#############################
## Testing for correctness ##
#############################

# Recreate table 1 from S. Tarantola et. al (2006)
let
        sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
        sobol_function1 = create_function_of_sobol(sobol_function_params1)

        frequencies = [11, 35]
        group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]

        res1000 = RBD_FAST(sobol_function1, 1000, group_assignments, frequencies)
        println("N=1000: ", round.(res1000, digits=3), "\n")

        res10000 = RBD_FAST(sobol_function1, 10000, group_assignments, frequencies)
        println("N=10000: ", round.(res10000, digits=3), "\n")
end


##
# Recreate table 2 from S. Tarantola et. al (2006)
let
        sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
        sobol_function1 = create_function_of_sobol(sobol_function_params1)

        frequencies = [11, 35]
        group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]
        num_experiments = 50

        results_samplesize500 = zeros(num_experiments, 8)
        results_samplesize1000 = zeros(num_experiments, 8)
        results_samplesize2000 = zeros(num_experiments, 8)

        for i in 1:num_experiments
            results_samplesize500[i,:] = RBD_FAST(sobol_function1, 500, group_assignments, frequencies)
            results_samplesize1000[i,:] = RBD_FAST(sobol_function1, 1000, group_assignments, frequencies)
            results_samplesize2000[i,:] = RBD_FAST(sobol_function1, 2000, group_assignments, frequencies)
        end

        standard_devs_ss500 = mapslices(Statistics.std, results_samplesize500, dims=1)
        means_ss500 = mapslices(Statistics.mean, results_samplesize500, dims=1)
        println("Means (ss=500): ", round.(means_ss500, digits=3))
        println("Standard Dev. (ss=500): ", round.(standard_devs_ss500, digits=3), "\n")

        standard_devs_ss1000 = mapslices(Statistics.std, results_samplesize1000, dims=1)
        means_ss1000 = mapslices(Statistics.mean, results_samplesize1000, dims=1)
        println("Means (ss=1000): ", round.(means_ss1000, digits=3))
        println("Standard Dev. (ss=1000): ", round.(standard_devs_ss1000, digits=3), "\n")

        standard_devs_ss2000 = mapslices(Statistics.std, results_samplesize2000, dims=1)
        means_ss2000 = mapslices(Statistics.mean, results_samplesize2000, dims=1)
        println("Means (ss=2000): ", round.(means_ss2000, digits=3))
        println("Standard Dev. (ss=2000): ", round.(standard_devs_ss2000, digits=3), "\n")
end

##
# Recreate figure 2a from S. Tarantola et. al (2006)

let
        sobol_function_params_2a = cat([99 for i in 1:10], [0 for i in 1:10], dims=1)
        sobol_function2a = create_function_of_sobol(sobol_function_params_2a)

        frequencies = [11, 21, 27, 35, 39]
        group_assignments =  cat([1 for i in 1:4], [2 for i in 1:4], [3 for i in 1:4], [4 for i in 1:4], [5 for i in 1:4], dims=1)


        fig2a_res5000= RBD_FAST(sobol_function2a, 5000, group_assignments, frequencies)
        scatter(fig2a_res5000, legend=nothing,  title="RBD-FAST on Sobol Function \n (10 non-important, 10 important factors)")
        hline!([0.0199])
        xaxis!("Model Factors")
        savefig("./18337FinalProject/plots/rbd_fast_fig2a.png")
end

##
# Recreate figure 2b from S. Tarantola et. al (2006)

let
        sobol_function_params_2b = cat([99 for i in 1:5], [0 for i in 1:15], dims=1)
        sobol_function2b = create_function_of_sobol(sobol_function_params_2b)

        frequencies = [1]
        group_assignments =  cat([1 for i in 1:20], dims=1)

        fig2a_res5000= RBD_FAST(sobol_function2b, 10000, group_assignments, frequencies)
        scatter(fig2a_res5000, legend=nothing, title="RBD-FAST on Sobol Function \n (5 non-important, 15 important factors)")
        hline!([0.0045])
        xaxis!("Model Factors")
        savefig("./18337FinalProject/plots/rbd_fast_fig2b.png")

end
## Verify convergence with increasing sample sizes

let
        sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
        sobol_function1 = create_function_of_sobol(sobol_function_params1)

        frequencies = [11, 35]
        group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]

        true_effects = [0.7160, 0.1790, 0.0240, 0.0072, 0.0001, 0.0001, 0.0001, 0.0001]
        sample_sizes = [i*100 for i in 5:20]
        divergences = zeros(length(sample_sizes))
        diverges_std = zeros(length(sample_sizes))
        num_experiments = 1000

        divergences_for_sample_size = zeros(num_experiments)
        for (index, sample_size) in enumerate(sample_sizes)
             println(index, " out of ", length(sample_sizes))

             for i in 1:num_experiments
                     estimated_effects = RBD_FAST(sobol_function1, sample_size, group_assignments, frequencies)
                     divergence =  sqrt(sum((estimated_effects .- true_effects).^2))
                     divergences_for_sample_size[i] = divergence
             end
            divergences[index] = mean(divergences_for_sample_size)
        end
        plot(sample_sizes, divergences,  markershape=:circle, legend=nothing, title=string("Average Euclidean Distance from Analytical Value with \n Increasing Sample Size \n(", num_experiments, " experiments per sample)"))
        # scatter!(sample_sizes, divergences)
        xaxis!("Sample Size")
        yaxis!("Euclidean Distance")
        savefig("./18337FinalProject/plots/rbd_fast_sample_size_increase.png")
end


##################
## Benchmarking ##
##################

# Compare parallelized/multi threaded RBD-FAST to non-parallelized RBD-FAST

let
        sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
        sobol_function1 = create_function_of_sobol(sobol_function_params1)

        frequencies = [11, 35]
        group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]

        sample_sizes = 2 .^ [i for i in 9:14]

        non_multi_threaded_times = zeros(length(sample_sizes))
        non_multi_threaded_allocs = zeros(length(sample_sizes))
        non_multi_threaded_memory = zeros(length(sample_sizes))
        multi_threaded_times = zeros(length(sample_sizes))
        multi_threaded_allocs = zeros(length(sample_sizes))
        multi_threaded_memory = zeros(length(sample_sizes))

        for (index, sample_size) in enumerate(sample_sizes)
                println(index, " out of ", length(sample_sizes))
                non_multi_threaded_benchmark = @benchmark _RBD_FAST_non_multi($sobol_function1, $sample_size, $group_assignments, $frequencies)
                multi_threaded_benchmark = @benchmark RBD_FAST($sobol_function1, $sample_size, $group_assignments, $frequencies)

                non_multi_threaded_times[index] = time(non_multi_threaded_benchmark)
                non_multi_threaded_allocs[index] = allocs(non_multi_threaded_benchmark)
                non_multi_threaded_memory[index] = memory(non_multi_threaded_benchmark)
                multi_threaded_times[index] = time(multi_threaded_benchmark)
                multi_threaded_allocs[index] = allocs(multi_threaded_benchmark)
                multi_threaded_memory[index] = memory(multi_threaded_benchmark)
        end

        plot(
                sample_sizes,
                non_multi_threaded_times,
                legend=:topleft,
                markershape=:circle,
                label="Non-multithreaded",
                xaxis=:log,
                xticks = (sample_sizes, sample_sizes),
                title="RBD-FAST Time Benchmarks \n on Sobol Function"
        )
        yaxis!("Time (ns)")
        xaxis!("Number of Samples")
        plot!(sample_sizes, multi_threaded_times, markershape=:circle, label="Multithreaded")
        savefig("./18337FinalProject/plots/rbd_fast_times.png")

        plot(
                sample_sizes,
                non_multi_threaded_memory,
                legend=:topleft,
                markershape=:circle,
                label="Non-multithreaded",
                xaxis=:log,
                xticks = (sample_sizes, sample_sizes),
                title="RBD-FAST Memory Benchmarks \n on Sobol Function"
        )
        yaxis!("Memory (MiB)")
        xaxis!("Number of Samples")
        plot!(sample_sizes, multi_threaded_memory, markershape=:circle, label="Multithreaded")
        savefig("./18337FinalProject/plots/rbd_fast_memory.png")


        plot(
                sample_sizes,
                non_multi_threaded_allocs,
                legend=:topleft,
                markershape=:circle,
                label="Non-multithreaded",
                xaxis=:log,
                xticks = (sample_sizes, sample_sizes),
                title="RBD-FAST Allocation Benchmarks \n on Sobol Function"
        )
        yaxis!("Number of Allocations")
        xaxis!("Number of Samples")
        plot!(sample_sizes, multi_threaded_allocs, markershape=:circle, label="Multithreaded")
        savefig("./18337FinalProject/plots/rbd_fast_allocations.png")
end

## Flame graph

let
        sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
        sobol_function1 = create_function_of_sobol(sobol_function_params1)

        frequencies = [11, 35]
        group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]

        @profile for i in 1:10 RBD_FAST(sobol_function1, 500, group_assignments, frequencies) end
        Juno.profiler()
        Profile.clear()
end
