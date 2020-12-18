using Plots, BenchmarkTools, Profile

include("../src/delta.jl")
include("../src/util.jl")

#############################
## Testing for correctness ##
#############################

# Test independence
X = hcat([i for i in 1:9], [1 for i in 1:9])
Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]
deltas, adjusted_deltas, adjusted_deltas_low, adjusted_deltas_high = delta_moment(X, Y, num_resamples=500)

###################################
# Try to recreate figure 3 from Plischke et al. (2013)

sample_sizes = 2 .^ [i for i in 9:14]
fig3_deltas = zeros(length(sample_sizes))
fig3_adjusted_deltas = zeros(length(sample_sizes))
fig3_adj_delta_95low = zeros(length(sample_sizes))
fig3_adj_delta_95high = zeros(length(sample_sizes))

for (index, sample_size) in enumerate(sample_sizes)
        println(index, " out of ", length(sample_sizes))
        X, Y = create_ishigami_samples(sample_size)
        res = delta_moment(X[:,4], Y, Ygrid_length=2048, num_classes=10)
        fig3_deltas[index], fig3_adjusted_deltas[index], fig3_adj_delta_95low[index], fig3_adj_delta_95high[index] = res[1][1], res[2][1], res[3][1], res[4][1]
end


plot(
        sample_sizes,
        fig3_deltas,
        label="Deltas",
        markershape=:circle,
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        legend=:topright,
        title="Ishigami Function Inactive Variable Sensitivity Measures"
)
yaxis!("Delta Sensitivity Measure")
xaxis!("Sample Size")
plot!(
        sample_sizes,
        fig3_adjusted_deltas,
        label="Bias Corrected Deltas",
        markershape=:circle,
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        legend=:topright,
)
savefig("./18337FinalProject/plots/delta_ishigami.png")

# Trend is correct, but the scale doesn't match that in the paper.

##
# Try to recreate estimated sensitives for the ishigami function

analytical_values = [0.208, 0.391, 0.156, 0.060] # from Plischke et al. (2013)
sample_size = 4096
X, Y = create_ishigami_samples(sample_size)
delta, delta_adjusted, delta_low, delta_high = delta_moment(X, Y, Ygrid_length=110, num_classes=10)

(delta_low .< analytical_values).== (analytical_values .< delta_high)

#################
## Benchmarking
#################
# Multi-Threading By Factors Time Advantage

sample_sizes = 2 .^ [i for i in 1:14]

non_multi_threaded_times = zeros(length(sample_sizes))
non_multi_threaded_allocs = zeros(length(sample_sizes))
non_multi_threaded_memory = zeros(length(sample_sizes))
multi_threaded_times = zeros(length(sample_sizes))
multi_threaded_allocs = zeros(length(sample_sizes))
multi_threaded_memory = zeros(length(sample_sizes))

for (index, sample_size) in enumerate(sample_sizes)
        println(index, " out of ", length(sample_sizes))

        X, Y = create_ishigami_samples(sample_size)
        non_multi_threaded_benchmark = @benchmark _delta_moment_non_multi($X, $Y)
        multi_threaded_benchmark = @benchmark delta_moment($X, $Y)

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
        title="Delta Moment Estimation Time Benchmarks \n on Ishigami Function"
)
yaxis!("Time (ns)")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_times, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_times.png")

plot(
        sample_sizes,
        non_multi_threaded_memory,
        legend=:topleft,
        markershape=:circle,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="Delta Moment Estimation Memory Benchmarks \n on Ishigami Function"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_memory, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_memory.png")


plot(
        sample_sizes,
        non_multi_threaded_allocs,
        legend=:topleft,
        markershape=:circle,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="Delta Moment Estimation Allocation Benchmarks \n on Ishigami Function"
)
yaxis!("Number of Allocations")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_allocs, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_allocations.png")


## Flame graph

@profile for i in 1:10 delta_moment(X, Y) end
Juno.profiler()
Profile.clear()
