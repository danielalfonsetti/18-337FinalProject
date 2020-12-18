
using BenchmarkTools, Profile, Gadfly, Cairo, Fontconfig, Distributions

include("../src/EASI.jl")
include("../src/util.jl")

#############################
## Testing for correctness ##
#############################

# Recreate figure 4 in https://www.sciencedirect.com/science/article/pii/S0951832009002579

ps1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(ps1)

# get analytical values
vis = [1/(3*(1+a)^2) for a in ps1] # conditional variances
v = prod([vi+1 for vi in vis]) - 1 # total variances.
analytic_values = vis ./ v # variance fractions


runs = 150
sample_sizes = [100, 300, 1000, 3000, 10000, 30000]

s_runs1 = zeros(runs, length(sample_sizes))
s_corrected_runs1 = zeros(runs, length(sample_sizes))
s_runs4 = zeros(runs, length(sample_sizes))
s_corrected_runs4 = zeros(runs, length(sample_sizes))
for (sample_size_index, sample_size) in enumerate(sample_sizes)
    for run in 1:runs

        X = rand(Uniform(0,1), sample_size, length(ps1))

        Y = zeros(sample_size)
        for i in 1:sample_size
            Y[i] = sobol_function1(@view X[i,:])
        end

        s, s_corrected = EASI(X, Y, 6)

        s_runs1[run, sample_size_index] = s[1]
        s_corrected_runs1[run, sample_size_index] = s_corrected[1]

        s_runs4[run, sample_size_index] = s[4]
        s_corrected_runs4[run, sample_size_index] = s_corrected[4]
    end
end


Gadfly.push_theme(Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt))

x = repeat(string.(sample_sizes), inner=150)
y = reshape(s_corrected_runs1, 150*6)

p = Gadfly.plot(
    x=x,
    y=y,
    yintercept=[analytic_values[1]],
    Geom.hline(style=:solid, color=["red"]),
    Geom.boxplot,
    Guide.title("Sobol S1 Estimate (Corrected) for EASI"),
    Guide.xlabel("Sample Size"),
    Guide.ylabel("Fraction of Variance")
)
Gadfly.draw(PNG("./18337FinalProject/plots/EASI_sobol_s1_corrected.png"), p)


y = reshape(s_corrected_runs4, 150*6)
p = Gadfly.plot(
    x=x,
    y=y,
    yintercept=[analytic_values[4]],
    Geom.hline(style=:solid, color=["red"]),
    Geom.boxplot,
    Guide.title("Sobol S4 Estimate (Corrected) for EASI"),
    Guide.xlabel("Sample Size"),
    Guide.ylabel("Fraction of Variance")
)
Gadfly.draw(PNG("./18337FinalProject/plots/EASI_sobol_s4_corrected.png"), p)

y = reshape(s_runs1, 150*6)
p = Gadfly.plot(
    x=x,
    y=y,
    yintercept=[analytic_values[1]],
    Geom.hline(style=:solid, color=["red"]),
    Geom.boxplot,
    Guide.title("Sobol S1 Estimate for EASI"),
    Guide.xlabel("Sample Size"),
    Guide.ylabel("Fraction of Variance")
)
Gadfly.draw(PNG("./18337FinalProject/plots/EASI_sobol_s1.png"), p)

y = reshape(s_runs4, 150*6)
p = Gadfly.plot(
    x=x,
    y=y,
    yintercept=[analytic_values[4]],
    Geom.hline(style=:solid, color=["red"]),
    Geom.boxplot,
    Guide.title("Sobol S4 Estimate for EASI"),
    Guide.xlabel("Sample Size"),
    Guide.ylabel("Fraction of Variance")
)
Gadfly.draw(PNG("./18337FinalProject/plots/EASI_sobol_s4.png"), p)

##################
## Benchmarking ##
##################

# Compare parallelized/multi threaded EASI to non-parallelized EASI

sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(sobol_function_params1)

sample_sizes = 2 .^ [i for i in 9:14]

non_multi_threaded_times = zeros(length(sample_sizes))
non_multi_threaded_allocs = zeros(length(sample_sizes))
non_multi_threaded_memory = zeros(length(sample_sizes))
multi_threaded_times = zeros(length(sample_sizes))
multi_threaded_allocs = zeros(length(sample_sizes))
multi_threaded_memory = zeros(length(sample_sizes))

for (index, sample_size) in enumerate(sample_sizes)
        println(index, " out of ", length(sample_sizes))

        X = rand(Uniform(0,1), sample_size, length(ps1))
        Y = zeros(sample_size)

        non_multi_threaded_benchmark = @benchmark _EASI_non_multi($X, $Y)
        multi_threaded_benchmark = @benchmark EASI($X, $Y)

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
        xaxis=:log,
        markershape=:circle,
        label="Non-multithreaded",
        xticks = (sample_sizes, sample_sizes),
        title="EASI Time Benchmarks\non Sobol Function"
)
yaxis!("Time (ns)")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_times, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/EASI_times.png")

plot(
        sample_sizes,
        non_multi_threaded_memory,
        legend=:topleft,
        markershape=:circle,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="EASI Memory Benchmarks\non Sobol Function"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_memory, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/EASI_memory.png")

plot(
        sample_sizes,
        non_multi_threaded_allocs,
        legend=:topleft,
        markershape=:circle,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="EASI Allocation Benchmarks\non Sobol Function"
)
yaxis!("Number of Allocations")
xaxis!("Number of Samples")
plot!(sample_sizes, multi_threaded_allocs, markershape=:circle, label="Multithreaded")
savefig("./18337FinalProject/plots/EASI_allocations.png")


## Flame graph
@profile for i in 1:1000 EASI(X, Y) end
Juno.profiler()
Profile.clear()
