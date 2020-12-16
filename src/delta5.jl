using KernelDensity, StatsBase, Trapz, Random, Distributions, Plots, BenchmarkTools
# https://www.sciencedirect.com/science/article/pii/S0377221712008995
# https://github.com/SALib/SALib/blob/master/src/SALib/analyze/delta.py

function calc_delta(Xi, Y, Ygrid, class_cutoffs)

        # Make sure Ys are not identical, otherwise KDE will be undefined.
        # If Ys are identical, then we know X and Y are independent, so return 0
        if all(Y[1] .== Y)
                return 0
        end

        N = length(Y) # Number of simulations
        @assert length(Xi) == N # Length of Y should equal length of X

         # Model pdf of Y using KDE, kde uses normal kernel by default
        fy = kde(Y, Ygrid).density # eq 23.1

        # Get probability of each y in Ygrid.
        x_rank = competerank(Xi) # Does what scipy.stats rankdata does. If tie, all tied values get same rank.
        d_hat = 0 # the delta estimator.

        # Iterate over each class
        weighted_class_seps = zeros(length(class_cutoffs)-1)
        for j in 1:length(class_cutoffs)-1

                # get X and Y indicies for samples that are in this class (where
                # class designation is based on the X value)
                condition(x) = (x > class_cutoffs[j]) ==  (x <= class_cutoffs[j+1])
                in_class_indices = findall(condition, x_rank)
                number_in_class = length(in_class_indices)

                # Get the the subset of Y values in this class
                in_class_Y = Y[in_class_indices]

                # get the separation between the total y pdf and the condition y pdf induced by this class
                # fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated conditional distribution of y (using normal kernel)
                # curve_diff = abs.(fy.density .- fyc.density) # eq 24

                if ~all(t == in_class_Y[1] for t in in_class_Y)
                        # get the separation between the total y pdf and the condition y pdf induced by this class
                        fyc = kde(in_class_Y, Ygrid).density # eq 23.2 - Estimated condition distribution of y
                        curve_diff = abs.(fy .- fyc) # eq 24
                else
                        curve_diff = fy
                end

                # Use trapezoidal rule to estimate the difference between the curves.
                class_separation = trapz(Ygrid, curve_diff) # eq 25

                # Increment estimator
                weighted_class_seps[j] = number_in_class * class_separation # eq 26
                # d_hat += number_in_class * class_separation # eq 26
        end

        d_hat = sum(weighted_class_seps)
        d_hat = d_hat/(2*N)
        return d_hat
end



function calc_delta2(fy, fyc, Xi, Y, Ygrid, class_cutoffs)

        # Make sure Ys are not identical, otherwise KDE will be undefined.
        # If Ys are identical, then we know X and Y are independent, so return 0
        if all(Y[1] .== Y)
                return 0
        end

        N = length(Y) # Number of simulations
        @assert length(Xi) == N # Length of Y should equal length of X

         # Model pdf of Y using KDE, kde uses normal kernel by default
        fy .= kde(Y, Ygrid).density # eq 23.1

        # Get probability of each y in Ygrid.
        x_rank = competerank(Xi) # Does what scipy.stats rankdata does. If tie, all tied values get same rank.
        d_hat = 0 # the delta estimator.

        # Iterate over each class
        weighted_class_seps = zeros(length(class_cutoffs)-1)
        for j in 1:length(class_cutoffs)-1

                # get X and Y indicies for samples that are in this class (where
                # class designation is based on the X value)
                condition(x) = (x > class_cutoffs[j]) ==  (x <= class_cutoffs[j+1])
                in_class_indices = findall(condition, x_rank)
                number_in_class = length(in_class_indices)

                # Get the the subset of Y values in this class
                in_class_Y = Y[in_class_indices]

                # get the separation between the total y pdf and the condition y pdf induced by this class
                # fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated conditional distribution of y (using normal kernel)
                # curve_diff = abs.(fy.density .- fyc.density) # eq 24

                if ~all(t == in_class_Y[1] for t in in_class_Y)
                        # get the separation between the total y pdf and the condition y pdf induced by this class
                        fyc .= kde(in_class_Y, Ygrid).density # eq 23.2 - Estimated condition distribution of y
                        curve_diff = abs.(fy .- fyc) # eq 24
                else
                        curve_diff = fy
                end

                # Use trapezoidal rule to estimate the difference between the curves.
                class_separation = trapz(Ygrid, curve_diff) # eq 25

                # Increment estimator
                weighted_class_seps[j] = number_in_class * class_separation # eq 26
                # d_hat += number_in_class * class_separation # eq 26
        end

        d_hat = sum(weighted_class_seps)
        d_hat = d_hat/(2*N)
        return d_hat
end


function calc_delta_multithread(Xi, Y, Ygrid, class_cutoffs)

        # Make sure Ys are not identical, otherwise KDE will be undefined.
        # If Ys are identical, then we know X and Y are independent, so return 0
        if all(Y[1] .== Y)
                println("Warning - Ys are identical. Returning 0 for delta estimate.")
                return 0
        end

        N = length(Y) # Number of simulations
        @assert length(Xi) == N # Length of Y should equal length of X

         # Model pdf of Y using KDE, kde uses normal kernel by default
        fy = kde(Y, Ygrid).density # eq 23.1

        # Get probability of each y in Ygrid.
        x_rank = competerank(Xi) # Does what scipy.stats rankdata does. If tie, all tied values get same rank.
        d_hat = 0 # the delta estimator.

        # Iterate over each class
        weighted_class_seps = zeros(length(class_cutoffs)-1)
        Threads.@threads for j in 1:length(class_cutoffs)-1

                # get X and Y indicies for samples that are in this class (where
                # class designation is based on the X value)
                condition(x) = (x > class_cutoffs[j]) ==  (x <= class_cutoffs[j+1])
                in_class_indices = findall(condition, x_rank)
                number_in_class = length(in_class_indices)

                # Get the the subset of Y values in this class
                in_class_Y = Y[in_class_indices]

                # get the separation between the total y pdf and the condition y pdf induced by this class
                # fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated conditional distribution of y (using normal kernel)
                # curve_diff = abs.(fy.density .- fyc.density) # eq 24

                if ~all(t == in_class_Y[1] for t in in_class_Y)
                        # get the separation between the total y pdf and the condition y pdf induced by this class
                        fyc = kde(in_class_Y, Ygrid).density # eq 23.2 - Estimated condition distribution of y
                        curve_diff = abs.(fy .- fyc) # eq 24
                else
                        curve_diff = fy
                end

                # Use trapezoidal rule to estimate the difference between the curves.
                class_separation = trapz(Ygrid, curve_diff) # eq 25

                # Increment estimator
                weighted_class_seps[j] = number_in_class * class_separation # eq 26
                # d_hat += number_in_class * class_separation # eq 26
        end

        d_hat = sum(weighted_class_seps)
        d_hat = d_hat/(2*N)
        return d_hat
end


function bias_reduced_delta(Y, Ygrid, X, m, num_resamples, conf_level)
        """Plischke et al. 2013 bias reduction technique (eqn 30)"""
        d = zeros(num_resamples)
        d_hat = calc_delta(X, Y)
end


function delta_moment_analyze(X_matrix, Y; num_resamples=500, conf_level=0.95, seed=nothing, Ygrid_length=2048, num_classes=nothing)

        if seed != nothing
                Random.seed!(seed)
        end

        N = length(Y)

        # Create number of classes and class cutoffs.
        if num_classes == nothing
                exp = (2 / (7 + tanh((1500 - N) / 500)))
                M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
        else
                M = num_classes
        end
        class_cutoffs =  range(0, N, length=M+1) # class cutoffs.

        # quadrature points.
        # Length should be a power of 2 for efficient FFT in kernel density estimates.
        Ygrid = range(minimum(Y), maximum(Y), length=Ygrid_length)

        dims = size(X_matrix)
        if length(dims) == 2
                num_factors = dims[2]
        else
                num_factors = 1
        end

        deltas = zeros(num_factors)
        adjusted_deltas = zeros(num_factors)
        adjusted_deltas_conf = zeros(num_factors)
        adjusted_deltas_low = zeros(num_factors)
        adjusted_deltas_high = zeros(num_factors)
        for factor_num in 1:num_factors
                Xi = view(X_matrix, :, factor_num)

                delta = calc_delta(Xi, Y, Ygrid, class_cutoffs)
                deltas[factor_num] = delta


                # eq. 30, bias reduction via bootstrapping.
                d = zeros(num_resamples)
                r = rand(1:N, num_resamples, N)
                for i in 1:num_resamples
                        r_i = r[i, :]
                        d[i] = calc_delta(Xi[r_i], Y[r_i], Ygrid, class_cutoffs)
                end

                d = 2 * delta .- d

                adjusted_deltas[factor_num] = mean(d)
                diff = quantile(Normal(0.0, 1.0), 0.5+conf_level/2)*std(d)/(sqrt(num_resamples))
                adjusted_deltas_low[factor_num] = adjusted_deltas[factor_num] - diff
                adjusted_deltas_high[factor_num] = adjusted_deltas[factor_num] + diff
        end

        return deltas, adjusted_deltas, adjusted_deltas_low, adjusted_deltas_high
end



function delta_moment_analyze2(X_matrix, Y; num_resamples=500, conf_level=0.95, seed=nothing, Ygrid_length=2048, num_classes=nothing)

        if seed != nothing
                Random.seed!(seed)
        end

        N = length(Y)

        # Create number of classes and class cutoffs.
        if num_classes == nothing
                exp = (2 / (7 + tanh((1500 - N) / 500)))
                M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
        else
                M = num_classes
        end
        class_cutoffs =  range(0, N, length=M+1) # class cutoffs.

        # quadrature points.
        # Length should be a power of 2 for efficient FFT in kernel density estimates.
        Ygrid = range(minimum(Y), maximum(Y), length=Ygrid_length)

        dims = size(X_matrix)
        if length(dims) == 2
                num_factors = dims[2]
        else
                num_factors = 1
        end

        deltas = zeros(num_factors)
        adjusted_deltas = zeros(num_factors)
        adjusted_deltas_conf = zeros(num_factors)
        adjusted_deltas_low = zeros(num_factors)
        adjusted_deltas_high = zeros(num_factors)

        fy = zeros(Ygrid_length)
        fyc = zeros(Ygrid_length)
        for factor_num in 1:num_factors
                Xi = view(X_matrix, :, factor_num)

                delta = calc_delta2(fy, fyc, Xi, Y, Ygrid, class_cutoffs)
                deltas[factor_num] = delta


                # eq. 30, bias reduction via bootstrapping.
                d = zeros(num_resamples)
                r = rand(1:N, num_resamples, N)
                for i in 1:num_resamples
                        r_i = r[i, :]
                        # d[i] = calc_delta2(fy, fyc, view(Xi, r_i), view(Y, r_i), Ygrid, class_cutoffs)
                        d[i] = calc_delta2(fy, fyc, Xi[r_i], Y[r_i], Ygrid, class_cutoffs)

                end

                d = 2 * delta .- d

                adjusted_deltas[factor_num] = mean(d)
                diff = quantile(Normal(0.0, 1.0), 0.5+conf_level/2)*std(d)/(sqrt(num_resamples))
                adjusted_deltas_low[factor_num] = adjusted_deltas[factor_num] - diff
                adjusted_deltas_high[factor_num] = adjusted_deltas[factor_num] + diff
        end

        return deltas, adjusted_deltas, adjusted_deltas_low, adjusted_deltas_high
end



function delta_moment_analyze_mutlithreaded(X_matrix, Y; num_resamples=500, conf_level=0.95, seed=nothing, Ygrid_length=2048, num_classes=nothing)

        if seed != nothing
                Random.seed!(seed)
        end

        N = length(Y)

        # Create number of classes and class cutoffs.
        if num_classes == nothing
                exp = (2 / (7 + tanh((1500 - N) / 500)))
                M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
        else
                M = num_classes
        end
        class_cutoffs =  range(0, N, length=M+1) # class cutoffs.

        # quadrature points.
        # Length should be a power of 2 for efficient FFT in kernel density estimates.
        Ygrid = range(minimum(Y), maximum(Y), length=Ygrid_length)

        dims = size(X_matrix)
        if length(dims) == 2
                num_factors = dims[2]
        else
                num_factors = 1
        end

        deltas = zeros(num_factors)
        adjusted_deltas = zeros(num_factors)
        adjusted_deltas_conf = zeros(num_factors)
        adjusted_deltas_low = zeros(num_factors)
        adjusted_deltas_high = zeros(num_factors)
        Threads.@threads for factor_num in 1:num_factors
                Xi = view(X_matrix, :, factor_num)

                delta = calc_delta(Xi, Y, Ygrid, class_cutoffs)
                deltas[factor_num] = delta


                # eq. 30, bias reduction via bootstrapping.
                d = zeros(num_resamples)
                r = rand(1:N, num_resamples, N)
                for i in 1:num_resamples
                        r_i = r[i, :]
                        d[i] = calc_delta(Xi[r_i], Y[r_i], Ygrid, class_cutoffs)
                end

                d = 2 * delta .- d

                adjusted_deltas[factor_num] = mean(d)
                diff = quantile(Normal(0.0, 1.0), 0.5+conf_level/2)*std(d)/(sqrt(num_resamples))
                adjusted_deltas_low[factor_num] = adjusted_deltas[factor_num] - diff
                adjusted_deltas_high[factor_num] = adjusted_deltas[factor_num] + diff
        end

        return deltas, adjusted_deltas, adjusted_deltas_low, adjusted_deltas_high
end



function ishigami(X)
        """
                Ishigami function. Takes a vector, returns a scalar.
                Saltelli et al., 2004
        """
        X1, X2, X3, X4 = X
        Y = sin(X1) + 7*sin(X2)^2 + 0.1*X3^4*sin(X1)
        return Y
end


function create_ishigami_samples(sample_size)

        X = rand(Uniform(-pi, pi), sample_size, 4)
        Y = zeros(sample_size)
        for sample_num in 1:sample_size
                ishigami_input = X[sample_num, :]
                Y[sample_num] = ishigami(ishigami_input)
        end
        return X, Y
end

##
# Multi-Threading By Factors Time Advantage

sample_sizes = 2 .^ [i for i in 1:14]

non_multi_threaded_times = zeros(length(sample_sizes))
non_multi_threaded_allocs = zeros(length(sample_sizes))
non_multi_threaded_memory = zeros(length(sample_sizes))
multi_threaded_times = zeros(length(sample_sizes))
multi_threaded_allocs = zeros(length(sample_sizes))
multi_threaded_memory = zeros(length(sample_sizes))

for (index, sample_size) in enumerate(sample_sizes)
        println(index)
        X, Y = create_ishigami_samples(sample_size)
        non_multi_threaded_benchmark = @benchmark delta_moment_analyze($X, $Y)
        multi_threaded_benchmark = @benchmark delta_moment_analyze_mutlithreaded($X, $Y)

        non_multi_threaded_times[index] = time(non_multi_threaded_benchmark)
        non_multi_threaded_allocs[index] = allocs(non_multi_threaded_benchmark)
        non_multi_threaded_memory[index] = memory(non_multi_threaded_benchmark)
        multi_threaded_times[index] = time(multi_threaded_benchmark)
        multi_threaded_allocs[index] = allocs(multi_threaded_benchmark)
        multi_threaded_memory[index] = memory(multi_threaded_benchmark)
end


scatter(
        sample_sizes,
        non_multi_threaded_times,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="Delta Moment Estimation Time Benchmarks \n on Ishigami Function"
)
yaxis!("Time (ns)")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_times, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_times.png")

scatter(
        sample_sizes,
        non_multi_threaded_memory,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="Delta Moment Estimation Memory Benchmarks \n on Ishigami Function"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_memory, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_memory.png")


scatter(
        sample_sizes,
        non_multi_threaded_allocs,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="Delta Moment Estimation Allocation Benchmarks \n on Ishigami Function"
)
yaxis!("Number of Allocations")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_allocs, label="Multithreaded")
savefig("./18337FinalProject/plots/delta_allocations.png")



##

###################################
# Test independence
X = hcat([i for i in 1:9], [1 for i in 1:9])
Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]
deltas, adjusted_deltas, adjusted_deltas_low, adjusted_deltas_high = delta_moment_analyze(X, Y)



###################################
# Try to recreate figure 3. # https://www.sciencedirect.com/science/article/pii/S0377221712008995#e0040
sample_sizes = 2 .^ [i for i in 9:14]
fig3_deltas = zeros(length(sample_sizes))
fig3_adjusted_deltas = zeros(length(sample_sizes))
fig3_adj_delta_95low = zeros(length(sample_sizes))
fig3_adj_delta_95high = zeros(length(sample_sizes))
for (index, sample_size) in enumerate(sample_sizes)
        X, Y = create_ishigami_samples(sample_size)
        res = delta_moment_analyze(X[:,4], Y, Ygrid_length=110, num_classes=10)
        fig3_deltas[index], fig3_adjusted_deltas[index], fig3_adj_delta_95low[index], fig3_adj_delta_95high[index] = res[1][1], res[2][1], res[3][1], res[4][1]
end


scatter(
        sample_sizes,
        fig3_deltas,
        label="Deltas",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        legend=:topright,
        title="Ishigami Function Inactive Variable Sensitivity Measures"
)
yaxis!("Delta Sensitivity Measure")
xaxis!("Sample Size")
scatter!(
        sample_sizes,
        fig3_adjusted_deltas,
        label="Bias Corrected Deltas",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        legend=:topright,
)


# Trend is correct, but the scale is off. Possible due to

###################################

# Try to recreate estimated sensitives for the ishigami function [0.208, 0.391, 0.156, 0.060] from the paper
sample_size = 512
X, Y = create_ishigami_samples(sample_size)
delta, delta_adjusted, delta_low, delta_high = delta_moment_analyze(X, Y, Ygrid_length=110, num_classes=10)


(delta_low .<  [0.208, 0.391, 0.156, 0.060]).== ([0.208, 0.391, 0.156, 0.060] .< delta_high)

############################
using Profile, BenchmarkTools


# sample_sizes = 2 .^ [i for i in 9:14]
sample_sizes = 2 .^ [i for i in 1:9]
benchmark_time = zeros(length(sample_sizes))
benchmark_memory = zeros(length(sample_sizes))
benchmark_allocs = zeros(length(sample_sizes))
for (index, sample_size) in enumerate(sample_sizes)
        X, Y = create_ishigami_samples(sample_size)
        benchmark_res = @benchmark delta_moment_analyze(X, Y, Ygrid_length=110, num_classes=10)
        benchmark_time[index] = time(benchmark_res)
        benchmark_memory[index] = memory(benchmark_res)
        benchmark_allocs[index] = allocs(benchmark_res)
end

# Profiling



@profile for i in 1:10 delta_moment_analyze(X, Y, Ygrid_length=110, num_classes=10) end
Juno.profiler()
Profile.clear()


benchmark_res = @benchmark delta_moment_analyze(X, Y, Ygrid_length=110, num_classes=100)


##

# Factor of 2 comparison
sample_size = 2^14
X, Y = create_ishigami_samples(sample_size)

@time delta_moment_analyze(X, Y, Ygrid_length=110, num_classes=10, num_resamples=10)

@time delta_moment_analyze2(X, Y, Ygrid_length=110, num_classes=10, num_resamples=10)



using Profile
@profile for i in 1:10 delta_moment_analyze2(X, Y, Ygrid_length=110, num_classes=10, num_resamples=10) end
Juno.profiler()
Profile.clear()

print("HI")
