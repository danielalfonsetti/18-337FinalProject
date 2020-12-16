# https://www.sciencedirect.com/science/article/pii/S0951832005001444
# Page 168 in Global Sensitivty Analysis - The Primer

# https://dsp.stackexchange.com/questions/15680/how-to-interpret-results-of-fourier-transform-by-using-python
# https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/fast.html#analyze

using FFTW, Random, Statistics, StatsBase
allsame(x) = all(y -> y == first(x), x)


function RBD_FAST(model, N, group_encoding, frequencies, M=6)
    """
    param: model
        A scalar valued function. Sensitivity indices are computed
        between each of this function's K inputs and its scalar output.

    param: N
        Number of samples

    param: group_encoding
        A list of numbers representing which factor is in which group.
        The group number for the i-th factor is encoded at group group_encoding[i]
        The length of this list (i.e. the number of factors) should be equal to
        the number of inputs that the "model" takes.

        There should be the same number of factors within each group.
        An error will be thrown otherwise.

    param: frequencies

        Group 1 frequency should be specified at frequencies[1],
        Group 2 frequency should be specified at frequencies[2], etc.

        The number of frequencies provided should equal the number of groups.
        An error will be thrown otherwise.

    param: N
        Number of harmonics to consider during power spectral density analysis.
    """
    @assert N >= 2*M*maximum(frequencies) # Nyquist–Shannon inequality

    K = length(group_encoding)
    group_to_group_count = countmap(group_encoding)

    @assert allsame(collect(values(group_to_group_count))) "Number of factors in each group must be equal."
    group_size = group_to_group_count[1]

    groups = keys(group_to_group_count)
    @assert length(groups) == length(frequencies) "The number of groups must equal the number of frequencies"

    for group in groups
        @assert group <= length(frequencies) "Group numbers must correspond to an index provided in 'frequencies'"
    end

    # Give a unique permutation to each group.
    perms_matrix = zeros(Int64, group_size, N)
    for i in 1:group_size
        perms_matrix[i,:] = randperm(N)
    end

    # We need to store a dictionary that counts
    # how many times we have come across any factor within a particular group.
    # This allows us to provide a unique permutation to each factor within
    # the same group.
    group_id_to_count_so_far = Dict(i => 0 for i in keys(group_to_group_count))

    # Initalize matrix containing range of values of the parametric variable
    # along each column (factor).
    s0_matrix = zeros(N, K)
    for i in 1:K
        s0_matrix[:,i] = collect(-pi:2*pi/(N-1):pi)
    end

    # Compute inputs
    s_matrix = zeros(N, K) # number of samples X number of factors
    x_matrix = zeros(N, K)

    for i in 1:K

        # Get the frequency associated with the group that the current factor is in.
        group_id = group_encoding[i]
        ω = frequencies[group_id]

        # Within each group, each factor in that group gets a different permutation
        # (even though all factors within that group get the same frequency)
        group_id_to_count_so_far[group_id] += 1
        perm = @view perms_matrix[group_id_to_count_so_far[group_id],:]

        # store permutation of parametric variable
        s_matrix[:,i] = @view s0_matrix[:,i][perm]

        # Generate corresponding x values with
        # imposed frequency ω based on permutation of the parametric variable.
        x_matrix[:,i] = 0.5.+asin.(sin.(ω .* @view s_matrix[:,i]))./pi
    end

    # Compute outputs
    y_out = zeros(N)
    for i in 1:N
        y_out[i] = model(@view x_matrix[i,:])
    end

    # Iterate over factors
    spectrum = zeros(N)  # preallocate spectrum
    sensitivites = zeros(K)
    for i in 1:K

        group_id = group_encoding[i]
        ω = frequencies[group_id]

        s = @view s_matrix[:,i]
        ranks = sortperm(s)

        # Order Ys by how they would occur if they were
        # monotonically increasing as the
        # parametric variable s (not its permutation) increased.
        y_reordered = @view y_out[ranks]

        # get the power density spectrum
        spectrum .= ((abs.(fft(y_reordered))).^2)./N;
        V=sum(@view spectrum[2:N])
        # get the harmonics we want to consider for the partial variance.
        harmonics = [ω*i for i=1:M]
        V1 = 2*sum(@view spectrum[harmonics.+1])
        Si = V1/V

        # unskew the sensitivies
        lambda = 2*M/N
        sensitivites[i] = Si - (lambda / (1 - lambda)) * (1-Si)
    end

    return sensitivites
end


function RBD_FAST_multi_threading(model, N, group_encoding, frequencies, M=6)
    """
    param: model
        A scalar valued function. Sensitivity indices are computed
        between each of this function's K inputs and its scalar output.

    param: N
        Number of samples

    param: group_encoding
        A list of numbers representing which factor is in which group.
        The group number for the i-th factor is encoded at group group_encoding[i]
        The length of this list (i.e. the number of factors) should be equal to
        the number of inputs that the "model" takes.

        There should be the same number of factors within each group.
        An error will be thrown otherwise.

    param: frequencies

        Group 1 frequency should be specified at frequencies[1],
        Group 2 frequency should be specified at frequencies[2], etc.

        The number of frequencies provided should equal the number of groups.
        An error will be thrown otherwise.

    param: N
        Number of harmonics to consider during power spectral density analysis.
    """
    @assert N >= 2*M*maximum(frequencies) # Nyquist–Shannon inequality

    K = length(group_encoding)
    group_to_group_count = countmap(group_encoding)

    @assert allsame(collect(values(group_to_group_count))) "Number of factors in each group must be equal."
    group_size = group_to_group_count[1]

    groups = keys(group_to_group_count)
    @assert length(groups) == length(frequencies) "The number of groups must equal the number of frequencies"

    for group in groups
        @assert group <= length(frequencies) "Group numbers must correspond to an index provided in 'frequencies'"
    end

    # Give a unique permutation to each group.
    perms_matrix = zeros(Int64, group_size, N)
    for i in 1:group_size
        perms_matrix[i,:] = randperm(N)
    end

    # We need to store a dictionary that counts
    # how many times we have come across any factor within a particular group.
    # This allows us to provide a unique permutation to each factor within
    # the same group.
    group_id_to_count_so_far = Dict(i => 0 for i in keys(group_to_group_count))

    # Initalize matrix containing range of values of the parametric variable
    # along each column (factor).
    s0_matrix = zeros(N, K)
    for i in 1:K
        s0_matrix[:,i] = collect(-pi:2*pi/(N-1):pi)
    end

    # Compute inputs
    s_matrix = zeros(N, K) # number of samples X number of factors
    x_matrix = zeros(N, K)

    Threads.@threads for i in 1:K

        # Get the frequency associated with the group that the current factor is in.
        group_id = group_encoding[i]
        ω = frequencies[group_id]

        # Within each group, each factor in that group gets a different permutation
        # (even though all factors within that group get the same frequency)
        group_id_to_count_so_far[group_id] += 1
        perm = @view perms_matrix[group_id_to_count_so_far[group_id],:]

        # store permutation of parametric variable
        s_matrix[:,i] = @view s0_matrix[:,i][perm]

        # Generate corresponding x values with
        # imposed frequency ω based on permutation of the parametric variable.
        x_matrix[:,i] = 0.5.+asin.(sin.(ω .* @view s_matrix[:,i]))./pi
    end

    # Compute outputs
    y_out = zeros(N)
    Threads.@threads for i in 1:N
        y_out[i] = model(@view x_matrix[i,:])
    end

    # Iterate over factors
    spectrum = zeros(N)  # preallocate spectrum
    sensitivites = zeros(K)
    Threads.@threads for i in 1:K

        group_id = group_encoding[i]
        ω = frequencies[group_id]

        s = @view s_matrix[:,i]
        ranks = sortperm(s)

        # Order Ys by how they would occur if they were
        # monotonically increasing as the
        # parametric variable s (not its permutation) increased.
        y_reordered = @view y_out[ranks]

        # get the power density spectrum
        spectrum .= ((abs.(fft(y_reordered))).^2)./N;
        V=sum(@view spectrum[2:N])
        # get the harmonics we want to consider for the partial variance.
        harmonics = [ω*i for i=1:M]
        V1 = 2*sum(@view spectrum[harmonics.+1])
        Si = V1/V

        # unskew the sensitivies
        lambda = 2*M/N
        sensitivites[i] = Si - (lambda / (1 - lambda)) * (1-Si)
    end

    return sensitivites
end

function create_function_of_sobol(ps)
    """
    Create a function of sobol, parameterized by the the vector ps.
    The function of sobol is widely used as a sensitivity analysis benchmark.
    https://www.sciencedirect.com/science/article/pii/S0951832005001444
    """

    g(x, p) = (abs(4*x-2)+p)/(1+p) # for 0 <= Xi <= 1 and ai >= 0
    gs = [k(x)=g(x,p) for p in ps] # create functions

    function f(xs)
        @assert length(gs) == length(xs)
        total = 1
        for (index, (g, x)) in enumerate(zip(gs, xs))
            total *= g(x)
        end
        return total
    end
    return f
end


function random_function(x)
    return rand()
end

##
#### Recreate table 1 in https://www.sciencedirect.com/science/article/pii/S0951832005001444#tbl2

sobol_function_params1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
frequencies = [11, 35]
group_assignments = [1, 1, 1, 1, 2, 2, 2, 2]
sobol_function1 = create_function_of_sobol(sobol_function_params1)
num_experiments = 50
results_samplesize500 = zeros(num_experiments, 8)
results_samplesize1000 = zeros(num_experiments, 8)
results_samplesize2000 = zeros(num_experiments, 8)

for i in 1:num_experiments
    results_samplesize500[i,:] = RBD_FAST_multi(sobol_function1, 500, group_assignments, frequencies)
    results_samplesize1000[i,:] = RBD_FAST_multi(sobol_function1, 1000, group_assignments, frequencies)
    results_samplesize2000[i,:] = RBD_FAST_multi(sobol_function1, 2000, group_assignments, frequencies)
end

standard_devs_ss500 = mapslices(Statistics.std, results_samplesize500, dims=1)
means_ss500 = mapslices(Statistics.mean, results_samplesize500, dims=1)
println("Means (ss=500): ", means_ss500)
println("Standard Dev. (ss=500): ", standard_devs_ss500)

standard_devs_ss1000 = mapslices(Statistics.std, results_samplesize1000, dims=1)
means_ss1000 = mapslices(Statistics.mean, results_samplesize1000, dims=1)
println("Means (ss=1000): ", means_ss1000)
println("Standard Dev. (ss=1000): ", standard_devs_ss1000)

standard_devs_ss2000 = mapslices(Statistics.std, results_samplesize2000, dims=1)
means_ss2000 = mapslices(Statistics.mean, results_samplesize2000, dims=1)
println("Means (ss=2000): ", means_ss2000)
println("Standard Dev. (ss=2000): ", standard_devs_ss2000)

##





##

using Profile
@profile for i in 1:10 RBD_FAST_multi(sobol_function1, 500, group_assignments, frequencies) end
Juno.profiler()
Profile.clear()

@time RBD_FAST_multi(sobol_function1, 2000, group_assignments, frequencies)





using Profile
@profile for i in 1:10 RBD_FAST_multi(sobol_function1, 500, group_assignments, frequencies) end
Juno.profiler()
Profile.clear()


using BenchmarkTools

sample_sizes = 2 .^ [i for i in 9:14]

non_multi_threaded_times = zeros(length(sample_sizes))
non_multi_threaded_allocs = zeros(length(sample_sizes))
non_multi_threaded_memory = zeros(length(sample_sizes))
multi_threaded_times = zeros(length(sample_sizes))
multi_threaded_allocs = zeros(length(sample_sizes))
multi_threaded_memory = zeros(length(sample_sizes))

for (index, sample_size) in enumerate(sample_sizes)
        print(index)
        non_multi_threaded_benchmark = @benchmark RBD_FAST($sobol_function1, $sample_size, $group_assignments, $frequencies)
        multi_threaded_benchmark = @benchmark RBD_FAST_multi_threading($sobol_function1, $sample_size, $group_assignments, $frequencies)

        non_multi_threaded_times[index] = time(non_multi_threaded_benchmark)
        non_multi_threaded_allocs[index] = allocs(non_multi_threaded_benchmark)
        non_multi_threaded_memory[index] = memory(non_multi_threaded_benchmark)
        multi_threaded_times[index] = time(multi_threaded_benchmark)
        multi_threaded_allocs[index] = allocs(multi_threaded_benchmark)
        multi_threaded_memory[index] = memory(multi_threaded_benchmark)
end

using Plots

scatter(
        sample_sizes,
        non_multi_threaded_times,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="RBD-FAST Time Benchmarks \n on Sobol Function"
)
yaxis!("Time (ns)")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_times, label="Multithreaded")
savefig("./18337FinalProject/plots/rbd_fast_times.png")

scatter(
        sample_sizes,
        non_multi_threaded_memory,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="RBD-FAST Memory Benchmarks \n on Sobol Function"
)
yaxis!("Memory (MiB)")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_memory, label="Multithreaded")
savefig("./18337FinalProject/plots/rbd_fast_memory.png")


scatter(
        sample_sizes,
        non_multi_threaded_allocs,
        legend=:topleft,
        label="Non-multithreaded",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        title="RBD-FAST Allocation Benchmarks \n on Sobol Function"
)
yaxis!("Number of Allocations")
xaxis!("Number of Samples")
scatter!(sample_sizes, multi_threaded_allocs, label="Multithreaded")
savefig("./18337FinalProject/plots/rbd_fast_allocations.png")



# N = 10
# s0 = collect(-pi:2*pi/(N-1):pi)
# perm = randperm(N)
# s = s0[perm]
#
#
# ranks = sortperm(s)
# res = s[ranks]
