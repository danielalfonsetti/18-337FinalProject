# https://www.sciencedirect.com/science/article/pii/S0951832005001444
# Page 168 in Global Sensitivty Analysis - The Primer


# https://dsp.stackexchange.com/questions/15680/how-to-interpret-results-of-fourier-transform-by-using-python

# https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/fast.html#analyze

using FFTW, Random, Statistics, StatsBase


allsame(x) = all(y -> y == first(x), x)

function conditional_print(x, verbose)
    if verbose
        println!(x)
    end
end

function RBD_FAST_multi(model, N, group_encoding, frequencies, M=6)
    """
    N = Number of samples
    K = Number of factors
    G = Number of groups
        RBD is applied independently within each group,
        meaning each group gets a different frequency.

    frequencies: should be a separate frequency for each group
    Group 1 frequency should be specified at frequencies[1], etc.
    Group 2 frequency should be specified at frequencies[2], etc.
    """

    # println("Running RBD_Fast_Multi")
    @assert N >= 2*M*maximum(frequencies) # Nyquist–Shannon

    K = length(group_encoding)
    group_sizes = countmap(group_encoding)
    @assert allsame(collect(values(group_sizes)))
    group_size = group_sizes[1]


    perms_matrix = zeros(Int64, group_size, N)
    for i in 1:group_size
        perms_matrix[i,:] = randperm(N)
    end


    group_id_to_count_so_far = Dict(i => 0 for i in keys(group_sizes))


    s0_matrix = zeros(N, K)
    for i in 1:K
        s0_matrix[:,i] = collect(-pi:2*pi/(N-1):pi)
    end

    # Compute inputs
    s_matrix = zeros(N, K)
    x_matrix = zeros(N, K)

    for i in 1:K
        group_id = group_encoding[i]
        ω = frequencies[group_id]


        group_id_to_count_so_far[group_id] += 1
        perm = perms_matrix[group_id_to_count_so_far[group_id],:]

        # println("----")
        # println(i)
        # println(perm)
        # println(ω)

        s_matrix[:,i] = s0_matrix[:,i][perm]
        x_matrix[:,i] = 0.5.+asin.(sin.(ω .* s_matrix[:,i]))./pi  # Assumes x values are uniformly distributed
    end

    # println("s_matrix: ", s_matrix)
    # println("x_matrix: ", x_matrix)

    # Compute outputs
    y_out = zeros(N)
    for i in 1:N
        y_out[i] = model(x_matrix[i,:])
    end
    # println("Y out: ", y_out)

    # Iterate over factors
    sensitivites = zeros(K)
    # println("Iterating over factors")
    for i in 1:K

        group_id = group_encoding[i]
        ω = frequencies[group_id]


        s = s_matrix[:,i]
        ranks = sortperm(s)

        # Order Ys by how they would occur if they were monotonically
        # increasing as s increased.
        y_reordered = y_out[ranks]

        spectrum=((abs.(fft(y_reordered))).^2)./N; # get the power density spectrum
        V=sum(spectrum[2:N])
        harmonics = [ω*i for i=1:M] # get the harmonics we want to consider for the partial variance.
        V1 = 2*sum(spectrum[harmonics.+1])
        Si = V1/V

        sensitivites[i]  = Si
        # unskew
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

num_samples = 1000
num_harmonics = 3
input_groups = [1]
frequencies = [1]
res = RBD_FAST_multi(random_function, num_samples, input_groups, frequencies, num_harmonics)
res[1]


# ps1 = [0, 0.5, 3, 9, 99, 99]
# sobol_function1 = create_function_of_sobol(ps1)
# res = RBD_FAST_multi(sobol_function1, 9, [1, 1, 1, 2, 2, 2], [1, 2], 2)
# res



ps1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(ps1)
num_experiments = 50
experiments = zeros(num_experiments, 8)
for i in 1:num_experiments
    experiments[i,:] = RBD_FAST_multi(sobol_function1, 10000, [1, 1, 1, 1, 2, 2, 2, 2], [11, 35])
end

standard_devs = mapslices(Statistics.std, experiments, dims=1)
means = mapslices(Statistics.mean, experiments, dims=1)
println("Means: ", means)



ps1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(ps1)
num_experiments = 50
experiments = zeros(num_experiments, 8)
for i in 1:num_experiments
    experiments[i,:] = RBD_FAST_multi(sobol_function1, 10000, [1, 1, 1, 1, 1, 1, 1, 1], [1])
end

standard_devs = mapslices(Statistics.std, experiments, dims=1)
means = mapslices(Statistics.mean, experiments, dims=1)
println("Means: ", means)


# ------


ps1 = cat([0 for i in 1:4], [1 for i in 1:4], [9, 9], [99 for i in 1:90], dims=1)
sobol_function1 = create_function_of_sobol(ps1)
num_experiments = 50
experiments = zeros(num_experiments, 8)
for i in 1:num_experiments
    experiments[i,:] = RBD_FAST_multi(sobol_function1, 10000, [1, 1, 1, 1, 1, 1, 1, 1], [1])
end

standard_devs = mapslices(Statistics.std, experiments, dims=1)
means = mapslices(Statistics.mean, experiments, dims=1)
println("Means: ", means)
