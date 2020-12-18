"""
Code based on the theory presented in:
    Saltelli, A. (2008). Global sensitivity analysis: The primer. Chichester: Wiley, pp. 167-169.
and
    S. Tarantola, D. Gatelli and T. Mara (2006)
    "Random Balance Designs for the Estimation of First Order Global Sensitivity Indices",
    Reliability Engineering and System Safety, 91:6, 717-727
"""

using FFTW, Random, Statistics, StatsBase, Distributions
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

    return:
        A vector of length equal to the number of factors, where the element at the ith index is
        the estimated fractional contribution of the ith factor to the
        overall output variance.
    """

    @assert maximum(frequencies) <= (N-1)/(2*M)

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
        # perms_matrix[i,:] = 1:N
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
        # s0_matrix[:,i] = collect(-pi:2*pi/(N-1):pi)
        s0_matrix[:,i] = range(-pi, stop=pi, length=N)
        # s0_matrix[:,i] = rand(Uniform(-pi, pi), N)
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
        x_matrix[:,i] = 0.5.+asin.(sin.(ω .* @view s_matrix[:,i]))./pi # search function
    end

    # Compute outputs
    y_out = zeros(N)
    Threads.@threads for i in 1:N
        y_out[i] = model(@view x_matrix[i,:])
    end

    # Iterate over factors
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


        ft = (fft(y_reordered))[2:(N ÷ 2)]
        ys = abs2.(ft) .* inv(N)
        V = 2*sum(ys)
        Vi = 2*sum(ys[(1:M)*Int(ω)])
        Si = Vi/V

        # unskew the sensitivies
        lambda = 2*M/N
        sensitivites[i] = Si - (lambda / (1 - lambda)) * (1-Si)
    end

    return sensitivites
end


# This functon is an exact copy of the function above,
# except its not multi threaded, and is simply being used for benchmarking comparisons.
function _RBD_FAST_non_multi(model, N, group_encoding, frequencies, M=6)

    @assert maximum(frequencies) <= (N-1)/(2*M)

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
        # perms_matrix[i,:] = 1:N
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
        # s0_matrix[:,i] = collect(-pi:2*pi/(N-1):pi)
        s0_matrix[:,i] = range(-pi, stop=pi, length=N)
        # s0_matrix[:,i] = rand(Uniform(-pi, pi), N)
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
        x_matrix[:,i] = 0.5.+asin.(sin.(ω .* @view s_matrix[:,i]))./pi # search function
    end

    # Compute outputs
    y_out = zeros(N)
    Threads.@threads for i in 1:N
        y_out[i] = model(@view x_matrix[i,:])
    end

    # Iterate over factors
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


        ft = (fft(y_reordered))[2:(N ÷ 2)]
        ys = abs2.(ft) .* inv(N)
        V = 2*sum(ys)
        Vi = 2*sum(ys[(1:M)*Int(ω)])
        Si = Vi/V

        # unskew the sensitivies
        lambda = 2*M/N
        sensitivites[i] = Si - (lambda / (1 - lambda)) * (1-Si)
    end

    return sensitivites
end
