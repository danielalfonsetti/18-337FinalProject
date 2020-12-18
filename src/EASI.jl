"""
Code based on the theory presented in
    - Elmar Plischke (2010) "An effective algorithm for computing global
      sensitivity indices (EASI) Reliability Engineering & System Safety",
      95:4, 354-360. doi:10.1016/j.ress.2009.11.005

and the python implementation of EASI in python's Sensitivty Analysis Library ("SALib")
"""

using FFTW

function _permute_outputs(X::AbstractArray, Y::AbstractArray)
    """
    Triangular shape permutation of the precomputed inputs
    """
   permutation_index = sortperm(X) # non-mutating
   result = cat(permutation_index[1:2:end], reverse(permutation_index[2:2:end]), dims=1)
   return @view Y[result]
end


function _compute_first_order(permuted_outputs, M, N)
    ft = (fft(permuted_outputs))[2:(N ÷ 2)]
    ys = abs2.(ft) .* inv(N)
    V = 2*sum(ys)
    Vi = 2*sum(ys[(1:M)])
    Si = Vi/V
end


function _unskew_S1(S1::Number, M::Integer, N::Integer)
    """
    Unskew the sensivity index
    (Jean-Yves Tissot, Clémentine Prieur (2012) "Bias correction for the
    estimation of sensitivity indices based on random balance designs.",
    Reliability Engineering and System Safety, Elsevier, 107, 205-213.
    doi:10.1016/j.ress.2012.06.010)
    """
    λ = (2 * M)/N
    return S1 - (λ/(1-λ))*(1-S1)
end


function EASI(X::AbstractArray, Y::AbstractArray, M::Integer=10)
    """
    param: X
        pre-sampled model inputs.
        Each row corresponds to one input vector.
    param: Y
        A vector of precomputed model outputs.
    param: M
        Maximum harmonic of the input frequency
        (which is always w=1 in this method) for which the output
        power spectrum is analyzed for.
    return:
        A 2-tuple of vectors of size K, where K is the number of factors
        (i.e. the number of columns in X).

        The first vector contains the estimated fraction of variance for
        each model input contributes to the output.

        The second vector is a bias-crrected version of the first.
    """

    # K is the number of variables, N is the number of simulations
    N, K = size(X, 1), size(X, 2)
    sensitivites = zeros(K)
    sensitivites_c = zeros(K)

    Threads.@threads for i in 1:K
        Xi = @view X[:, i]

        Y_reordered = _permute_outputs(Xi, Y)
        S1 = _compute_first_order(Y_reordered, M, N)

        S1_C = _unskew_S1(S1, M, N) # get bias-corrected version
        sensitivites[i] = S1
        sensitivites_c[i] = S1_C
    end

    return sensitivites, sensitivites_c
end



# This functon is an exact copy of the function above,
# except its not multi threaded, and is simply being used for benchmarking comparisons.
function _EASI_non_multi(X::AbstractArray, Y::AbstractArray, M::Integer=10)
    """
    param: X
        pre-sampled model inputs.
        Each row corresponds to one input vector.
    param: Y
        A vector of precomputed model outputs.
    param: M
        Maximum harmonic of the input frequency
        (which is always w=1 in this method) for which the output
        power spectrum is analyzed for.
    return:
        A 2-tuple of vectors of size K, where K is the number of factors
        (i.e. the number of columns in X).

        The first vector contains the estimated fraction of variance for
        each model input contributes to the output.

        The second vector is a bias-crrected version of the first.
    """

    # K is the number of variables, N is the number of simulations
    N, K = size(X, 1), size(X, 2)
    sensitivites = zeros(K)
    sensitivites_c = zeros(K)

    for i in 1:K
        Xi = @view X[:, i]

        Y_reordered = _permute_outputs(Xi, Y)
        S1 = _compute_first_order(Y_reordered, M, N)

        S1_C = _unskew_S1(S1, M, N) # get bias-corrected version
        sensitivites[i] = S1
        sensitivites_c[i] = S1_C
    end

    return sensitivites, sensitivites_c
end
