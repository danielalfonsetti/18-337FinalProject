# https://www.sciencedirect.com/science/article/pii/S0951832009002579
# Code based on: https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/rbd_fast.html#analyze

using FFTW


function permute_outputs(X, Y)
    """
    Triangular shape permutation of the precomputed inputs

    References
    ----------
    .. [2] Elmar Plischke (2010) "An effective algorithm for computing global
          sensitivity indices (EASI) Reliability Engineering & System Safety",
          95:4, 354-360. doi:10.1016/j.ress.2009.11.005

    """

    # rank (rank 1 = smallest)
   permutation_index = sortperm(X) # non-mutating
   result = cat(permutation_index[1:2:end], reverse(permutation_index[2:2:end]), dims=1)
   return Y[result]
end


function compute_first_order(permuted_outputs, M)
    """
    References
    ----------
    .. [2] Elmar Plischke (2010) "An effective algorithm for computing global
          sensitivity indices (EASI) Reliability Engineering & System Safety",
          95:4, 354-360. doi:10.1016/j.ress.2009.11.005

    """
    power_spectrum = abs.(fft(permuted_outputs)).^2

    # Estimated total variance
    V = sum(power_spectrum[2:end]) # don't take 0 frequency

    # Estimated conditional variance
    VC = sum(power_spectrum[2:M+1]) # don't take 0 frequency

    S = VC/V
end


function unskew_S1(S1, M, N)
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


function my_analyze(X, Y, M=10)
    """

    Parameters :
        X: precomputed model inputs. (matrix)
        Y: model outputs
        M: Maximum harmonic of the input frequency
        (which is always w=1 in this method) for which the output
        power spectrum is analyzed for.
    Returns :

    """

    # K is the number of variables, N is the number of simulations
    N, K = size(X, 1), size(X, 2)
    sensitivites = zeros(K)
    sensitivites_c = zeros(K)
    for i in 1:K
        Y_reordered = permute_outputs(X[:, i], Y)
        S1 = compute_first_order(Y_reordered, M)
        S1_C = unskew_S1(S1, M, N) # get bias-corrected version
        sensitivites[i] = S1
        sensitivites_c[i] = S1_C
    end
end



ps1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(ps1)

runs = 150
for N in [100, 300, 1000, 3000, 10000, 30000]
    for run in 1:runs
        X = rand(N, 8)
    end
end

# X = 1:1:10
# spectrum = abs.(fft(X)).^2
# V1 = sum(spectrum[1:M])
# V = sum(spectrum[1:end])
#
# W = V1/V
#
#
# abs.compute_first_order(X, 1)
#
#
# Y = X.*10
# my_analyze(X, Y)
#
#
# res = periodogram(X)

#
# M = 6
# N = 10
# s0 = collect(-pi:2*pi/N:pi)
# # each factor i gets its own random permutation
# # s = s0[randperm(N)]
# perm = [5 10 8 2 6 7 4 9 1 3]
# s = vec(s0[perm])
#
# # uniform distribution between 0 and 1.
# # The function used can be different for each X.
# x=0.5.+asin.(sin.(s))./pi
# y = x .+ 10 # Apply the model
#
#
#  # rank (rank 1 = smallest). Why are we not ordering by x?
# ranks = sortperm(s)
#
# # reorder the Ys such that the corresponding values of the x
# # would be in increasing order. (Note: Ys don't necessarily have to be in increasing order)
# y_reordered = y[ranks]
#
# using FFTW
# spectrum=(abs.(fft(y_reordered))).^2/N; # get the fourier transform of the reordered Ys
#
# V1=2*sum(spectrum[2:M+1])
# V=sum(spectrum[2:N])
#
# S1 = V1/V

 # e.g. if ranks is [1 3 2],
 # then element in position 1 of s gets moved to position 1
 # element in position of 2 of s gets moved to position 3
 # element in position 3 of s gets moved to position 2
