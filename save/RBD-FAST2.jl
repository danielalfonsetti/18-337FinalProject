# https://www.sciencedirect.com/science/article/pii/S0951832005001444
# Page 168 in Global Sensitivty Analysis - The Primer

using FFTW

function RBD_FAST(model, N, M=6)
    """
    M = number of harmonics
    N = Number of samples to take
    """
    M = 6

    s0 = collect(-pi:2*pi/N:pi)

    # each factor i gets its own random permutation
    s = s0[randperm(N)]

    # uniform distribution between 0 and 1.
    # The function used can be different for each X.
    x=0.5.+asin.(sin.(s))./pi

    y = model(x) # Apply the model to each X vector

    ranks = sortperm(s)

    # reorder the Ys such that the corresponding values of the s
    # would be in increasing order.
    # (Note: This  does not necessarily mean the Ys are in increasing order)
    y_reordered = y[ranks]

    spectrum=(abs.(fft(y_reordered))).^2/N; # get the power density spectrum

    V1=2*sum(spectrum[2:M+1])
    V=sum(spectrum[2:N])

    S1 = V1/V
    return S1
end



function create_function_of_sobol(ps)
    """
    Create a function of sobol, parameterized by the the vector of ps.
    The function of sobol is widely used as a sensitivity analysis benchmark.
    https://www.sciencedirect.com/science/article/pii/S0951832005001444
    """

    g(x, p) = (abs(4*x-2)+p)/(1+p) # for 0 <= Xi <= 1 and ai >= 0
    gs = [k(x)=g(x,p) for p in ps] # create functions

    function f(xs)
        @assert length(gs) == length(xs)
        total = 1
        for (index, (g, x)) in enumerate(zip(gs, xs))
            println(index)
            total *= g(x)
        end
        return total
    end
    return f
end

ps1 = [0, 1, 4.5, 9, 99, 99, 99, 99]
sobol_function1 = create_function_of_sobol(ps)

sobol_function1([1, 2, 3, 4, 5, 6, 7, 8])


RBD_FAST(sobol_function1, 1000, M=6)


#
# ps2 = [0, 0, 0, 0, 0]
# sobol_function2 = create_function_of_sobol(ps2)
# res = sobol_function2([1, 100, 1, 1, 1])
#
#



# for (index, (val1, val2)) in enumerate(zip([1,2], [3,4]))
#     println("---")
#     println(index)
#     println(val1)
#     println(val2)
# end
#
# model(x) = x + 100
#
# RBD_FAST(model, 10)

 # e.g. if ranks is [1 3 2],
 # then element in position 1 of s gets moved to position 1
 # element in position of 2 of s gets moved to position 3
 # element in position 3 of s gets moved to position 2
