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





# num_experiments = 50
# experiments = zeros(num_experiments, 8)
# for i in 1:num_experiments
#     experiments[i,:] = RBD_FAST_multi(sobol_function1, 10000, [1, 1, 1, 1, 2, 2, 2, 2], [11, 35])
# end
