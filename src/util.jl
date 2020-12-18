function create_function_of_sobol(ps)
    """
    Create a function of sobol, parameterized by the the vector ps.
    The function of sobol is widely used as a sensitivity analysis benchmark.
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

function ishigami(X)
    """
        Ishigami function. Takes a vector input, returns a scalar.
    """
    X1, X2, X3, X4 = X
    Y = sin(X1) + 7*sin(X2)^2 + 0.1*X3^4*sin(X1)
    return Y
end


function create_ishigami_samples(sample_size)
    """
        Sample input vectors for the ishigami function,
        evaluate the ishigami function on them,
        and returns a vector of the scalar outputs.
    """
    X = rand(Uniform(-pi, pi), sample_size, 4)
    Y = zeros(sample_size)
    for sample_num in 1:sample_size
        ishigami_input = X[sample_num, :]
        Y[sample_num] = ishigami(ishigami_input)
    end

    return X, Y
end
