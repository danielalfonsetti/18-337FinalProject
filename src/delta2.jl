using KernelDensity, StatsBase, Trapz, Distributions

# https://www.sciencedirect.com/science/article/pii/S0377221712008995
# https://github.com/SALib/SALib/blob/master/src/SALib/analyze/delta.py


function calc_delta(X, Y)
        # Make sure Ys are not identical, otherwise KDE will be undefined.
        # If Ys are identical, then we know X and Y are independent, so return 0
        if all(Y[1] .== Y)
                println("Warning - Ys are identical. Returning 0 for delta estimate.")
                return 0
        end

        N = length(Y) # Number of simulations
        @assert length(Xi) == N # Length of Y should equal length of X

        # Create number of classes and class cutoffs.
        exp = (2 / (7 + tanh((1500 - N) / 500)))
        M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
        class_cutoffs =  range(0, N, length=M+1) # class cutoffs.


        # quadrature points.
        # Length should be a power of 2 for efficient FFT in kernel density estimates.
        Ygrid = range(minimum(Y), maximum(Y), length=100)
        Ygrid_collect = collect(Ygrid)
         # Model pdf of Y using KDE, kde uses normal kernel by default

        println(Ygrid[1:5])
        fy = kde(Y, Ygrid) # eq 23.1
        println(fy.density[1:5])


        # Get probability of each y in Ygrid.
        x_rank = competerank(Xi) # Does what scipy.stats rankdata does. If tie, all tied values get same rank.
        d_hat = 0 # the delta estimator.

        # Iterate over each class
        for j in 1:length(class_cutoffs)-1

                # get X and Y indicies for samples that are in this class (where
                # class designation is based on the X value)
                condition(x) = (x > class_cutoffs[j]) ==  (x <= class_cutoffs[j+1])
                in_class_indices = findall(condition, x_rank)
                number_in_class = length(in_class_indices)

                # Get the the subset of Y values in this class
                in_class_Y = Y[in_class_indices]

                # get the separation between the total y pdf and the condition y pdf induced by this class
                fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated condition distribution of y
                curve_diff = abs.(fy.density .- fyc.density) # eq 24

                # Use trapezoidal rule to estimate the difference between the curves.
                class_separation = trapz(Ygrid, curve_diff) # eq 25

                # Increment estimator
                d_hat += number_in_class * class_separation # eq 26
        end

        d_hat = d_hat/(2*N)
        return d_hat
end


Xi = [-2.39677772,  -2.39884787, 1.60521748, -2.64049256, 1.89295698]
Y = [5.79724193, 2.71923391, 1.1565645 , 0.25716799, 4.40692813]
println(calc_delta(Xi, Y))



# delta1 = 0.208
# delta2 = 0.391
# delta3 = 0.156
# delta4 = 0.060


function ishigami(X)
        """
                Ishigami function. Takes a vector, returns a scalar.
                Saltelli et al., 2004
        """
        X1, X2, X3, X4 = X
        part1 = 0.1*X3^4*sin(X1)
        Y = sin(X1) + 7*sin(X2)^2 + part1
        return Y
end



function test_function(sample_size)
        # sample_sizes = 2 .^ [i for i in 1:10]
        # delta_results = zeros(10, 4) #
        # for (i, sample_size) in enumerate(sample_sizes)
        # end

        X = rand(Uniform(-pi, pi), sample_size, 4)
        Y = zeros(sample_size)
        for sample_num in 1:sample_size
                ishigami_input = X[sample_num, :]
                Y[sample_num] = ishigami(ishigami_input)
        end

        deltas = zeros(4)
        for column in 1:4
                Xi = @view X[:, column]
                deltas[column] = calc_delta(Xi, Y)
        end
        return deltas
end

function test_function2(sample_size)
        # sample_sizes = 2 .^ [i for i in 1:10]
        # delta_results = zeros(10, 4) #
        # for (i, sample_size) in enumerate(sample_sizes)
        # end

        X = rand(Uniform(-pi, pi), sample_size, 4)
        Y = zeros(sample_size)
        for sample_num in 1:sample_size
                ishigami_input = X[sample_num, :]
                Y[sample_num] = ishigami(ishigami_input)
        end

        return calc_delta(X[:,4], Y, 110, 10)
end



# Try to recreate figure 3. # https://www.sciencedirect.com/science/article/pii/S0377221712008995#e0040
sample_sizes = 2 .^ [i for i in 9:14]
fig3_deltas = zeros(length(sample_sizes))
fig3_adj_deltas = zeros(length(sample_sizes))
fig3_delta_conf = zeros(length(sample_sizes))
for (index, sample_size) in enumerate(sample_sizes)
        fig3_deltas[index] = test_function2(sample_size)
end


scatter(
        sample_sizes,
        fig3_deltas,
        label="Deltas",
        xaxis=:log,
        xticks = (sample_sizes, sample_sizes),
        legend=:topright,
)



using KernelDensity
Y = [i^2 for i in 1:100]
grid = range(minimum(Y), maximum(Y), length=20)
fy = kde(Y, grid).density
print(fy[1:5])


using Plots
plot(kde(Y, grid).density)
