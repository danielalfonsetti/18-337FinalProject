using KernelDensity, StatsBase, Trapz, Random, Distributions

using NumericalIntegration

"""
The code here is based on the theory presented in

        Plischke, E., E. Borgonovo, and C. L. Smith (2013). "Global
        sensitivity measures from given data." European Journal of
        Operational Research, 226(3):536-550

and the python implementation of delta-moment in python's Sensitivty Analysis Library ("SALib")
"""

function calc_delta(Xi, Y, Ygrid, class_cutoffs)

        # Make sure Ys are not identical, otherwise KDE will be undefined.
        # If Ys are identical, then we know X and Y are independent, so return 0
        if all(Y[1] .== Y)
                println("Warning - Ys are identical. Returning 0 for delta estimate.")
                return 0
        end

        N = length(Y) # Number of simulations
        @assert length(Xi) == N # Length of Y should equal length of X


        k = KernelDensity.kde(Y) # defaults are kernel = normal and bandwidth = Silverman which match SALib
        fy = pdf(k, Ygrid)
        model_input_ranks = competerank(Xi)

        d_hat = 0
        for j = 1:length(class_cutoffs) - 1
            ix = findall((model_input_ranks .> class_cutoffs[j]) .& (model_input_ranks .<= class_cutoffs[j + 1]))
            nm = length(ix)

            tmp = Y[ix]
            if length(tmp) == 0
                    continue
            end
            k = KernelDensity.kde(Y[ix]) # defaults are kernel = normal and bandwidth = Silverman which match SALib
            fyc = pdf(k, Ygrid)
            d_hat += (nm / (2 * N)) * integrate(abs.(fy - fyc), sort(Ygrid, rev = true))
        end
        return d_hat


        #  # Model pdf of Y using KDE, kde uses normal kernel by default
        # fy = pdf(kde(Y), Ygrid) # eq 23.1
        #
        # # Get probability of each y in Ygrid.
        # x_ranks = ordinalrank(Xi) # Does what scipy.stats rankdata does. If tie, all tied values get same rank.
        # d_hat = 0 # the delta estimator.
        #
        # # Iterate over each class
        # weighted_class_seps = zeros(length(class_cutoffs)-1)
        # Threads.@threads for j in 1:length(class_cutoffs)-1
        #
        #         # get X and Y indicies for samples that are in this class (where
        #         # class designation is based on the X value)
        #         # condition(x) = (x .> class_cutoffs[j]) .&  (x .<= class_cutoffs[j+1])
        #         # in_class_indices = findall(condition, x_rank)
        #         condition(x) = (x > class_cutoffs[j]) ==  (x <= class_cutoffs[j+1])
        #         in_class_indices = findall(condition, x_ranks)
        #         number_in_class = length(in_class_indices)
        #
        #         # Get the the subset of Y values in this class
        #         in_class_Y = Y[in_class_indices]
        #
        #         # get the separation between the total y pdf and the condition y pdf induced by this class
        #         # fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated conditional distribution of y (using normal kernel)
        #         # curve_diff = abs.(fy.density .- fyc.density) # eq 24
        #
        #         # if ~all(t == in_class_Y[1] for t in in_class_Y)
        #         #         # get the separation between the total y pdf and the condition y pdf induced by this class
        #         #         # fyc = kde(in_class_Y, Ygrid).density # eq 23.2 - Estimated condition distribution of y
        #         #         # curve_diff = abs.(fy .- fyc) # eq 24
        #         # else
        #         #         curve_diff = fy
        #         # end
        #         if length(in_class_Y) == 0
        #                 continue
        #         end
        #
        #         fyc = pdf(kde(in_class_Y), Ygrid)
        #         class_separation = integrate(abs.(fy-fyc), sort(Ygrid, rev= true))
        #
        #         # Use trapezoidal rule to estimate the difference between the curves.
        #         # class_separation = trapz(Ygrid, curve_diff) # eq 25
        #
        #         # Increment estimator
        #         weighted_class_seps[j] = number_in_class * class_separation # eq 26
        #         # d_hat += number_in_class * class_separation # eq 26
        # end
        #
        # d_hat = sum(weighted_class_seps)
        # d_hat = d_hat/(2*N)
        # return d_hat
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
