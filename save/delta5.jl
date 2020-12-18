using KernelDensity, StatsBase, Trapz

# https://www.sciencedirect.com/science/article/pii/S0377221712008995
# https://github.com/SALib/SALib/blob/master/src/SALib/analyze/delta.py

# Xi = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.3, 0.1]
# Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]

Xi = [0.9, 0.8, 0.7, 0.4, 0.2, 0.3, 0.1, 0.6, 0.5]
Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]

@assert length(Xi) == length(Y)

N = length(Y) # Number of simulations
exp = (2 / (7 + tanh((1500 - N) / 500)))
M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
class_cutoffs =  range(0, N, length=M+1) # class cutoffs.
Ygrid = range(minimum(Y), maximum(Y), length=2048) # quadrature points. Length should be a power of 2 for FFT in kernel density estimates.


 # Model pdf of Y using KDE, kde uses normal kernel by default
fy = kde(Y, Ygrid) # eq 23.1

# Get probability of each y in Ygrid.
x_rank = competerank(Xi) # Does what scipy.stats rankdata does. What happens if Xis are tied??
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
        global d_hat += number_in_class * class_separation # eq 26
end

d_hat = d_hat/(2*N)



# if ~all(t == Y_ix[1] for t in Y_ix)
#         # get the separation between the total y pdf and the condition y pdf induced by this class
#         fyc = kde(in_class_Y, Ygrid) # eq 23.2 - Estimated condition distribution of y
#         curve_diff = abs.(fy.density .- fyc.density) # eq 24
# else
#         curve_diff = abs_fy.density
# end
