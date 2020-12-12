
using KernelDensity

# https://www.sciencedirect.com/science/article/pii/S0377221712008995
# https://github.com/SALib/SALib/blob/master/src/SALib/analyze/delta.py

function calc_delta(Y)
        """Plischke et al. (2013) delta index estimator (eqn 26) for d_hat."""

        # equal frequency partition
        N = length(Y) # Number of simulations
        exp = (2 / (7 + tanh((1500 - N) / 500)))
        M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
        m =  range(0, N, length=M+1) # class cutoffs
        Ygrid = range(minimum(Y), maximum(Y), length=100)

        N = length(Y)
        fy = kde(Y) # using normal kernel by default
        return fy

end

##############
#
# https://www.sciencedirect.com/science/article/pii/S0377221712008995
#

# Xi = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.3, 0.1]
# Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]

Xi = [0.9, 0.8, 0.7, 0.4, 0.2, 0.3, 0.1, 0.6, 0.5]
Y = [100, 200, 300, 400, 500, 600, 700, 800, 900]

# equal frequency partition
N = length(Y) # Number of simulations
exp = (2 / (7 + tanh((1500 - N) / 500)))
M = Integer(round(min(Integer(ceil(N^exp)), 48))) # Number of classes
m =  range(0, N, length=M+1) # class cutoffs
Ygrid = range(minimum(Y), maximum(Y), length=100)


# Number of simulations
N = length(Y)

 # Model pdf of Y using KDE.
 # kde uses normal kernel by default, and 2048 points
fy = kde(Y) # eq 23.1

# Get probability of each y in Ygrid.


# Xi[sortperm(Xi)]
x_rank = sortperm(Xi)

d_hat = 0


# Get rankings in the first set of Xs
j = 1
########
# in for loop
# for j in 1:length(m)-1
########
ranks_in_class = x_rank[(x_rank .> m[j]) .== (x_rank .<= m[j+1])]

number_in_class = length(ix)

# Get the Y that correspond to the chosen rankings (using indexing)
Y_ix = Y[ranks_in_class]
Y_ix


if ~all(t == Y_ix[1] for t in Y_ix)
        fyc = kde(Y_ix) # eq 23.2 - Estimated condition distribution of y
        curve_diff = abs(fy - fyc)
else
        curve_diff = abs_fy

end


 # Use trapezoidal rule to estimate the difference between the curves.
 # eq (25)
class_separation = trapz(curve_diff)
d_hat += nm * class_separation

########
# after for loop
########

d_hat *= (1/(2*N))

# https://github.com/JuliaStats/KernelDensity.jl


## -

# Iterate over classes.
# For each class, get the rankins
d_hat = 0
for j in 1:length(m)-1
        println(j)
        ix = x_rank[(x_rank .> m[j]) .== (x_rank .< m[j+1])]
        class_amount = sum(x_rank)
        class_Y = Y[ix]
end
