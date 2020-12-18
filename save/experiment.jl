using KernelDensity
Y = [i^2 for i in 1:100]
my_grid = range(minimum(Y), maximum(Y), length=100)
fy = kde(Y, my_grid).density
fy[1:5]



fy2 =  pdf(kde(Y), my_grid)
fy2[1:5]


using KernelDensityEstimate
p100 = kde!([randn(50);10.0.+2*randn(50)])


# https://github.com/JuliaStats/KernelDensity.jl
# https://github.com/JuliaRobotics/KernelDensityEstimate.jl
# https://github.com/panlanfeng/KernelEstimator.jl




using KernelDensity
Y = [i for i in 1:10]
grid = range(minimum(Y), maximum(Y), length=10)
fy = pdf(kde(Y), grid)
print(fy)
