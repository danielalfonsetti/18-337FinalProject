using KernelDensity
Y = [i^2 for i in 1:100]
grid = range(minimum(Y), maximum(Y), length=100)
fy = kde(Y, grid).density
print(fy[1:5])
