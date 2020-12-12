# https://github.com/Sarah-Ju/RBD-Fast
# https://www.sciencedirect.com/science/article/pii/S0951832005001444
# Page 168
# https://www.sciencedirect.com/science/article/pii/S0951832005001444

# https://salib.readthedocs.io/en/latest/api.html#rbd-fast-random-balance-designs-fourier-amplitude-sensitivity-test
# https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/rbd_fast.html#analyze
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html
# https://en.wikipedia.org/wiki/Spectral_density
# https://en.wikipedia.org/wiki/Spectral_density_estimation
# https://www.scopus.com/record/display.uri?eid=2-s2.0-84953649167&origin=inward&txGid=e0075e43de1008dc516f212917e43189

# values up to:  (N−1)/2M

M = 3
N = 1000
ω = 1
@assert ω <= (N-1)/(2*M)

s0 = collect(-pi:2*pi/N:pi)
# each factor i gets its own random permutation
s = s0[randperm(N)]
# perm = [5 10 8 2 6 7 4 9 1 3]
# s = vec(s0[perm])

# uniform distribution between 0 and 1.
# The function used can be different for each X.
x=0.5.+asin.(sin.( ω .* s))./pi
y = x .+ 10 # Apply the model
 # rank (rank 1 = smallest). Why are we not ordering by x?
ranks = sortperm(s)

# reorder the Ys such that the corresponding values of the x
# would be in increasing order. (Note: Ys don't necessarily have to be in increasing order)
y_reordered = y[ranks]

using FFTW
# f = fft(y_reordered)
# signal_power = f[1:()]

# sp = (abs)
spectrum = ((abs.(fft(y_reordered))) .^ 2) ./ N
V=sum(spectrum[2:N])
harmonics = [ω*i for i=1:M]
V1 = 2*sum(spectrum[harmonics.+1])
S1 = V1/V
S1


# Sp = abs.(y_fft[2:], int((N + 1) / 2)))

# Sum the power at harmonics of omega:


 # e.g. if ranks is [1 3 2],
 # then element in position 1 of s gets moved to position 1
 # element in position of 2 of s gets moved to position 3
 # element in position 3 of s gets moved to position 2
