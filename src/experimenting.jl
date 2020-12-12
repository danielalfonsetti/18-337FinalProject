

f(x)  = x + 10


function h(k, x, y)
    return k.(x) .+ y
end



h(f, [2, 2], [10, 30])
