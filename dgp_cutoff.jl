using Random, Distributions

function generateX(n=1000)
    x = rand(Uniform(0, 5), n)
    sort!(x)
    y = zeros(n)
    h = 1
    y[x .> 4] = 

end
