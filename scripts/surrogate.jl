using GLMakie
using LinearAlgebra


## Define the function to be approximated, plus the data & parameter domains

function nlfunc(x, μ)
    return @. (1-x)*cos(3π * μ * (x+1)) * exp(-1*(1+x) * μ)
end

# dimensionality of function output, i.e. the grid size
n = 100
# number of samples of the output space
s = 51

xgrid = range(-1, 1, length=n)
μgrid = range(1, π, length=s)


## Collect data
X = zeros(n, s)

for (i, μ) in enumerate(μgrid)
    X[:, i] = transpose(nlfunc.(xgrid, μ))
end

# subtract the mean componentwise
# Xc = X .- mean(X, dims=2)


## Function visualization
fig = Figure()
ax = Axis(fig[1,1])
for (i, μ) in enumerate(μgrid)
    lines!(ax, xgrid, X[:,i], label="μ = $μ")
end
axislegend(ax)
display(fig)


## POD, i.e. SVD
F = svd(X)
coeff = Diagonal(F.S) * F.Vt
U = F.U # * Diagonal(F.S)
Nsv = length(F.S)


# calculate the projection error for using only the first x basis vectors
rerror = zeros(Nsv)
currentSum = 0.0
for i in Nsv:-1:1
    rerror[i] = currentSum
    currentSum += F.S[i]^2
end

figfact = Figure()
ax1 = Axis(figfact[1,1], title="Singular values of $s snapshots",
        yscale=log10)
ax2 = Axis(figfact[1,2], title="Error")
ax3 = Axis(figfact[2,:], title="POD/DEIM basis vectors")

# show all of the singular values
lines!(ax1, F.S, label="Singular values")
lines!(ax2, rerror, label="Error")

# show the first d POD basis vectors
d = 6
for i in 1:d
    lines!(ax3, xgrid, U[:,i], label="Basis vector $i")
end
axislegend(ax3)
display(figfact)


## Display POD approximation

k = 30
Uk = U[:,1:k]
μtest = 1.17

figtest = Figure(title="POD approximations using k=$k basis vectors")
ax = Axis(figtest[1,1])

lines!(ax, xgrid, nlfunc.(xgrid, μtest), label="Exact")
lines!(ax, xgrid, Uk * transpose(Uk) * nlfunc.(Uk * transpose(Uk) * xgrid, μtest), label="Approx.")
axislegend(ax)

display(figtest)



## DEIM, i.e. nonlinear function approximation

# the (max) number of interpolation indices to use
m = 30
PT = zeros(m, n)

for l in 1:m
    if l == 1
        r = U[:,1]
    else
        # TODO use a linear system solver?
        # TODO use views?
        # TODO use permutations/row selection instead of multiplying by PT?
        c = (PT[1:l-1, :] * U[:, 1:l-1]) \ (PT[1:l-1, :] * U[:, l])
        r = U[:,l] - U[:, 1:l-1] * c
    end

    # find the index of the maximum element in r
    _, pl = findmax(abs.(r))
    PT[l,pl] = 1.0
end

# test whether there are m unique indices in PT
begin
    vals, indices = findmax(PT, dims=2)
    @assert length(unique(idx[2] for idx in indices)) == size(PT, 1) # = m
end


## Offline construction of the nonlinear function

# TODO alternative to direct inversion?
Um = U[:,1:m]
mat1 = Um * inv(PT * Um)
mat2 = PT * Um # TODO this only works for nonlinear functions that are evaluated componentwise


## Visualize reconstructions based on m interpolation points

μtest = 1.17

figtest = Figure(title="DEIM approximations using m=$m interpolation points")
ax = Axis(figtest[1,1])
lines!(ax, xgrid, nlfunc.(xgrid, μtest), label="Exact")
lines!(ax, xgrid, mat1 * nlfunc.(mat2 * transpose(Um) * xgrid, μtest), label="Approx.")
axislegend(ax)
display(figtest)

