using ModelingToolkit
using MethodOfLines
using Symbolics
using NonlinearSolve
using DifferentialEquations
using LinearAlgebra
using GLMakie


## Load the model definition. This defines a PDESystem from ModelingToolkit
# and discretizes it
include("modeldef.jl")


## High-fidelity simulation

# define which parameter values to evaluate (based on a uniform grid)
Nvals = 50
# activation energy
EVals = range(5.5e11, 1.5e13, length=Nvals)
# pre-exponential coefficient
AVals = range(1.5e3, 9.5e3, length=Nvals)
# the pairs of act. energy and pre-exp. coeff. to evaluate
paramVals = Iterators.product(EVals, AVals)

# store high-fidelity results
# HFResults = Dict()
@assert typeof(problem.u0) <: Vector
HForder = length(problem.u0)
HFdata = zeros(HForder, length(paramVals))

# TODO look into reason for allocations:
# ~32GB of allocations for 2500 simulations of dim 1300
@time for (i, (EVal, AVal)) in enumerate(paramVals)
    # TODO tell somebody (slack/git) about `solve` not working
    # when the parameters are set like this:
    # paramProb = remake(problem, p=[E => EVal, A => AVal])
    paramProb = remake(problem, p=[AVal, EVal])
    paramSolution = solve(paramProb, NewtonRaphson())

    # HFResults[(EVal, AVal)] = paramSolution
    HFdata[:,i] = paramSolution.u
end


## Start generating surrogate model

# First POD

F = svd(HFdata)
# number of basis vectors to use
k = 100

# now DEIM, i.e. nonlinear function approximation

# the (max) number of interpolation indices to use
m = 30
PT = zeros(m, HForder)

for l in 1:m
    if l == 1
        r = F.U[:,1]
    else
        # TODO use a linear system solver?
        # TODO use views?
        # TODO use permutations/row selection instead of multiplying by PT?
        c = (PT[1:l-1, :] * F.U[:, 1:l-1]) \ (PT[1:l-1, :] * F.U[:, l])
        r = F.U[:,l] - F.U[:, 1:l-1] * c
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

# Offline construction of the nonlinear function

# TODO alternative to direct inversion?
Um = @view F.U[:,1:m]
mat1 = Um * inv(PT * Um)
mat2 = PT * Um # TODO this only works for nonlinear functions that are evaluated componentwise


## Find the discretized set of equations

system, timespan = symbolic_discretize(pdesys, discretization)
# simplified = structural_simplify(system) # not necessary (?) since it's a nonlinear system
discreteEqs = equations(system)

states = system.states
statesvec = [state for state in states]
RHSs = [eq.rhs for eq in discreteEqs] # TODO can you assume that the LHS is zero?

# TODO check whether the RHSs are linear
# if islinear(RHSs, statesvec)
#     jac = sparsejacobian(RHSs, statesvec)
#     nonlinear = nothing
# else
jac, nonlinear = semilinear_form(RHSs, statesvec)
# end


## Find out which components are used by the surrogate model

_, indices = findmax(PT, dims=2)
# the indices of the outputs of the nonlinearity that will be computed
outputIndices = sort!(unique(idx[2] for idx in indices))
# the indices of the inputs of the nonlinearity that are needed for calculations
inputVars = Set()
inputIndices = []
for idx in outputIndices
    expr = nonlinear[idx]
    vars = get_variables(expr)

    for v in vars
        # only register output variables, not parameters
        # TODO is this type too restrictive?
        if typeof(v) <: Term{Real, Nothing}
            push!(inputVars, v)
        end
    end
end

# is this even necessary? think about how you translate between high and lower order solutions
# TODO figure out how to compare inputVars elements to statesvec elements
# (their types are different)
# statesvec -> Num
# inputVars -> Term{Real, Nothing}
# both are displayed as YF[16,5] in the terminal

## Extract results (needs to be changed if edge aligned is used instead of center aligned discretization)

discrete_x = xmin:dx:xmax
discrete_y = ymin:dy:ymax

Nx = floor(Int64, (xmax - xmin) / dx) + 1
Ny = floor(Int64, (ymax - ymin) / dy) + 1

@variables T[1:Nx,1:Ny](t)

solT = map(1:length(sol.t)) do k
    # TODO extend this functionality to the mass fractions
    return reshape([sol[T[(i-1)*Ny+j]][k] for i in 1:Nx for j in 1:Ny], (Nx,Ny))
end


## Show results

# TODO extend this section to apply to the msas fractions
# find highest and lowest temperatures throughout entire simulation
# TODO rewrite using extrema function?
Tmin = minimum(map(minimum, solT))
Tmax = maximum(map(maximum, solT))

fig = Figure()
timeObs = Observable(1)
dataObs = @lift(solT[$timeObs])
# plot the initial frame
ax, hm = heatmap(fig[1, 1], discrete_x, discrete_y, dataObs, 
                colorrange=(Tmin, Tmax))
Colorbar(fig[1, 2], hm)

record(fig, "heatdiffusion2.mp4", 1:length(sol.t); framerate=8) do timestamp
    # update the time, this automatically updates the plot through the dataObs
    timeObs[] = timestamp
end
