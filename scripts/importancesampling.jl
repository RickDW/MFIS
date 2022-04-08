# this program demonstrates importance sampling for a simple univariate distribution
using Distributions
using GLMakie
using LaTeXStrings


## Represent normal dist. using latex

function reprDist(dist::Normal)
    return L"N(\cdot | \mu = %$(dist.μ), \sigma = %$(dist.σ))"
end



## Define nominal and sampling distributions + limit state function

nomDist = Normal(0, 1) # target distribution
sampDists = Dict(
    :sampling1 => Normal(4, 2),
    :sampling2 => Normal(-4, 4),
)

boundary = 3
limstate(x) = x > boundary # true indicates failure


## Define importance weight function


## Calculate true failure probability
pfail = 1 - cdf(nomDist, boundary)


## Perform importance sampling

function simpleIS(nominalD::Distribution, samplingD::Distribution, func, Nsamples, Ntrials=1)
    # generate samples from the sampling distribution
    samples = rand(samplingD, (Nsamples, Ntrials))

    # calculate importance weights
    weights = pdf.(nominalD, samples) ./ pdf.(samplingD, samples)

    # calculate the quantity of interest, i.e. whatever you want to take the expectation of
    evaluations = func.(samples)

    # calculate the approximate expectation per trial (using IS)
    # returns a 1xNtrials matrix
    values = mean(evaluations .* weights, dims=1)

    return values
end


Ntrials = 20
NsamplesArr = [100, 500, 1000, 5000, 10000, 50000, 100000]
results = Dict(dist => zeros(Ntrials, length(NsamplesArr)) for dist in keys(sampDists))

for (i, Nsamples) in enumerate(NsamplesArr)
    for distKey in keys(sampDists)
        estimates = simpleIS(nomDist, sampDists[distKey], limstate, Nsamples, Ntrials)
        results[distKey][:, i] = estimates
    end
end



#
## Plot the distributions
#


plotProps = Dict(
    :nominal => Dict{Symbol,Any}(
        :color => Cycled(1),
        :label => reprDist(nomDist)
    ),
    :sampling1 => Dict{Symbol,Any}(
        :color => Cycled(2),
        :label => reprDist(sampDists[:sampling1])
    ),
    :sampling2 => Dict{Symbol,Any}(
        :color => Cycled(3),
        :label => reprDist(sampDists[:sampling2])
    )
)


fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98))

# axdens will show the pdf's of the densities
axdens = Axis(fig[1,1], xlabel="value", ylabel="density")

# axexper1/2 will show a representative sample of IS runs
axexper1 = Axis(fig[2,1], xlabel="N_samples", ylabel="p_fail", xscale = log10)
axexper2 = Axis(fig[2,2], xlabel="N_samples", xscale=log10)
linkyaxes!(axexper1, axexper2)

# keep track of plotted items so they can be added to the legend
plotItems = []
addItem(it) = push!(plotItems, it)


xplot = -7:0.1:7

# probability densities
lines!(
    axdens, xplot, nomDist, 
    label="nominal dist.", 
    color=plotProps[:nominal][:color]
) |> addItem
lines!(
    axdens, xplot, sampDists[:sampling1], 
    label="good sampling dist.", 
    color=plotProps[:sampling1][:color]
) |> addItem
lines!(
    axdens, xplot, sampDists[:sampling2], 
    label="bad sampling dist.", 
    color=plotProps[:sampling2][:color]
) |> addItem
# vlines!(axdens, boundary, label="failure boundary", color=Cycled(4))
vspan!(
    axdens, boundary, maximum(xplot), 
    label="failure region", color=(:pink, 0.3)
) |> addItem
xlims!(axdens, extrema(xplot))

for (sampD, axis) in zip([:sampling1, :sampling2], [axexper1, axexper2])
    # turn the results matrices into vectors such that
    # (scatterx[i], scattery[i]) = (N_samples_i, p_fail_i)
    NsamplesMat = repeat(transpose(NsamplesArr), Ntrials, 1)
    scatterx = reshape(NsamplesMat, (:))
    scattery = reshape(results[sampD], (:)) # stack columsn one under the other

    # plot the vectors
    scatter!(
        axis, scatterx, scattery, color=plotProps[sampD][:color],
        markersize=3
    )
    pfailline = hlines!(axis, pfail, color=Cycled(5), label="true failure prob.") |> addItem
end
# remove last failure lines to prevent duplicates
# workaround for `unique=true` not working right for the legend
pop!(plotItems)

hideydecorations!(axexper2) # hide y axis decorations

# add a legend for the entire figure
# unique=true removes repeated legend items (i.e. the true failure prob. line)
leg = Legend(fig[1,2], plotItems, [it.label[] for it in plotItems], unique=true)
leg.tellwidth = false

# add a super title
title = Label(fig[1,:,Top()], "Failure probability estimation using Importance Sampling",
    textsize=20, valign=:bottom, padding=(0,0,10,0))
    

display(fig)
save("ISdemo.png", fig)