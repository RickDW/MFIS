# Multi-Fidelity Importance Sampling

This repository will soon contain an implementation of multi-fidelity importance sampling (MFIS) in Julia. This method is taken from the paper ["Multi-fidelity Importance Sampling" by Peherstorfer et al.](https://www.sciencedirect.com/science/article/pii/S004578251500393X).

Importance sampling allows you to speed up certain statistical estimators by decreasing their variance. To reduce the variance of these estimators and get better results, you need to define a sampling distribution. MFIS constructs one using a surrogate model.

## Surrogate generation

Surrogate models can be constructed in countless different ways but can be time consuming to implement. This package will allow you to automatically construct surrogates from partial differential equation systems using [ModelingToolkit.jl](mtk.sciml.ai). This allows you to use all of the other analysis and simulation tools from the SciML stack without much hassle.

## Importance sampling

For a quick overview of importance sampling and how MFIS fits into it you can take a look at the presentation pdf in this repository.