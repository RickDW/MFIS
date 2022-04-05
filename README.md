# Multi-Fidelity Importance Sampling

This repository will soon contain an implementation of multi-fidelity importance sampling (MFIS) in Julia. This method is taken from the paper ["Multi-fidelity Importance Sampling" by Peherstorfer et al.](https://www.sciencedirect.com/science/article/pii/S004578251500393X).

Importance sampling allows you to speed up certain statistical estimators by decreasing their variance. To reduce the variance of these estimators and get better results, you need to define a sampling distribution. MFIS constructs one using a surrogate model.

## Surrogate generation

Surrogate models can be constructed in countless different ways but can be time consuming to implement. This package will allow you to automatically construct surrogates from partial differential equation systems using [ModelingToolkit.jl](https://mtk.sciml.ai). This allows you to use all of the other analysis and simulation tools from the SciML stack without much hassle. 

The type of surrogate models that will be supported are Proper Orthogonal Decomposition (POD) and Discrete Empirical Interpolation Method (DEIM). POD allows you to reduce the order of linear systems such as ODEs and root finding problems. It can handle nonlinearities but not with the same computational efficiency as DEIM. Combining POD and DEIM results in a relatively simple method that applies to lots of system types.

## Importance sampling

For a quick overview of importance sampling and how MFIS fits into it you can take a look at the presentation pdf in this repository.
