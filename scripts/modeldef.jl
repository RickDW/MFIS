# This program creates a symbolically defined model of a hydrogen flame
using ModelingToolkit
using MethodOfLines
using DomainSets
using NonlinearSolve


## Define variables and derivative operators

# the two spatial dimensions
@parameters x y A E
# the mass fractions MF (respectively H₂, O₂, and H₂O), temperature,
# and the source terms of the mass fractions and temperature
@variables YF(..) YO(..) YP(..) T(..) 
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

# laplace operator
∇²(u) = Dxx(u) + Dyy(u)


## Define domains

xmin = ymin = 0.0 # cm
xmax = 1.8 # cm
ymax = 0.9 # cm

domains = [
    x ∈ Interval(xmin, xmax),
    y ∈ Interval(ymin, ymax),
]


## Define parameters

# turn the reaction on/off (without this, only diffusion and convection are left)
reaction = true
# diffusivity coefficient (for temperature and mass fractions)
κ =  2.0 #* cm^2/s
# constant (divergence-free) velocity field
U = [50.0, 0.0] #.* cm/s
# density of the mixture
ρ = 1.39e-3 #* g/cm^3
# molecular weights (respectively H₂, O₂, and H₂O)
W = [2.016, 31.9, 18.0] #.* g/mol
# stoichiometric coefficients (positive for reactants, negative for product)
ν = [2, 1, -2]
# heat of the reaction
Q = 9800 #* K
# universal gas constant
R = 8.314472 * 100 #* 1e-2 J/mol/K -> because J = Nm = 100 N cm
# y coordinates of the inlet area on the left side of the domain
inlety = (0.3, 0.6) # mm

# pre-exponential factor in source term
A_default = 5.5e11 # dimensionless
# activation energy in source term
E_default = 5.5e13 * 100 # 1e-2 J/mol -> J = Nm = 100 N cm


## Define the model

function sourceTerm(species::Int)
    return -ν[species] * (W[species] / ρ) * 
        (ρ * YF(x,y) / W[1]) ^ ν[1] *
        (ρ * YO(x,y) / W[2]) ^ ν[2] * 
        A * exp(-E / R / T(x,y))
end


if reaction
    # diffusion, convection, and source terms for the mass fractions and temperature
    eqs = [
        0.0 ~ κ * ∇²( YF(x,y) ) - U[1] * Dx( YF(x,y) ) - U[2] * Dy( YF(x,y) ) +
                sourceTerm(1)

        0.0 ~ κ * ∇²( YO(x,y) ) - U[1] * Dx( YO(x,y) ) - U[2] * Dy( YO(x,y) ) +
                sourceTerm(2)

        0.0 ~ κ * ∇²( YP(x,y) ) - U[1] * Dx( YP(x,y) ) - U[2] * Dy( YP(x,y) ) +
                sourceTerm(3)

        0.0 ~ κ * ∇²( T(x,y) )  - U[1] * Dx( T(x,y) )  - U[2] * Dy( T(x,y) ) +
                Q * sourceTerm(3)
    ]
else
    # diffusion and convection terms for the mass fractions and temperature
    eqs = [
        0.0 ~ κ * ∇²( YF(x,y) ) - U[1] * Dx( YF(x,y) ) - U[2] * Dy( YF(x,y) )
        0.0 ~ κ * ∇²( YO(x,y) ) - U[1] * Dx( YO(x,y) ) - U[2] * Dy( YO(x,y) )
        0.0 ~ κ * ∇²( YP(x,y) ) - U[1] * Dx( YP(x,y) ) - U[2] * Dy( YP(x,y) )
        0.0 ~ κ * ∇²( T(x,y) )  - U[1] * Dx( T(x,y) )  - U[2] * Dy( T(x,y) )
    ]
end

function atInlet(x,y)
    return (inlety[1] < y) * (y < inlety[2])
end

bcs = [
    # this is a steady problem, i.e. d./dt = 0, so no initial conditions (?)
    # -> not sure whether MTK will allow this to be modeled this way

    #
    # boundary conditions
    #

    # left side
    T(xmin,y) ~ atInlet(xmin,y) * 950 + (1-atInlet(xmin,y)) *  300, # K
    YF(xmin,y) ~ atInlet(xmin,y) * 0.0282, # mass fraction
    YO(xmin,y) ~ atInlet(xmin,y) * 0.2259,
    YP(xmin,y) ~ 0.0,

    # right side
    Dx( T(xmax,y) ) ~ 0.0, # K
    Dx( YF(xmax,y) ) ~ 0.0, # mass fraction
    Dx( YO(xmax,y) ) ~ 0.0,
    Dx( YP(xmax,y) ) ~ 0.0,

    # top
    Dy( T(x,ymax) ) ~ 0.0, # K
    Dy( YF(x,ymax) ) ~ 0.0, # mass fraction
    Dy( YO(x,ymax) ) ~ 0.0,
    Dy( YP(x,ymax) ) ~ 0.0,

    # bottom
    Dy( T(x,ymin) ) ~ 0.0, # K
    Dy( YF(x,ymin) ) ~ 0.0, # mass fraction
    Dy( YO(x,ymin) ) ~ 0.0,
    Dy( YP(x,ymin) ) ~ 0.0,
]

@named pdesys = PDESystem(
    eqs, # partial differential equations
    bcs, # initial/boundary conditions
    domains, # domain of the independent variables (i.e. spatial/time)
    [x,y], # independent variables
    [YF(x,y), YO(x,y), YP(x,y), T(x,y)], # dependent variables
    [A => A_default, E => E_default], # parameters
) 


## Discretize the system

# TODO explain reduced node count
# Ngrid = (72, 36)
Ngrid = (floor(72/3), floor(36/3))

# Ngrid = (4, 4)
dx = (xmax-xmin)/Ngrid[1]
dy = (ymax-ymin)/Ngrid[2]
order = 2

discretization = MOLFiniteDifference(
    [x=>dx, y=>dy], approx_order=order
)

# this creates an ODEProblem or a NonlinearProblem, depending on whether the system
# is a steady/nonsteady PDE
problem = discretize(pdesys, discretization)
