# Intensity

Defines InhomogeneousExponential, a Sampleable{Univariate, Continuous}
(see [Distributions.jl](https://juliastats.org/Distributions.jl/stable/)) that
generates random positive numbers having a exponential distribution with a
nontrivial hazard function.

The hazard function is specified through an intensity, which is given on a grid,
with an interpolation scheme used to extrapolate or interpolate. Functionality
is provided to allow the construction of the intensity out of pre-specified
survival probabilities.

## Usage

```julia
using Intensity, Random, Plots
gr()

# Start from observed survival probabilities at time pegs
ts = [0.0, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00]
obs_surv_prob = [1.0, 0.99, 0.95, 0.88, 0.80, 0.75, 0.70, 0.60, 0.45]
sp = Intensity.SurvivalProbability(ts, obs_surv_prob)
# From survival probabilities, calibrate a piecewise linear intensity matching it.
dc = Intensity.calibrate_intensity(sp, Intensity.PWL)
# Create the InhomogeneousExponential
ie = InhomogeneousExponential(dc, 0:0.01:5)
# Sample and verify that the empirical distribution matches the survival probability
# computed on a grid
rng = MersenneTwister(31)
NPTS = 3000
dts = rand(rng, ie, NPTS)
surv_prob = Intensity.get_survival_probability(ie.intensity, ie.grid)
emp_surv_prob = [1 - mean(dts .< t) for t in ie.grid]
plt = plot()
plot!(plt, ie.grid, surv_prob; label="Th.")
plot!(plt, ie.grid, emp_surv_prob; label="emp.")
```
