module Intensity

using Distributions, Random

export InhomogeneousExponential

EPSILON_DIGIT = 10

AVR = AbstractVector{<:Real}
AV = AbstractVector

abstract type AbstractTermStructure end

abstract type AbstractInterExtrapolation end
abstract type PWL <: AbstractInterExtrapolation end
abstract type Constant <: AbstractInterExtrapolation end

"""
    TermStructure <: AbstractTermStructure

basic struct to store scalar quantities possessing a nontrivial term structure.
    Terms are stored as floating point numbers.
    An interpolation scheme has to be specified

    Fields:
    - terms::AbstractVector{<:Real}
    - values::AbstractVector{<:Real}
    - interpolation::Type{<:AbstractInterExtrapolation}

    _N.B._ : terms have to be chronologically ordered
"""
struct TermStructure <: AbstractTermStructure
    terms::AVR
    values::AVR
    interpolation::Type{<:AbstractInterExtrapolation}

    function TermStructure(terms::AVR, values::AVR, interpolation::Type{<:AbstractInterExtrapolation})
        @assert sort(terms) == terms "Terms should be time-ordered"
        return new(terms, values, interpolation)
    end
end


function bracket_term(terms::AVR, target::Float64)
    lowidxs = findall(terms .<= target)
    return isempty(lowidxs) ? 0 : maximum(lowidxs)
end



"""
interextrapolate
    arguments
        - ::Type{<:AbstractInterExtrapolation},
        - terms::AV,
        - qtys::AV,
        - term::Float64,

    given the terms and the quantities corresponding to each of them, determine
    the interpolated/extrapolated value according to the desired scheme.

    _Notes_
    - Unless otherwise specified, the initial time is 0.0.
    - For PWL, if the first time peg (0.0) is omitted, the initial value is assumed to
    be zero
    - For Constant, the interpolation is continuous on the left (i.e. values
    are assumed to be on (a_i, a_{i+1}])
"""
function interextrapolate(
    ::Type{<:AbstractInterExtrapolation},
    terms::AV,
    qtys::AV,
    term::Float64)
    return error("interextrapolate not implemented.")
end


function interextrapolate(::Type{<:PWL}, terms::AV, qtys::AV, term::Float64)
    lowidx = bracket_term(terms, term)
    if lowidx == 0
        t_low = 0.0
        t_high = terms[1]
        qty_low = 0
        qty_high = qtys[1]
    else
        t_low = terms[lowidx]
        qty_low = qtys[lowidx]
        if lowidx < length(terms)
            t_high = terms[lowidx + 1]
            qty_high = qtys[lowidx + 1]
        else
            t_high = NaN
            qty_high = qtys[end]
        end
    end
    delta = qty_high - qty_low
    if delta == 0.0
        return qty_low
    else
        slope = delta / (t_high - t_low)
        dt = term - t_low
        return delta / (t_high - t_low) * (term - t_low) + qty_low
    end
end


function interextrapolate(::Type{Constant}, terms::AV, qtys::AV, term::Float64)
    highidxs = findall(terms .>= term)
    return isempty(highidxs) ? qtys[end] : qtys[minimum(highidxs)]
end


function interextrapolate(ts::AbstractTermStructure, term::Float64)
    @assert term >= 0
    return interextrapolate(ts.interpolation, ts.terms, ts.values, term)
end


function interextrapolate(ts::AbstractTermStructure, terms::AVR)
    return map(x -> interextrapolate(ts.interpolation, ts.terms, ts.values, x), terms)
end

"""
    SurvivalProbability

    A simple TermStructure enforcing the semantics of survival probabilities,
    without interpolation
"""
struct SurvivalProbability <: AbstractTermStructure
    terms::AVR
    probabilities::AVR

    function SurvivalProbability(terms::AVR, probs::AVR)
        sort(terms) != terms && throw(DomainError("Terms must be sorted"))
        length(terms) != length(probs) && throw(DomainError("Terms and probabilities have different length"))
        any(probs[2:end] .> probs[1:end-1]) && throw(DomainError("Survival probabilities must be decreasing"))
        any(probs .< 0.0) && throw(DomainError("Survival probabilities must be positive"))
        any(probs .> 1.0) && throw(DomainError("Survival probabilities must be not greater than one"))
        return new(terms, probs)
    end
end


"""
    DeterministicIntensity

    implement deterministic intensity semantics, with interpolation/extrapolation
    method specified
"""
mutable struct DeterministicIntensity <: AbstractTermStructure
    terms::AVR
    values::AVR
    interpolation::Type{<:AbstractInterExtrapolation}

    function DeterministicIntensity(
        terms::AVR,
        intensities::AVR,
        interpolation::Type{<:AbstractInterExtrapolation},
    )
        any(intensities .< 0) && throw(DomainError("All intensities have to be positive"))
        length(terms) != length(intensities) && throw(DomainError("Terms and intensities have different lengths"))
        return new(terms, intensities, interpolation)
    end
end


DeterministicIntensity(terms, intensities) = DeterministicIntensity(terms, intensities, PWL)


"""
    get_rolling_probabilities

    for each term of the DeterministicIntensity, compute the survival probability
"""
function get_rolling_probabilities(di::DeterministicIntensity)
    p = 1.0
    terms = vcat([0.0], di.terms)
    deltas = terms[2:end] - terms[1:end-1]
    intensities = interextrapolate(di, terms)
    if di.interpolation === Constant
        areas = deltas .* intensities[2:end]
    elseif di.interpolation === PWL
        areas = deltas .* (intensities[2:end] + intensities[1:end-1]) * 0.5
    end
    probs = zeros(size(di.terms))
    for (i, area) in enumerate(areas)
        p *= exp(-area)
        probs[i] = p
    end
    return probs
end



"""
    get_survival_probability

    args
    - di::DeterministicIntensity
    - term::Float64
"""
function get_survival_probability(di::DeterministicIntensity, term::Float64)
    lowidx = bracket_term(di.terms, term)
    probs = get_rolling_probabilities(di)
    if lowidx == 0
        base = 0.0
        base_int = 0.0
        p = 1.0
    else
        base = di.terms[lowidx]
        base_int = di.values[lowidx]
        p = probs[lowidx]
    end
    high_int = interextrapolate(di, term)
    delta = term - base
    height = di.interpolation == Constant ? high_int : (base_int + high_int) * 0.5
    area = delta * height
    return p * exp(-area)
end

function get_survival_probability(di::DeterministicIntensity, terms::AbstractArray{Float64})
    return map(x -> get_survival_probability(di, x), terms)
end


"""
    calibrate_intensity(sp::SurvivalProbability, interpolation::Type{<:AbstractInterExtrapolation})

    args
    - sp::SurvivalProbability,
    - interpolation
"""
function calibrate_intensity(sp::SurvivalProbability, interpolation::Type{<:AbstractInterExtrapolation})
    terms = sp.terms
    probs = sp.probabilities
    intensities = Float64[]
    new_terms = Float64[]
    if terms[1] == 0
        @assert probs[1] == 1 "survival probability at time zero must be one"
        terms = terms[2:end]
        probs = probs[2:end]
    end
    logp_prev = 0.0
    base = 0.0
    base_int = 0.0
    for i in 1:length(terms)
        term = terms[i]
        delta = term - base
        area = logp_prev - log(probs[i])
        height = area / delta
        if interpolation == PWL
            new_int = 2 * height - base_int
        else
            new_int = height
        end
        new_int = round(new_int; digits=EPSILON_DIGIT)
        @assert new_int >= 0.0 "Intensity has to be positive: $new_int, $height, $base_int, $base"
        base_int = new_int
        push!(new_terms, term)
        push!(intensities, new_int)
        base = term
        logp_prev = log(probs[i])
    end
    return DeterministicIntensity(new_terms, intensities, interpolation)
end


function make_grid(start::Float64, stop::Float64, npts::Integer; dates=nothing)
    grid = range(start, stop; length=npts)
    if dates != nothing
        grid = sort(union(grid, dates))
    end
    return grid
end



struct InhomogeneousExponential <: Sampleable{Univariate, Continuous}
    intensity::DeterministicIntensity
    grid::AVR
end


"""
    rand(rng::AbstractRNG, ie::InhomogeneousExponential)

Sample using Lewis'[1] thinning algorithm

References:
[1] Lewis, P. A. W. and G. S. Shedler 1979.
    `Simulation of nonhomogeneous Poisson processes by thinning.' Naval Res. Logistics Quart, 26:403– 413.
"""
function Base.rand(rng::AbstractRNG, ie::InhomogeneousExponential)
    t = []
    s = [0.0]
    λs = interextrapolate(ie.intensity, ie.grid)
    λmax = maximum(λs)
    while length(t) == 0
        ξ = - log(rand(rng)) / λmax
        s_new = s[end] + ξ
        push!(s, s_new)
        λnew = interextrapolate(ie.intensity, s_new)
        θ = λnew / λmax
        ζ = rand(rng)
        ζ <= θ && push!(t, s_new)
    end
    return first(t)
end



end
