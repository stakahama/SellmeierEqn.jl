module SellmeierEqn

export MATERIAL, Parameters, Coefficients, Equation, ParametersTuple, Bounds, refidx

"""
Coefficients for the Sellmeier equation.

n² = A + ∑ (B λ²) / (λ²+ C)

Sometimes the coefficients are written with the notation `λ₀²` instead of `C`, and the values of `λ₀` are provided.`A` == 1 unless otherwise specified. Sources:

- Ramer and Lendl, doi:[10.1002/9780470027318.a9287](https://doi.org/10.1002/9780470027318.a9287), 2006.
- Weber, *Handbook of Optical Materials*, CRC Press, 2003. Chapter 1.3.4, "Dispersion Formulas for Refractive Index".

"""
struct Coefficients
    A::Float64
    B::Vector{Float64}
    C::Vector{Float64}

    # Inner constructor
    function Coefficients(A::Float64, B::Vector{<:Real}, C::Vector{<:Real})
        length(B) == length(C) || throw(DimensionMismatch("B and C must have equal lengths"))
        new(A, Float64.(B), Float64.(C))
    end    
end

Coefficients(B::Vector{Float64}, C::Vector{Float64}) = Coefficients(1.0, B::Vector{Float64}, C::Vector{Float64})

"""
Sellmeier equation with coefficients (appropriate to use when the equations do not follow the canonical form). Sources:

- Ramer and Lendl, doi:[10.1002/9780470027318.a9287](https://doi.org/10.1002/9780470027318.a9287), 2006.
- Weber, *Handbook of Optical Materials*, CRC Press, 2003. Chapter 1.3.4, "Dispersion Formulas for Refractive Index".

"""
struct Equation{F}
    f::F
end

struct Bounds
    lower::Real
    upper::Real
end


const ParametersTuple = NamedTuple{<:Any, <:Tuple{Vararg{Union{Coefficients, Equation}}}}

"""
Parameters is a container for Coefficients, Equation, or a ParametersTuple (a NamedTuple containing the former two).
"""
struct Parameters
    params::Union{Coefficients, Equation, ParametersTuple}
    bounds::Bounds
end

"""
    isin(λ::Real, b::Bounds)

Test if λ is in the interval b.

To use in vectorized form:
```
isin.(λ, Ref(b))
```
"""
function isin(λ::Real, b::Bounds)
    mod(searchsortedlast([b.lower, b.upper], λ), 2) == 0
end

"""
Optical constants. (Need to populate this dictionary.)

- ordinary and extraordinary rays: o, e
- optic/crystallographic axes: x, y, z

Sources:

- Ramer and Lendl, doi:[10.1002/9780470027318.a9287](https://doi.org/10.1002/9780470027318.a9287), 2006.
- Weber, *Handbook of Optical Materials*, CRC Press, 2003. Chapter 1.3.4, "Dispersion Formulas for Refractive Index".
"""
const MATERIAL = Dict(
    :ZnSe => Parameters(
        Coefficients(
            [4.2980149, 0.62776557, 2.8955633],
            [0.1920630, 0.37878260, 46.994595].^2
        ),
        Bounds(0.55, 18)
    ),
    :Ge => Parameters(
        Coefficients(
            9.28156,
            [6.72880, 0.21307],
            [0.44105, 3870.1].^2
        ),
        Bounds(2, 12)
    ),
    :AgCl => Parameters(
        Coefficients(
            [2.062508, 0.9461465, 4.300785],
            [0.1039054, 0.2438691, 70.85723].^2
        ),
        Bounds(0.54, 21.0)
    ),
    :SiO2 => Parameters(
        (
            nₒ = Coefficients(
                [0.663044, 0.517852, 0.175912, 0.565380, 1.675299],
                [0.060, 0.106, 0.119, 8.844, 20.742].^2
            ),
            nₑ = Coefficients(
                [0.665721, 0.503511, 0.214792, 0.539173, 1.807613],
                [0.060, 0.106, 0.119, 8.792, 197.709].^2
            )
        ),
        Bounds(0.18, 0.72)
    ),
    :LiB3O5 => Parameters(
        (
            nx = Equation(
                λ -> 2.45768 + 0.0098877 * λ^2 / (λ^2 - 0.026095^2) - 0.013847 * λ^2,
            ),
            ny = Equation(
                λ -> 2.52500 + 0.017123 * λ^2 / (λ^2 - 0.0060517^2) - 0.0087838 * λ^2,
            ),
            nz = Equation(
                λ -> 2.58488 + 0.012737 * λ^2 / (λ^2 - 0.016293^2) - 0.016293 * λ^2,
            )
        ),
        Bounds(0.29, 1.06)
    )
)

"""

    refidx(λ::Union{Real, AbstractVector{<:Real}}, material::Symbol)
    refidx(λ::Union{Real, AbstractVector{<:Real}}, param::Parameters)
    refidx(λ::Union{Real, AbstractVector{<:Real}}, pt::ParametersTuple, b::Bounds)
    refidx(λ::AbstractVector{<:Real}, p::Union{Coefficients, Equation}, b::Bounds)
    refidx(λ::Real, p::Coefficients, b::Bounds)
    refidx(λ::Real, f::Equation, b::Bounds)

Calculate refractive index from Sellmeier's equation using pre-defined entry in material library, or provide a custom entry (as a `Parameters` struct). Canonical Sellmeier's equation:

n² = A + ∑ (B λ²) / (λ²+ C)

### Parameters
- `λ`: scalar or vector of wavelengths (in μm)
- other

### Returns
scalar or vector of refractive indices
"""
function refidx end

function refidx(λ::Union{Real, AbstractVector{<:Real}}, material::Symbol)
    param = get(MATERIAL, material) do
        throw(ArgumentError("Material $material not found"))
    end
    refidx(λ, param)
end

function refidx(λ::Union{Real, AbstractVector{<:Real}}, param::Parameters)
    refidx(λ, param.params, param.bounds)    
end


function refidx(λ::Union{Real, AbstractVector{<:Real}}, pt::ParametersTuple, b::Bounds)
    map(p -> refidx(λ, p, b), pt)
end

function refidx(λ::AbstractVector{<:Real}, p::Union{Coefficients, Equation}, b::Bounds)
    refidx.(λ, Ref(p), Ref(b))
end

function refidx(λ::Real, p::Coefficients, b::Bounds)
    if isin(λ, b) return missing end    
    (; A, B, C) = p
    n² = A + mapreduce((B, C) -> B * λ^2 / (λ^2 - C), +, B, C)
    √n²
end

function refidx(λ::Real, f::Equation, b::Bounds)
    if isin(λ, b) return missing end
    n² = f(λ)
    √n²
end

end # SellmeierEqn
