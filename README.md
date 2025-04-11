# SellmeierEqn

[![Build Status](https://github.com/stakahama/SellmeierEqn.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stakahama/SellmeierEqn.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package encodes the Sellmeier equation and parameters to estimate the refactive index of crystals. Sources:

- Ramer and Lendl, doi:[10.1002/9780470027318.a9287](https://doi.org/10.1002/9780470027318.a9287), 2006.
- Weber, *Handbook of Optical Materials*, CRC Press, 2003. Chapter 1.3.4, "Dispersion Formulas for Refractive Index".

Materials implemented so far.
```julia
keys(MATERIAL)
```

Can accept scalar or vector values for wavenumbers.
```julia
n = refidx(1e4 / 1000, :ZnSe)
```

```julia
using Plots
ν = range(4000, 400, length = 100)
n = refidx(1e4 ./ ν, :ZnSe)
plot(ν, n, xflip = true, legend = false)
```

The call above is identical to the following.
```julia
refidx(ν, MATERIAL[:ZnSe])
```
Therefore, users can provide their own parameters in the form of `MATERIAL[:ZnSe]`. Additionally, a struct of type `SellmeierFn` with a function and wavelength bounds can be provided if the equation does not conform to the conventional form.

When there are ordinary / extraordinary rays or optic axes defined for the material, the function will return a tuple of vectors.
```julia
refidx(.5, :SiO2)
refidx([.5, .6, .8], :SiO2)
```

Examples of parameter formats.
```julia
MATERIAL[:ZnSe]
MATERIAL[:SiO2]
MATERIAL[:LiB3O5]
```

Related packages:

- [Sellmeier.jl](https://github.com/jagot/Sellmeier.jl): clever implementation but form of equation is more limited
- [opticalmaterialspy](https://github.com/jtambasco/opticalmaterialspy): ...
