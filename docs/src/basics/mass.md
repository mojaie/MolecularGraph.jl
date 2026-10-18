# Mass and isotopes

This page covers fundamental operations related to molecular mass and isotopes.

- Molecular weight and exact mass
- Uncertainty
- Isotopic composition
- Simulated mass spectrum

These are based on the experimental isotopic composition data provided by NIST
([https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses)),
and are not intended for structure-based mass spectrometry analysis (e.g. prediction of
fragmentation patterns).

```@example mass
using MolecularGraph, Graphs
```


## Molecular weight and exact mass

[`standard_weight`](@ref)`(mol[, digits])` returns the standard molecular weight (relative molecular mass).

```@example mass
standard_weight(smilestomol("CCO"), 2)
```

Both [`exact_mass`](@ref)`(mol[, digits])` and [`monoiso_mass`](@ref)`(mol[, digits])` return
the exact mass of the given molecule. The difference is that `exact_mass` considers the atomic
mass specified by the isotope of each atom (e.g. `[2H]` in SMILES), whereas `monoiso_mass`
always returns the mass of the molecule consisting of the most abundant isotopes in terrestrial
sources (monoisotopic mass).

```@example mass
etohd6 = smilestomol("[2H]C([2H])([2H])C([2H])([2H])O[2H]")
```

```@example mass
standard_weight(etohd6, 2)
```

```@example mass
monoiso_mass(etohd6, 6)
```

```@example mass
exact_mass(etohd6, 6)
```

[`nominal_mass`](@ref)`(mol)` returns the monoisotopic mass rounded to an integer.

```@example mass
nominal_mass(etohd6)
```


## Uncertainty

The above functions have `_unc` versions (`standard_weight_unc`, `exact_mass_unc` and
`monoiso_mass_unc`) that return a tuple of the mass value and its uncertainty. The following
result means that the molecular weight of C2H6O is 46.069±0.005.

```@example mass
standard_weight_unc(smilestomol("CCO"))
```


## Isotopic composition

- [`isotopic_composition`](@ref)`(atomsymbol, num)` returns a vector of tuples of the total mass
  of the given number of atoms and its abundance for each isotopic combination.
- Isotopic combinations that have lower abundance than the `threshold` parameter are filtered
  out (default 0.001 = 0.1%).
- `isotopic_composition(mol)` returns the isotopic composition of the molecule.

```@example mass
isotopic_composition(:C, 100; threshold=0.01)
```

```@example mass
isotopic_composition(:C, 1000; threshold=0.01)
```

```@example mass
isotopic_composition(smilestomol("CCO"))
```


## Simulated mass spectrum

[`massspec_peaks`](@ref)`(mol)` returns the masses of the isotopic peaks and their relative
intensities (base peak = 100).

```@example mass
mol = smilestomol("c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]")
```

```@example mass
massspec_peaks(mol; threshold=0.01)
```

[`simulate_massspec`](@ref)`(mol; threshold, resolution, rate)` returns a matrix of a simulated
mass spectrum (each row is a data point of mass and intensity), which is calculated from the
isotopic composition by superposing Gaussian peaks with the given `resolution`, sampled every
`rate` mass units.

```@example mass
data = simulate_massspec(mol)
size(data)
```

The spectrum can be plotted with any plotting library. The following example uses Plots.jl
(not executed here).

```julia
using Plots

data = simulate_massspec(mol)
plot(data[:, 1], data[:, 2], size=(320, 240), leg=false, xlabel="Mass", ylabel="Intensity")
```
