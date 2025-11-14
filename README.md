# Monotonicity of the first Dirichlet eigenvalue of a regular polygon

This repository contains the code for the computer assisted parts of
the proof in the paper [Monotonicity of the first Dirichlet eigenvalue
of a regular polygon](TODO).

The results of the paper are presented in different
[Pluto.jl](https://plutojl.org/) notebooks, found in the
[`proof`](proof) directory. These notebooks are responsible for
generating all the numbers and figures that appear in the paper. It is
possible to view the results of notebooks without running any code by
opening the corresponding html-files found in the [`proof`](proof)
directory, they can be opened in any browser such as Firefox. The
notebooks contains proofs for the following lemmas and propositions in
the paper:
- Lemma 2.6 (`lemma_2_6.jl`)
- Lemma 2.7 (`lemma_2_7.jl`)
- Lemma 2.10 (`lemma_2_10.jl`)
- Lemma 2.13 (`lemma_2_13.jl`)
- Lemma 2.15 (`lemma_2_15.jl`)
- Lemma 2.16 (`lemma_2_16.jl`)
- Lemma 2.21 (`lemma_2_21.jl`)
- Lemma 2.22 (`lemma_2_22.jl`)
- Proposition 3.1 (`proposition_3_1.jl`)
- Lemma C.1 (`lemma_C_1.jl`)
- Lemma C.3(`lemma_C_3.jl`)

For some of the shorted proofs (for example Lemma 2.7) most of the
implementation is given directly in the notebook. For other ones most
of the implementation is given in the [`src`](src) directory.

## Reproducing the proof
The proofs were generated with Julia version 1.11.7. This repository
contains the same `Manifest.toml` file as was used when running the
proofs, this allows installing exactly the same versions of the Julia
packages.

You can see if the package is working properly by running the tests.
You can do this by starting Julia from this directory and running

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

This will likely take some time the first time you run it since it has
to download and compile all the packages. If the tests run
successfully then you should be good to go!

To run any of the notebooks you first need to start Pluto, starting Julia
from this directory you can run
``` julia
using Pkg
Pkg.activate(".")
using Pluto
Pluto.run()
```

which should open a Pluto tab in your browser. Now you can open the
notebooks inside the `proof` directory through this and it should run
the proof.

## Notes about implementation
The implementation is split into two parts, one for large $N$
corresponding to Section 2 in the paper and one for small $N$
corresponding to Section 3 in the paper. The code for the large $N$
parts is found in [`src/large_N`](src/large_N) and for the small $N$
in [`src/small_N`](src/small_N).

### Large $N$


### Small $N$
The code for this part computes enclosures of the first eigenvalue of
regular polygons with $N$ sides using the [Method of Particular
Solutions](https://doi.org/10.1137/S0036144503437336). The
implementation consists of the following files:

- [`src/small_N/Polar.jl`](src/small_N/Polar.jl): Contains code for
  representing points in polar coordinates.
- [`src/small_N/RegularPolygon.jl`](src/small_N/RegularPolygon.jl):
  Contains code for representing and computing with regular polygons.
  Such as code for computing vertices and angles, sampling points from
  the boundary and the interior and converting between coordinate
  systems adapted for expansions around the vertices and the center.
- [`src/small_N/MPS/Eigenfunction.jl`](src/small_N/MPS/Eigenfunction.jl):
  Contains code for representing the approximate eigenfunction of a
  regular polygon. The main type is `Eigenfunction`, which internally
  consists of a `VertexExpansion` representing the expansions at the
  vertices of the polygon and an `InteriorExpansion` representing the
  expansion at the center of the polygon.
- [`src/small_N/MPS/sigma.jl`](src/small_N/MPS/sigma.jl): Contains
  code for computing the $\sigma(\lambda)$ function that is in the MPS.
- [`src/small_N/MPS/mps.jl`](src/small_N/MPS/mps.jl): Contains code
  for applying the MPS, in practice code for minimizing
  $\sigma(\lambda)$.
- [`src/small_N/MPS/enclosing/maximum_boundary.jl`](src/small_N/MPS/enclosing/maximum_boundary.jl):
  Contains code for enclosing the maximum of an approximate
  eigenfunction on the boundary of its domain.
- [`src/small_N/MPS/enclosing/norm.jl`](src/small_N/MPS/enclosing/norm.jl):
  Contains code for lower bounding the norm of an approximate
  eigenfunction on its domain.
- [`src/small_N/MPS/enclosing/eigenvalue.jl`](src/small_N/MPS/enclosing/eigenvalue.jl):
  Contains code for computing an enclosure of an eigenvalue given an
  approximate eigenvalue and eigenfunction.
- [`src/small_N/precomputed_eigenvalues.jl`](src/small_N/precomputed_eigenvalues.jl):
  Contains utility functions for loaded precomputed approximate
  eigenvalues taken from [the
  repository](https://github.com/David-Berghaus/master-thesis-data)
  associated to [Computation of Laplacian eigenvalues of
  two-dimensional shapes with dihedral
  symmetry](http://dx.doi.org/10.1007/s10444-024-10138-3).
