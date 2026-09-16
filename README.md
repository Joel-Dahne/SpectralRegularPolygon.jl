# Monotonicity of the first Dirichlet eigenvalue of regular polygons

This repository contains the code for the computer-assisted parts of
the proof in the paper [Monotonicity of the first Dirichlet eigenvalue
of regular polygons](https://doi.org/10.48550/arXiv.2601.16285).

The results of the paper are presented in different
[Pluto.jl](https://plutojl.org/) notebooks, found in the
[`proofs`](proofs) directory. These notebooks are responsible for
generating all the numbers and figures that appear in the paper. It is
possible to view the results of notebooks without running any code by
opening the corresponding html-files found in the [`proofs`](proofs)
directory, they can be opened in any browser such as Firefox. The
notebooks contain proofs for the following lemmas and propositions in
the paper:
- [Lemma 2.6](proofs/lemma_2_6.jl)
- [Lemma 2.7](proofs/lemma_2_7.jl)
- [Lemma 2.10](proofs/lemma_2_10.jl)
- [Lemma 2.13](proofs/lemma_2_13.jl)
- [Lemma 2.15](proofs/lemma_2_15.jl)
- [Lemma 2.16](proofs/lemma_2_16.jl)
- [Corollary 2.18](proofs/corollary_2_18.jl)
- [Lemma 2.21](proofs/lemma_2_21.jl)
- [Lemma 2.22](proofs/lemma_2_22.jl)
- [Proposition 3.1](proofs/proposition_3_1.jl)
- [Lemma C.1](proofs/lemma_C_1.jl)
- [Lemma C.3](proofs/lemma_C_3.jl)

For some of the shorter proofs, for example Lemma 2.7, most of the
implementation is given directly in the notebook. For other ones, for
example Lemma 2.10, most of the implementation is given in the
[`src`](src) directory.

For most notebooks the only output is to the notebook itself. The only
exception is [Proposition 3.1](proofs/proposition_3_1.jl) which
generates the figures that appear in the paper and also stores the
computed eigenvalues to
[`proofs/data/eigenvalues.csv`](proofs/data/eigenvalues.csv). The
stored eigenvalues can be loaded with

``` julia
using Arblib, CSV
data = CSV.File("proofs/data/eigenvalues.csv", types = [Int, String, String])
Ns = data.N
λs = Arblib.load_string.(Arb, data.λ_N_dump)
```

## Reproducing the proof
The proofs were generated with Julia version 1.11.8. This repository
contains the same `Manifest-v1.11.toml` file as was used when running
the proofs, this allows installing exactly the same versions of the
Julia packages.

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

To run any of the notebooks, you first need to start Pluto. Start
Julia from this directory and run:

``` julia
using Pkg
Pkg.activate(".")
using Pluto
Pluto.run()
```

which should open a Pluto tab in your browser. Now you can open the
notebooks inside the `proofs` directory through this and it should
allow you to run the proof.

## Notes about implementation
The implementation is split into two parts, one for large $N$,
corresponding to Section 2 in the paper, and one for small $N$,
corresponding to Section 3 in the paper. The code for the large $N$
parts is found in [`src/large_N`](src/large_N) and for the small $N$
in [`src/small_N`](src/small_N).

Some good information to have:
- The code and the notebooks with the proofs are written to be as
  readable as possible, but in general they assume that the reader is
  familiar with the associated paper.
- Many of the implemented functions contain documentation that
  explains what they do. These occasionally refer to equations or
  lemmas in the paper. For these references we use the notation "REF",
  e.g. "Equation REF(18)" refers to Equation 18 in the paper.
- Many of the functions have associated tests in [`test`](test) that
  serve to increase the confidence in the implementation. All of these
  tests can be run with `Pkg.test()` as discussed above.

### Large $N$
The code for this part handles computing estimates related to the
approximate eigenfunction used in Section 2. The implementation
consists of the following files:

- [`src/large_N/constants.jl`](src/large_N/constants.jl): Contains all
  of the explicit constants in the paper which we want to prove are
  bounding the different functions.
- [`src/large_N/special_functions.jl`](src/large_N/special_functions.jl):
  Contains implementations of some common special functions, such as
  the hypergeometric ${}_{2}F_{1}$ function. They are mostly direct
  wrappers of corresponding functions in FLINT.
- [`src/large_N/basic_functions.jl`](src/large_N/basic_functions.jl):
  Contains implementation of some basic functions, primarily related
  to handling powers and logarithms of arguments overlapping zero. For
  example it contains the function `logabspow(x, m, y)` for computing
  $\log^m(|x|)|x|^y$ in a way that works when the interval for $x$
  overlaps zero.
- [`src/large_N/basic_integrals.jl`](src/large_N/basic_integrals.jl):
  Contains implementation of integrals from Lemma REF(A.4) in the
  paper.
- [`src/large_N/fx_div_x.jl`](src/large_N/fx_div_x.jl): Contains code
  for computing enclosures of functions around removable singularities.
- [`src/large_N/TaylorModel.jl`](src/large_N/TaylorModel.jl): Contains
  an implementation of Taylor models that are discussed in Appendix
  REF(B.4) in the paper.
- [`src/large_N/polylog.jl`](src/large_N/polylog.jl): Contains
  implementation of different versions of polylogarithms, including
  standard polylogarithms $\mathrm{Li}_s$, Nielsen generalized
  polylogarithms $S_n$ and certain multiple polylogarithms appearing
  in the paper. Most details are discussed in Appendix REF(B.3) in the
  paper.
- [`src/large_N/section_2.jl`](src/large_N/section_2.jl): Contains
  implementations of various functions that appear in Section 2 of the
  paper.
- [`src/large_N/V.jl`](src/large_N/V.jl): Contains implementations
  related to evaluating and bounding the functions $V_l$ from Equation
  REF(13) in the paper.
- [`src/large_N/lemma_2_10.jl`](src/large_N/lemma_2_10.jl): Contains
  implementation related to computing the integrals that appear in
  Lemma REF(2.10).

In general the code follows a similar notation as in the paper and it
should hopefully be relatively straightforward to understand the
correspondence between the code and the paper. One important place
where the code differs is in how the parameter $N$ is handled. In the
paper most of the expressions and expansions are written in terms of
$N$, typically with terms involving powers of 1/N. For example we have

$$\lambda_{app} = \lambda \left(1 + \frac{4\zeta(3)}{N^3} + \frac{(12 - 2\lambda)\zeta(5)}{N^5}\right),$$

with $N$ taking values in $[N_0, \infty)$. On the computer it is in
general easier to handle intervals where both endpoints are finite.
For that reason we in the code typically write the expressions in
terms of $N^{-1}$ (which we in the code call `inv_N`) instead of $N$.
So $\lambda_{app}$ is then instead written as

$$\lambda_{app} = \lambda \left(1 + 4\zeta(3)(N^{-1})^3 + (12 - 2\lambda)\zeta(5)(N^{-1})^5\right),$$

with $N^{-1}$ taking values in $[0, N_0^{-1}]$. In code this would be
written as

``` julia
λ * (1 + 4zeta(3) * inv_N^3 + (12 - 2λ) * zeta(5) * inv_N^5)
```

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
  vertices of the polygon (see Equation REF(30) in the paper) and an
  `InteriorExpansion` representing the expansion at the center of the
  polygon (see Equation REF(31) in the paper).
- [`src/small_N/MPS/sigma.jl`](src/small_N/MPS/sigma.jl): Contains
  code for computing the $\sigma(\lambda)$ function that is used in
  the MPS.
- [`src/small_N/MPS/mps.jl`](src/small_N/MPS/mps.jl): Contains code
  for applying the MPS, in practice code for minimizing
  $\sigma(\lambda)$.
- [`src/small_N/MPS/enclosing/maximum_boundary.jl`](src/small_N/MPS/enclosing/maximum_boundary.jl):
  Contains code for enclosing the maximum of an approximate
  eigenfunction on the boundary of its domain.
- [`src/small_N/MPS/enclosing/norm.jl`](src/small_N/MPS/enclosing/norm.jl):
  Contains code for lower-bounding the norm of an approximate
  eigenfunction on its domain.
- [`src/small_N/MPS/enclosing/eigenvalue.jl`](src/small_N/MPS/enclosing/eigenvalue.jl):
  Contains code for computing an enclosure of an eigenvalue given an
  approximate eigenvalue and eigenfunction. It is based on Lemma
  REF(2.2) in the paper.
- [`src/small_N/precomputed_eigenvalues.jl`](src/small_N/precomputed_eigenvalues.jl):
  Contains utility functions for loading precomputed approximate
  eigenvalues taken from [the
  repository](https://github.com/David-Berghaus/master-thesis-data)
  associated to [Computation of Laplacian eigenvalues of
  two-dimensional shapes with dihedral
  symmetry](http://dx.doi.org/10.1007/s10444-024-10138-3).
