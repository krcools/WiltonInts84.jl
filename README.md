# WiltonInts84

[![Build Status](https://travis-ci.org/krcools/WiltonInts84.jl.svg?branch=master)](https://travis-ci.org/krcools/WiltonInts84.jl)
[![Coverage Status](https://coveralls.io/repos/github/krcools/WiltonInts84.jl/badge.svg?branch=master)](https://coveralls.io/github/krcools/WiltonInts84.jl?branch=master)

A package to compute

```math
\int_T |x-y|^n dy
```

where ``T`` is a planar triangle.

Meant as the basis for singularity extraction type strategies for the computation of near singular integrals as encountered in the acoustic and electromagnetic boundary element method.

The methods here are generalisations of those described in:

[1] D. Wilton, S. Rao, A. Glisson, D. Schaubert, O. Al-Bundak, and C. Butler, “Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains,” IEEE Transactions on Antennas and Propagation, vol. 32, no. 3, pp. 276–281, Mar. 1984.
