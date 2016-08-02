# WiltonInts84

[![Build Status](https://travis-ci.org/krcools/WiltonInts84.jl.svg?branch=master)](https://travis-ci.org/krcools/WiltonInts84.jl)
[![Coverage Status](https://coveralls.io/repos/github/krcools/WiltonInts84.jl/badge.svg?branch=master)](https://coveralls.io/github/krcools/WiltonInts84.jl?branch=master)

A package to compute integrals of powers of $R$ over flat triangles.

```math
\int_T |x-y|^n dy
```

Meant as the basis for singularity extraction type strategies for the computation of near singular integrals as encountered in the acoustic and electromagnetic boundary element method.

The methods here are generalisations of those described in:

[[1]](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1143304&tag=1) D. Wilton, S. Rao, A. Glisson, D. Schaubert, O. Al-Bundak, and C. Butler, “Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains,” IEEE Transactions on Antennas and Propagation, vol. 32, no. 3, pp. 276–281, Mar. 1984.

For details on how the techniques described here where generalised and the implementation, see the Jupyter [notebook](http://nbviewer.jupyter.org/github/krcools/WiltonInts84.jl/blob/master/docs/notebooks/Wilton_integrals_up_to_arbitrary_degree.ipynb) in this repo.

## Usage Example

```julia
using FixedSizeArrays

v1 = Vec(0.0, 0.0, 0.0)
v2 = Vec(1.0, 0.0, 0.0)
v3 = Vec(0.0, 1.0, 0.0)
n = normalize(cross(v1-v3,v2-v3))
x = (v1 + v2 + v3)/3 + 20n

I = WI.wiltonints(v1,v2,v3,x,Val{7})
```

*note*: `I[1]` will contain the integral of ``R^{-3}``, `I[2]` the integral of ``R^{-1}``, `I[i]` the integral of ``R^{i-3}`` for `i` larger or equal than `3`. The integral of ``R^{-2}`` is not computed. This is a special case that can only be expressed in terms of rather exotic special functions. Fortunately in boundary element methods this case never is required in the computation of interaction elements (not a coincidence I'm sure!). Users need to be aware however when indexing into the result array.


*note*: In the example above points as provide by `FixedSizeArrays` are used. The code itself however does not rely upon this and any object complying to the vaguely defined notion of point semantics should work.
